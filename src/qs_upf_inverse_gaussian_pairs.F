!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Fixed native Gaussian OCE pairs with complex Bloch density and operator adjoints.
MODULE qs_upf_inverse_gaussian_pairs
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              clist_type,&
                                              get_alist
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE gaussian_pair_type
      INTEGER :: left_atom = 0, right_atom = 0, cell(3) = 0, left_image(3) = 0, right_image(3) = 0
      INTEGER, ALLOCATABLE :: rows(:), cols(:)
      REAL(dp) :: factor = 0.0_dp
      REAL(dp), ALLOCATABLE :: left_h(:, :, :), left_s(:, :, :), right_h(:, :, :), right_s(:, :, :)
   END TYPE gaussian_pair_type

   TYPE qs_upf_inverse_gaussian_pairs_type
      PRIVATE
      INTEGER :: nao = 0, ncol = 0, nprimitive = 0, nprimitive_right = 0, center = 0
      LOGICAL :: derivatives = .FALSE., ordered = .FALSE.
      TYPE(mp_comm_type) :: comm
      TYPE(gaussian_pair_type), ALLOCATABLE :: pairs(:)
   END TYPE qs_upf_inverse_gaussian_pairs_type

   PUBLIC :: qs_upf_inverse_gaussian_pairs_type, qs_upf_inverse_gaussian_pairs_prepare
   PUBLIC :: qs_upf_inverse_gaussian_pairs_density, qs_upf_inverse_gaussian_pairs_integrate
   PUBLIC :: qs_upf_inverse_gaussian_pairs_has_derivatives

CONTAINS

! **************************************************************************************************
!> \brief Whether the owned native pair coefficients include all first geometry derivatives.
!> \param map Prepared pair map at the current geometry.
!> \return True when every retained pair has the required response coefficients.
! **************************************************************************************************
   FUNCTION qs_upf_inverse_gaussian_pairs_has_derivatives(map) RESULT(ready)
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: map
      LOGICAL                                            :: ready

      ready = map%derivatives
   END FUNCTION qs_upf_inverse_gaussian_pairs_has_derivatives

! **************************************************************************************************
!> \brief Freeze the actual distributed overlap/OCE pair pattern for one native Gaussian center.
!>        No density-dependent CPC screening or additional factorized pair products are used.
!> \param oce Borrowed current native one-center coefficients and available derivatives.
!> \param sab Borrowed native pairs. Cross-basis maps require the complete ordered cross-basis list.
!> \param kinds Borrowed current native kinds, including GAPW_1C bases.
!> \param center Selected original atom number, without an image shift.
!> \param center_kind Native kind index of that atom.
!> \param ao_first Global AO starting index of every original atom.
!> \param ao_last Global AO ending index of every original atom.
!> \param comm Communicator of the native distributed overlap pairs, not a point-row communicator.
!> \param map New owned pair map. Rebuild together with native OCE and overlap lists after geometry changes.
!> \param positions Original physical coordinates before native neighbor-list folding.
!> \param lattice Physical lattice vectors used by the native OCE lists.
!> \param right_oce Optional independent column-basis OCE owner. Selects ordered cross-basis pairs.
!> \param right_kinds Column-basis kinds, with the same atom-kind assignment as kinds.
!> \param right_ao_first Column AO starting indices, required with right_oce.
!> \param right_ao_last Column AO ending indices, required with right_oce.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_pairs_prepare(oce, sab, kinds, center, center_kind, ao_first, ao_last, comm, map, &
                                                    positions, lattice, right_oce, right_kinds, &
                                                    right_ao_first, right_ao_last)
      TYPE(oce_matrix_type), INTENT(IN), TARGET          :: oce
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(qs_kind_type), INTENT(IN), TARGET             :: kinds(:)
      INTEGER, INTENT(IN)                                :: center, center_kind, ao_first(:), &
                                                            ao_last(:)

      CLASS(mp_comm_type), INTENT(IN)                    :: comm
      TYPE(qs_upf_inverse_gaussian_pairs_type), INTENT(OUT) :: map
      REAL(dp), INTENT(IN)                               :: positions(:, :), lattice(3, 3)
      TYPE(oce_matrix_type), INTENT(IN), OPTIONAL, TARGET :: right_oce
      TYPE(qs_kind_type), INTENT(IN), OPTIONAL, TARGET   :: right_kinds(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: right_ao_first(:), right_ao_last(:)

      INTEGER                                            :: a, b, ca, cb, cell(3), count, iac, ibc, &
                                                            ikind, info, jkind, nkind, pass, ready
      INTEGER, ALLOCATABLE                               :: col_first(:), col_last(:)
      REAL(dp)                                           :: image_error, inverse_lattice(3, 3)
      TYPE(alist_type), POINTER                          :: alist_a, alist_b
      TYPE(clist_type), POINTER                          :: left, right
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(neighbor_list_iterator_p_type), POINTER       :: iterator(:)
      TYPE(oce_matrix_type), POINTER                     :: column_oce
      TYPE(qs_kind_type), POINTER                        :: column_kinds(:)

      nkind = SIZE(kinds)
      CPASSERT(center_kind >= 1 .AND. center_kind <= nkind)
      CPASSERT(center >= 1 .AND. center <= SIZE(ao_first))
      CPASSERT(SIZE(ao_first) == SIZE(ao_last))
      CPASSERT(ALL(SHAPE(positions) == [3, SIZE(ao_first)]))
      inverse_lattice = lattice
      CALL invmat(inverse_lattice, info)
      CPASSERT(info == 0)
      CPASSERT(ASSOCIATED(oce%intac))
      CALL get_qs_kind(kinds(center_kind), basis_set=basis, basis_type="GAPW_1C")
      CPASSERT(ASSOCIATED(basis))
      CALL get_paw_basis_info(basis, nsatbas=map%nprimitive)
      map%nao = MAXVAL(ao_last)
      map%ordered = PRESENT(right_oce)
      CPASSERT(map%ordered .EQV. PRESENT(right_kinds))
      CPASSERT(map%ordered .EQV. PRESENT(right_ao_first))
      CPASSERT(map%ordered .EQV. PRESENT(right_ao_last))
      column_oce => oce
      column_kinds => kinds
      IF (map%ordered) THEN
         CPASSERT(SIZE(right_kinds) == nkind)
         CPASSERT(SIZE(right_ao_first) == SIZE(ao_first))
         CPASSERT(SIZE(right_ao_last) == SIZE(ao_last))
         column_oce => right_oce
         column_kinds => right_kinds
         ALLOCATE (col_first, SOURCE=right_ao_first)
         ALLOCATE (col_last, SOURCE=right_ao_last)
      ELSE
         ALLOCATE (col_first, SOURCE=ao_first)
         ALLOCATE (col_last, SOURCE=ao_last)
      END IF
      CPASSERT(ASSOCIATED(column_oce%intac))
      CALL get_qs_kind(column_kinds(center_kind), basis_set=basis, basis_type="GAPW_1C")
      CPASSERT(ASSOCIATED(basis))
      CALL get_paw_basis_info(basis, nsatbas=map%nprimitive_right)
      map%ncol = MAXVAL(col_last)
      map%center = center
      CALL map%comm%set_handle(comm%get_handle())
      ready = 1
      DO pass = 1, 2
         count = 0
         CALL neighbor_list_iterator_create(iterator, sab)
         DO WHILE (neighbor_list_iterate(iterator) == 0)
            CALL get_iterator_info(iterator, ikind=ikind, jkind=jkind, iatom=a, jatom=b, cell=cell)
            iac = ikind + nkind*(center_kind - 1)
            ibc = jkind + nkind*(center_kind - 1)
            IF (.NOT. ASSOCIATED(oce%intac(iac)%alist)) CYCLE
            IF (.NOT. ASSOCIATED(column_oce%intac(ibc)%alist)) CYCLE
            CALL get_alist(oce%intac(iac), alist_a, a)
            CALL get_alist(column_oce%intac(ibc), alist_b, b)
            IF (.NOT. ASSOCIATED(alist_a) .OR. .NOT. ASSOCIATED(alist_b)) CYCLE
            DO ca = 1, alist_a%nclist
               left => alist_a%clist(ca)
               IF (left%catom /= center .OR. left%nsgf_cnt == 0) CYCLE
               DO cb = 1, alist_b%nclist
                  right => alist_b%clist(cb)
                  IF (right%catom /= center .OR. right%nsgf_cnt == 0) CYCLE
                  IF (ANY(cell + right%cell - left%cell /= 0)) CYCLE
                  count = count + 1
                  IF (pass == 2) THEN
                     ASSOCIATE (pair => map%pairs(count))
                        pair%left_atom = a
                        pair%right_atom = b
                        ! Both the AO Bloch phase and geometry images refer to the original
                        ! physical coordinates used by the inverse Gaussian fields.
                        pair%left_image = NINT(MATMUL(inverse_lattice, positions(:, a) - positions(:, center) + left%rac))
                        pair%right_image = NINT(MATMUL(inverse_lattice, positions(:, b) - positions(:, center) + right%rac))
                        image_error = NORM2(positions(:, a) - positions(:, center) + left%rac - &
                                            MATMUL(lattice, REAL(pair%left_image, dp)))
                        CPASSERT(image_error < 1.e-9_dp)
                        image_error = NORM2(positions(:, b) - positions(:, center) + right%rac - &
                                            MATMUL(lattice, REAL(pair%right_image, dp)))
                        CPASSERT(image_error < 1.e-9_dp)
                        pair%cell = pair%left_image - pair%right_image
                        pair%factor = 1.0_dp
                        IF (.NOT. map%ordered .AND. a == b) pair%factor = 0.5_dp
                        ALLOCATE (pair%rows(left%nsgf_cnt), pair%cols(right%nsgf_cnt))
                        pair%rows(:) = ao_first(a) + left%sgf_list(1:left%nsgf_cnt) - 1
                        pair%cols(:) = col_first(b) + right%sgf_list(1:right%nsgf_cnt) - 1
                        CPASSERT(ALL(pair%rows >= ao_first(a)))
                        CPASSERT(ALL(pair%rows <= ao_last(a)))
                        CPASSERT(ALL(pair%cols >= col_first(b)))
                        CPASSERT(ALL(pair%cols <= col_last(b)))
                        CALL copy_leg(left, a == center .AND. ALL(left%cell == 0), &
                                      map%nprimitive, pair%left_h, pair%left_s, ready)
                        CALL copy_leg(right, b == center .AND. ALL(right%cell == 0), &
                                      map%nprimitive_right, pair%right_h, pair%right_s, ready)
                     END ASSOCIATE
                  END IF
                  EXIT
               END DO
            END DO
         END DO
         CALL neighbor_list_iterator_release(iterator)
         IF (pass == 1) ALLOCATE (map%pairs(count))
      END DO
      CALL map%comm%min(ready)
      map%derivatives = ready == 1
   END SUBROUTINE qs_upf_inverse_gaussian_pairs_prepare

! **************************************************************************************************
!> \brief Native CPC products of a Hermitian density or an ordered rectangular complex matrix.
!>        AO occupations are included in matrix. The caller sums k weights once outside this map.
!> \param map ...
!> \param reduced_k Reduced Bloch vector, with native Gaussian phase exp(i*2*pi*k.T).
!> \param matrix Unweighted Hermitian AO density, or arbitrary ordered row/column coefficients.
!> \param hard Optional replicated compressed hard Gaussian coefficients, rectangular for an ordered map.
!> \param soft Matching soft coefficients, present with hard.
!> \param hard_response Optional fixed-matrix geometry/Bloch response. All response inputs are required.
!> \param soft_response ...
!> \param atom_response Cartesian response of all original atoms.
!> \param cell_response Physical lattice response.
!> \param reduced_k_response Reduced Bloch-vector response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_pairs_density(map, reduced_k, matrix, hard, soft, &
                                                    hard_response, soft_response, atom_response, &
                                                    cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: map
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: hard(:, :), soft(:, :), &
                                                            hard_response(:, :), &
                                                            soft_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp)                                        :: dphase, phase
      INTEGER                                            :: slot
      LOGICAL                                            :: derive
      REAL(dp)                                           :: dl(3), dr(3)

      CPASSERT(PRESENT(hard) .EQV. PRESENT(soft))
      CPASSERT(PRESENT(hard) .OR. PRESENT(hard_response))
      IF (PRESENT(hard)) CALL check_map(map, SHAPE(matrix), SHAPE(hard), SHAPE(soft))
      IF (.NOT. map%ordered) THEN
         IF (MAXVAL(ABS(matrix - CONJG(TRANSPOSE(matrix)))) > 1.e-10_dp*MAX(1.0_dp, MAXVAL(ABS(matrix)))) THEN
            CPABORT("Non-Hermitian Gaussian AO density")
         END IF
      END IF
      derive = PRESENT(hard_response)
      CPASSERT(derive .EQV. PRESENT(soft_response))
      CPASSERT(derive .EQV. PRESENT(atom_response))
      CPASSERT(derive .EQV. PRESENT(cell_response))
      CPASSERT(derive .EQV. PRESENT(reduced_k_response))
      IF (derive) THEN
         CPASSERT(map%derivatives)
         CALL check_map(map, SHAPE(matrix), SHAPE(hard_response), SHAPE(soft_response))
         CPASSERT(SIZE(atom_response, 1) == 3)
         CPASSERT(SIZE(atom_response, 2) >= map%center)
         hard_response = 0.0_dp
         soft_response = 0.0_dp
      END IF
      IF (PRESENT(hard)) THEN
         hard = 0.0_dp
         soft = 0.0_dp
      END IF
      DO slot = 1, SIZE(map%pairs)
         ASSOCIATE (pair => map%pairs(slot))
            phase = pair_phase(pair, reduced_k)
            IF (PRESENT(hard)) THEN
               CALL add_density(pair%left_h(:, :, 1), pair%right_h(:, :, 1), &
                                matrix(pair%rows, pair%cols), pair%factor*phase, hard, map%ordered)
               CALL add_density(pair%left_s(:, :, 1), pair%right_s(:, :, 1), &
                                matrix(pair%rows, pair%cols), pair%factor*phase, soft, map%ordered)
            END IF
            IF (derive) THEN
               CPASSERT(SIZE(atom_response, 2) >= pair%left_atom)
               CPASSERT(SIZE(atom_response, 2) >= pair%right_atom)
               dl = atom_response(:, pair%left_atom) - atom_response(:, map%center) - &
                    MATMUL(cell_response, REAL(pair%left_image, dp))
               dr = atom_response(:, pair%right_atom) - atom_response(:, map%center) - &
                    MATMUL(cell_response, REAL(pair%right_image, dp))
               dphase = CMPLX(0.0_dp, -twopi*SUM(reduced_k_response*pair%cell), dp)*phase
               CALL add_response(pair%left_h, pair%right_h, matrix(pair%rows, pair%cols), &
                                 pair%factor, phase, dphase, dl, dr, hard_response, map%ordered)
               CALL add_response(pair%left_s, pair%right_s, matrix(pair%rows, pair%cols), &
                                 pair%factor, phase, dphase, dl, dr, soft_response, map%ordered)
            END IF
         END ASSOCIATE
      END DO
      IF (PRESENT(hard)) THEN
         CALL map%comm%sum(hard)
         CALL map%comm%sum(soft)
      END IF
      IF (derive) THEN
         CALL map%comm%sum(hard_response)
         CALL map%comm%sum(soft_response)
      END IF
   END SUBROUTINE qs_upf_inverse_gaussian_pairs_density

! **************************************************************************************************
!> \brief Exact unweighted complex AO transpose of the retained Gaussian pair density map.
!> \param map ...
!> \param reduced_k ...
!> \param hard_adjoint Replicated real compressed primitive hard covector, including its radial measure.
!> \param soft_adjoint Matching soft covector, including its sign and radial measure.
!> \param adjoint New replicated Hermitian AO matrix or ordered rectangular Frobenius adjoint, without a mesh weight.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_pairs_integrate(map, reduced_k, hard_adjoint, soft_adjoint, adjoint)
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: map
      REAL(dp), INTENT(IN)                               :: reduced_k(3), hard_adjoint(:, :), &
                                                            soft_adjoint(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: adjoint(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :)
      INTEGER                                            :: i, j, slot
      REAL(dp), ALLOCATABLE                              :: bh(:, :), bs(:, :)

      CALL check_map(map, SHAPE(adjoint), SHAPE(hard_adjoint), SHAPE(soft_adjoint))
      ALLOCATE (bh(map%nprimitive, map%nprimitive_right), bs(map%nprimitive, map%nprimitive_right))
      IF (map%ordered) THEN
         bh(:, :) = hard_adjoint
         bs(:, :) = soft_adjoint
      ELSE
         bh(:, :) = 0.5_dp*(hard_adjoint + TRANSPOSE(hard_adjoint))
         bs(:, :) = 0.5_dp*(soft_adjoint + TRANSPOSE(soft_adjoint))
      END IF
      adjoint = CMPLX(0.0_dp, 0.0_dp, dp)
      DO slot = 1, SIZE(map%pairs)
         ASSOCIATE (pair => map%pairs(slot))
            ALLOCATE (block(SIZE(pair%rows), SIZE(pair%cols)))
            block(:, :) = pair%factor*CONJG(pair_phase(pair, reduced_k))*( &
                          MATMUL(MATMUL(pair%left_h(:, :, 1), bh), TRANSPOSE(pair%right_h(:, :, 1))) + &
                          MATMUL(MATMUL(pair%left_s(:, :, 1), bs), TRANSPOSE(pair%right_s(:, :, 1))))
            DO j = 1, SIZE(pair%cols)
               DO i = 1, SIZE(pair%rows)
                  adjoint(pair%rows(i), pair%cols(j)) = adjoint(pair%rows(i), pair%cols(j)) + block(i, j)
                  IF (.NOT. map%ordered) THEN
                     adjoint(pair%cols(j), pair%rows(i)) = adjoint(pair%cols(j), pair%rows(i)) + CONJG(block(i, j))
                  END IF
               END DO
            END DO
            DEALLOCATE (block)
         END ASSOCIATE
      END DO
      CALL map%comm%sum(adjoint)
   END SUBROUTINE qs_upf_inverse_gaussian_pairs_integrate

! **************************************************************************************************
!> \brief Copy the native hard home leg and soft off-center leg, with their actual derivative order.
!> \param ENTRY ...
!> \param home ...
!> \param n ...
!> \param hard ...
!> \param soft ...
!> \param ready ...
! **************************************************************************************************
   SUBROUTINE copy_leg(ENTRY, home, n, hard, soft, ready)
      TYPE(clist_type), INTENT(IN)                       :: entry
      LOGICAL, INTENT(IN)                                :: home
      INTEGER, INTENT(IN)                                :: n
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: hard(:, :, :), soft(:, :, :)
      INTEGER, INTENT(INOUT)                             :: ready

      INTEGER                                            :: nd

      CPASSERT(ASSOCIATED(entry%acint))
      CPASSERT(SIZE(entry%acint, 2) == n)
      ALLOCATE (hard(entry%nsgf_cnt, n, 4), soft(entry%nsgf_cnt, n, 4))
      hard = 0.0_dp
      soft = 0.0_dp
      nd = MIN(4, SIZE(entry%acint, 3))
      IF (nd /= 4 .AND. .NOT. home) ready = 0
      soft(:, :, 1:nd) = entry%acint(:, :, 1:nd)
      ! Native overlap derivatives act on the projector center. Convert to AO-minus-center displacement.
      IF (nd > 1) soft(:, :, 2:nd) = -entry%acint(:, :, 2:nd)
      IF (home) THEN
         CPASSERT(ASSOCIATED(entry%achint))
         hard(:, :, 1) = entry%achint(:, :, 1)
      ELSE
         hard(:, :, :) = soft
      END IF
   END SUBROUTINE copy_leg

! **************************************************************************************************
!> \brief ...
!> \param left ...
!> \param right ...
!> \param matrix ...
!> \param phase ...
!> \param density ...
!> \param ordered Keep a single ordered product instead of its Hermitian completion.
! **************************************************************************************************
   SUBROUTINE add_density(left, right, matrix, phase, density, ordered)
      REAL(dp), INTENT(IN)                               :: left(:, :), right(:, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :), phase
      REAL(dp), INTENT(INOUT)                            :: density(:, :)
      LOGICAL, INTENT(IN)                                :: ordered

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :)

      ALLOCATE (block(SIZE(density, 1), SIZE(density, 2)))
      block(:, :) = phase*MATMUL(TRANSPOSE(left), MATMUL(matrix, right))
      IF (ordered) THEN
         density = density + REAL(block, dp)
      ELSE
         density = density + REAL(block + CONJG(TRANSPOSE(block)), dp)
      END IF
   END SUBROUTINE add_density

! **************************************************************************************************
!> \brief ...
!> \param left ...
!> \param right ...
!> \param matrix ...
!> \param factor ...
!> \param phase ...
!> \param dphase ...
!> \param dl ...
!> \param dr ...
!> \param response ...
!> \param ordered Retain the rectangular product convention of the fixed map.
! **************************************************************************************************
   SUBROUTINE add_response(left, right, matrix, factor, phase, dphase, dl, dr, response, ordered)
      REAL(dp), INTENT(IN)                               :: left(:, :, :), right(:, :, :)
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :)
      REAL(dp), INTENT(IN)                               :: factor
      COMPLEX(dp), INTENT(IN)                            :: phase, dphase
      REAL(dp), INTENT(IN)                               :: dl(3), dr(3)
      REAL(dp), INTENT(INOUT)                            :: response(:, :)
      LOGICAL, INTENT(IN)                                :: ordered

      INTEGER                                            :: axis

      IF (dphase /= CMPLX(0.0_dp, 0.0_dp, dp)) THEN
         CALL add_density(left(:, :, 1), right(:, :, 1), matrix, factor*dphase, response, ordered)
      END IF
      DO axis = 1, 3
         IF (dl(axis) /= 0.0_dp) THEN
            CALL add_density(left(:, :, axis + 1), right(:, :, 1), matrix, factor*phase*dl(axis), response, ordered)
         END IF
         IF (dr(axis) /= 0.0_dp) THEN
            CALL add_density(left(:, :, 1), right(:, :, axis + 1), matrix, factor*phase*dr(axis), response, ordered)
         END IF
      END DO
   END SUBROUTINE add_response

! **************************************************************************************************
!> \brief ...
!> \param pair ...
!> \param reduced_k ...
!> \return ...
! **************************************************************************************************
   FUNCTION pair_phase(pair, reduced_k) RESULT(phase)
      TYPE(gaussian_pair_type), INTENT(IN)               :: pair
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp)                                        :: phase

      REAL(dp)                                           :: angle

      angle = -twopi*SUM(reduced_k*pair%cell)
      phase = CMPLX(COS(angle), SIN(angle), dp)
   END FUNCTION pair_phase

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param ao_shape ...
!> \param hard_shape ...
!> \param soft_shape ...
! **************************************************************************************************
   SUBROUTINE check_map(map, ao_shape, hard_shape, soft_shape)
      TYPE(qs_upf_inverse_gaussian_pairs_type), &
         INTENT(IN)                                      :: map
      INTEGER, INTENT(IN)                                :: ao_shape(2), hard_shape(2), soft_shape(2)

      CPASSERT(ALLOCATED(map%pairs))
      CPASSERT(ALL(ao_shape == [map%nao, map%ncol]))
      CPASSERT(ALL(hard_shape == [map%nprimitive, map%nprimitive_right]))
      CPASSERT(ALL(soft_shape == [map%nprimitive, map%nprimitive_right]))
   END SUBROUTINE check_map

END MODULE qs_upf_inverse_gaussian_pairs
