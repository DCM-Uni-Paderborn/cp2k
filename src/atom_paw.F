!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief PAW one-center corrections for spherical atomic Kohn-Sham calculations.
! **************************************************************************************************
MODULE atom_paw
   USE atom_types,                      ONLY: atom_basis_gridrep,&
                                              atom_basis_type,&
                                              atom_potential_type,&
                                              atom_type,&
                                              release_atom_basis,&
                                              upf_pseudo
   USE atom_upf,                        ONLY: upf_radial_weights
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: jacobi
   USE paw_onecenter,                   ONLY: paw_onecenter_hartree,&
                                              paw_onecenter_init,&
                                              paw_onecenter_type
   USE paw_xc,                          ONLY: paw_xc_evaluate,&
                                              paw_xc_init,&
                                              paw_xc_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_paw'

   TYPE atom_paw_type
      LOGICAL :: active = .FALSE.
      REAL(dp) :: core_energy = 0.0_dp
      INTEGER, ALLOCATABLE :: radial_index(:), angular_index(:)
      TYPE(paw_onecenter_type) :: onecenter
      TYPE(paw_xc_type) :: xcgrid
   END TYPE atom_paw_type

   PUBLIC :: atom_paw_type, atom_paw_init, atom_paw_add, atom_paw_guess

CONTAINS

! **************************************************************************************************
!> \brief Project occupied PAW reference waves into the atomic basis for a physical SCF start.
!>        The L2 projection uses the ordinary overlap; the resulting orbitals are orthonormalized
!>        in the PAW metric. Missing or linearly dependent reference states retain the default guess.
!> \param atom Atomic basis, overlap, original PAW data and requested occupations.
!> \param success Whether every occupied radial state received a usable reference wave.
! **************************************************************************************************
   SUBROUTINE atom_paw_guess(atom, success)
      TYPE(atom_type), POINTER                           :: atom
      LOGICAL, INTENT(OUT)                               :: success

      INTEGER                                            :: i, j, k, l, n, ns, pass
      REAL(dp)                                           :: norm, projection
      REAL(dp), ALLOCATABLE :: coefficients(:), diagonal_norm(:), eigenvalues(:), &
         eigenvectors(:, :), guess(:, :, :), overlap(:, :), rhs(:), weights(:)
      TYPE(atom_basis_type)                              :: radial

      success = .FALSE.
      IF (atom%potential%ppot_type /= upf_pseudo) RETURN
      IF (.NOT. atom%potential%upf_pot%is_paw) RETURN
      ASSOCIATE (pot => atom%potential%upf_pot, basis => atom%basis, integrals => atom%integrals)
         DO l = 0, atom%state%maxl_occ
            IF (COUNT(pot%lbeta == l .AND. pot%paw%occupations > 0.0_dp) < atom%state%maxn_occ(l)) RETURN
         END DO
         ALLOCATE (guess, MOLD=atom%orbitals%wfn)
         guess = 0.0_dp
         CALL atom_basis_gridrep(basis, radial, pot%r, pot%rab)
         ALLOCATE (weights, SOURCE=upf_radial_weights(pot)*pot%r)
         DO l = 0, atom%state%maxl_occ
            ns = atom%state%maxn_occ(l)
            IF (ns == 0) CYCLE
            n = basis%nbas(l)
            ALLOCATE (overlap(n, n), diagonal_norm(n), eigenvalues(n), eigenvectors(n, n), rhs(n), coefficients(n))
            DO i = 1, n
               diagonal_norm(i) = SQRT(integrals%ovlp(i, i, l))
            END DO
            DO j = 1, n
               overlap(:, j) = integrals%ovlp(:n, j, l)/(diagonal_norm(:)*diagonal_norm(j))
            END DO
            CALL jacobi(overlap, eigenvalues, eigenvectors)
            j = 0
            DO k = 1, pot%number_of_proj
               IF (pot%lbeta(k) /= l .OR. pot%paw%occupations(k) <= 0.0_dp) CYCLE
               j = j + 1
               IF (j > ns) EXIT
               rhs(:) = MATMUL(TRANSPOSE(radial%bf(:, :n, l)), weights*pot%paw%ps_wfc(:, k))/diagonal_norm
               coefficients = 0.0_dp
               DO i = 1, n
                  IF (eigenvalues(i) <= basis%eps_eig) CYCLE
                  coefficients(:) = coefficients + eigenvectors(:, i)*DOT_PRODUCT(eigenvectors(:, i), rhs)/eigenvalues(i)
               END DO
               guess(:n, j, l) = coefficients/diagonal_norm
               ! Reorthogonalize in the generalized metric before constructing a density.
               DO pass = 1, 2
                  DO i = 1, j - 1
                     projection = DOT_PRODUCT(guess(:n, i, l), MATMUL(integrals%metric(:n, :n, l), guess(:n, j, l)))
                     guess(:n, j, l) = guess(:n, j, l) - projection*guess(:n, i, l)
                  END DO
               END DO
               norm = DOT_PRODUCT(guess(:n, j, l), MATMUL(integrals%metric(:n, :n, l), guess(:n, j, l)))
               IF (norm <= MAX(basis%eps_eig, 1.e-12_dp)) THEN
                  CALL release_atom_basis(radial)
                  RETURN
               END IF
               guess(:n, j, l) = guess(:n, j, l)/SQRT(norm)
            END DO
            DEALLOCATE (overlap, diagonal_norm, eigenvalues, eigenvectors, rhs, coefficients)
         END DO
         CALL release_atom_basis(radial)
      END ASSOCIATE
      atom%orbitals%wfn = guess
      atom%orbitals%wfna = guess
      atom%orbitals%wfnb = guess
      success = .TRUE.
   END SUBROUTINE atom_paw_guess

! **************************************************************************************************
!> \brief Prepare immutable PAW data once per atomic SCF, leaving other potentials inactive.
!> \param potential Atomic pseudopotential, including original PAW data if present.
!> \param cache One-center data and the separate frozen-core energy constant.
! **************************************************************************************************
   SUBROUTINE atom_paw_init(potential, cache)
      TYPE(atom_potential_type), INTENT(IN)              :: potential
      TYPE(atom_paw_type), INTENT(OUT)                   :: cache

      TYPE(upf_projector_type)                           :: projectors

      IF (potential%ppot_type /= upf_pseudo) RETURN
      IF (.NOT. potential%upf_pot%is_paw) RETURN
      CPASSERT(ALLOCATED(potential%upf_pot%paw))
      CALL upf_projector_init(potential%upf_pot, MAXVAL(potential%upf_pot%lbeta), projectors)
      CALL paw_onecenter_init(potential%upf_pot, projectors, cache%onecenter)
      CALL paw_xc_init(potential%upf_pot, 50, cache%xcgrid)
      ALLOCATE (cache%radial_index, SOURCE=projectors%radial_index)
      ALLOCATE (cache%angular_index, SOURCE=projectors%angular_index)
      cache%core_energy = potential%upf_pot%paw%core_energy
      cache%active = .TRUE.
   END SUBROUTINE atom_paw_init

! **************************************************************************************************
!> \brief Add the PAW Hartree/XC energies and their atomic radial-basis derivatives.
!> \param cache Prepared one-center data.
!> \param basis Atomic radial basis with the original beta overlaps.
!> \param xc_section Atomic XC section, also used by the smooth-density functional.
!> \param pmata Restricted total or unrestricted alpha density, summed over magnetic states.
!> \param hartree Smooth Hartree energy on entry; complete Hartree energy on return.
!> \param exc Smooth XC energy on entry; complete XC energy on return.
!> \param jmat Hartree matrix, incremented in place.
!> \param xcmata Restricted or alpha XC matrix, incremented in place.
!> \param pmatb Optional beta density; must be present together with xcmatb.
!> \param xcmatb Optional beta XC matrix, incremented in place.
! **************************************************************************************************
   SUBROUTINE atom_paw_add(cache, basis, xc_section, pmata, hartree, exc, jmat, xcmata, pmatb, xcmatb)
      TYPE(atom_paw_type), INTENT(IN)                    :: cache
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      TYPE(section_vals_type), POINTER                   :: xc_section
      REAL(dp), INTENT(IN)                               :: pmata(:, :, 0:)
      REAL(dp), INTENT(INOUT)                            :: hartree, exc, jmat(:, :, 0:), &
                                                            xcmata(:, :, 0:)
      REAL(dp), INTENT(IN), OPTIONAL                     :: pmatb(:, :, 0:)
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: xcmatb(:, :, 0:)

      INTEGER                                            :: np, ns
      REAL(dp)                                           :: correction, density_cutoff, &
                                                            gradient_cutoff
      REAL(dp), ALLOCATABLE                              :: density(:, :, :), vh(:, :), vxc(:, :, :)
      TYPE(section_vals_type), POINTER                   :: functional

      IF (.NOT. cache%active) RETURN
      CPASSERT(PRESENT(pmatb) .EQV. PRESENT(xcmatb))
      CPASSERT(ALLOCATED(basis%uspp))
      CPASSERT(ASSOCIATED(xc_section))
      np = cache%onecenter%nproj
      ns = 1
      IF (PRESENT(pmatb)) ns = 2
      ALLOCATE (density(np, np, ns), vh(np, np), vxc(np, np, ns))
      CALL project_density(cache, basis, pmata, density(:, :, 1))
      IF (PRESENT(pmatb)) CALL project_density(cache, basis, pmatb, density(:, :, 2))
      CALL paw_onecenter_hartree(cache%onecenter, SUM(density, DIM=3), correction, vh)
      hartree = hartree + correction
      CALL integrate_potential(cache, basis, vh, jmat)
      functional => section_vals_get_subs_vals(xc_section, 'XC_FUNCTIONAL')
      CALL section_vals_val_get(xc_section, 'DENSITY_CUTOFF', r_val=density_cutoff)
      CALL section_vals_val_get(xc_section, 'GRADIENT_CUTOFF', r_val=gradient_cutoff)
      CALL paw_xc_evaluate(cache%xcgrid, cache%onecenter, functional, density, density_cutoff, gradient_cutoff, &
                           correction, vxc)
      exc = exc + correction
      CALL integrate_potential(cache, basis, vxc(:, :, 1), xcmata)
      IF (PRESENT(xcmatb)) CALL integrate_potential(cache, basis, vxc(:, :, 2), xcmatb)
   END SUBROUTINE atom_paw_add

! **************************************************************************************************
!> \brief Expand an occupation-summed spherical density into individual magnetic projector blocks.
!> \param cache Radial and magnetic projector indices.
!> \param basis Original beta overlaps with atomic radial basis functions.
!> \param pmat Atomic density summed over the 2*l+1 magnetic states.
!> \param density Full real-harmonic projector density.
! **************************************************************************************************
   SUBROUTINE project_density(cache, basis, pmat, density)
      TYPE(atom_paw_type), INTENT(IN)                    :: cache
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      REAL(dp), INTENT(IN)                               :: pmat(:, :, 0:)
      REAL(dp), INTENT(OUT)                              :: density(:, :)

      INTEGER                                            :: i, j, l, n, p, q

      density = 0.0_dp
      DO q = 1, cache%onecenter%nproj
         j = cache%radial_index(q)
         l = basis%uspp%lbeta(j)
         n = basis%nbas(l)
         DO p = 1, cache%onecenter%nproj
            IF (cache%angular_index(p) /= cache%angular_index(q)) CYCLE
            i = cache%radial_index(p)
            density(p, q) = DOT_PRODUCT(basis%uspp%beta(:n, i), &
                                        MATMUL(pmat(:n, :n, l), basis%uspp%beta(:n, j)))/REAL(2*l + 1, dp)
         END DO
      END DO
   END SUBROUTINE project_density

! **************************************************************************************************
!> \brief Exact transpose of the spherical density projection, including magnetic-state averaging.
!> \param cache Radial and magnetic projector indices.
!> \param basis Original beta overlaps with atomic radial basis functions.
!> \param potential Full projector-space energy derivative.
!> \param matrix Atomic radial-basis operator, incremented in place.
! **************************************************************************************************
   SUBROUTINE integrate_potential(cache, basis, potential, matrix)
      TYPE(atom_paw_type), INTENT(IN)                    :: cache
      TYPE(atom_basis_type), INTENT(IN)                  :: basis
      REAL(dp), INTENT(IN)                               :: potential(:, :)
      REAL(dp), INTENT(INOUT)                            :: matrix(:, :, 0:)

      INTEGER                                            :: i, j, k, l, n, p, q
      REAL(dp)                                           :: coefficient

      DO q = 1, cache%onecenter%nproj
         j = cache%radial_index(q)
         l = basis%uspp%lbeta(j)
         n = basis%nbas(l)
         DO p = 1, cache%onecenter%nproj
            IF (cache%angular_index(p) /= cache%angular_index(q)) CYCLE
            i = cache%radial_index(p)
            coefficient = potential(p, q)/REAL(2*l + 1, dp)
            DO k = 1, n
               matrix(:n, k, l) = matrix(:n, k, l) + coefficient*basis%uspp%beta(:n, i)*basis%uspp%beta(k, j)
            END DO
         END DO
      END DO
   END SUBROUTINE integrate_potential

END MODULE atom_paw
