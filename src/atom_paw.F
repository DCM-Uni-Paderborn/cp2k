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
   USE atom_types,                      ONLY: atom_basis_type,&
                                              atom_potential_type,&
                                              upf_pseudo
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
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

   PUBLIC :: atom_paw_type, atom_paw_init, atom_paw_add

CONTAINS

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
