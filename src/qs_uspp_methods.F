!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect native USPP geometry, density, and operators to Quickstep GPW.
! **************************************************************************************************
MODULE qs_uspp_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_get_info,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_constants,                 ONLY: do_ppl_analytic
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_unequal
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_derive,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_type,&
                                              qs_uspp_update
   USE sap_kind_types,                  ONLY: build_sap_ints,&
                                              release_sap_int,&
                                              sap_int_type
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
   USE uspp_projector_operators,        ONLY: uspp_projector_add_operator,&
                                              uspp_projector_collocate,&
                                              uspp_projector_create_onsite,&
                                              uspp_projector_integrate
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_methods'
   PUBLIC :: qs_uspp_prepare, qs_uspp_init_ao_matrix, qs_uspp_collocate, qs_uspp_integrate
   PUBLIC :: qs_uspp_core_matrix

CONTAINS

! **************************************************************************************************
!> \brief Prepare USPP centers, unweighted projector overlaps, and the generalized metric.
!>        Called after rebuilding the Gaussian overlap. The NC path allocates nothing.
!> \param ks_env KS environment with current kinds, particles and ABBA projector neighbor list.
!> \param gaussian_overlap Current Gaussian overlap; never modified or aliased by the metric.
! **************************************************************************************************
   SUBROUTINE qs_uspp_prepare(ks_env, gaussian_overlap)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(dbcsr_type), INTENT(IN)                       :: gaussian_overlap

      INTEGER                                            :: handle
      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CALL qs_uspp_init(uspp, qs_kind_set, particle_set)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_prepare', handle)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints, sap_ppnl=sap_ppnl, dft_control=dft_control)
      CPASSERT(.NOT. do_kpoints)
      CPASSERT(dft_control%qs_control%gpw)
      CPASSERT(.NOT. dft_control%qs_control%gapw .AND. .NOT. dft_control%qs_control%gapw_xc)
      CPASSERT(dft_control%qs_control%do_ppl_method == do_ppl_analytic)
      CPASSERT(.NOT. dft_control%use_kinetic_energy_density)
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=0)
      CALL qs_uspp_update(uspp, gaussian_overlap, sap_int, particle_set)
      CALL release_sap_int(sap_int)
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_prepare

! **************************************************************************************************
!> \brief Initialize a matrix on the union of its ordinary and USPP projector-mediated blocks.
!>        Called after ordinary neighbor-list allocation, before assigning matrix values.
!> \param ks_env KS environment with current USPP geometry data.
!> \param matrix Newly allocated AO matrix, zeroed after expanding its block pattern.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_ao_matrix(ks_env, matrix)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CPASSERT(ASSOCIATED(uspp%ao_pattern))
      CALL dbcsr_add(matrix, uspp%ao_pattern, 1.0_dp, 1.0_dp)
      CALL dbcsr_set(matrix, 0.0_dp)
   END SUBROUTINE qs_uspp_init_ao_matrix

! **************************************************************************************************
!> \brief Add USPP augmentation (or its spatial derivative) to a Quickstep density grid.
!>        The reciprocal and real representations are both updated. Ordinary kinds do nothing.
!> \param ks_env KS environment with the current projector overlaps.
!> \param density AO density for one spin, including occupations.
!> \param rho_g Existing reciprocal density, normalized as integral/cell volume.
!> \param rho_r Corresponding real density, overwritten from the updated reciprocal grid.
!> \param total_rho Optional signed electronic charge, decremented by the augmentation charge.
!> \param derivative Optional Cartesian derivative orders; total_rho must then be absent.
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate(ks_env, density, rho_g, rho_r, total_rho, derivative)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(IN)                       :: density
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: charge
      TYPE(pw_c1d_gs_type)                               :: augmentation_g
      TYPE(pw_r3d_rs_type)                               :: augmentation_r
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_collocate', handle)
      CPASSERT(ASSOCIATED(rho_g%pw_grid, rho_r%pw_grid))
      CALL check_grid(uspp, rho_g%pw_grid)
      CALL augmentation_g%create(rho_g%pw_grid)
      augmentation_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL uspp_projector_collocate(uspp%overlaps, density, uspp%augmentation, uspp%center_kinds, &
                                    uspp%position, rho_g%pw_grid%g, rho_g%pw_grid%vol, augmentation_g%array, charge)
      ! On even FULLSPACE grids, project the Nyquist coefficients onto a real grid function.
      ! This same real projection is implicit in the adjoint's FFT of a real potential.
      IF (rho_g%pw_grid%grid_span == FULLSPACE .AND. ANY(MODULO(rho_g%pw_grid%npts, 2) == 0)) THEN
         CALL augmentation_r%create(rho_g%pw_grid)
         CALL pw_transfer(augmentation_g, augmentation_r)
         CALL pw_transfer(augmentation_r, augmentation_g)
         CALL augmentation_r%release()
      END IF
      IF (PRESENT(derivative)) THEN
         CPASSERT(.NOT. PRESENT(total_rho))
         CALL pw_derive(augmentation_g, derivative)
      ELSE IF (PRESENT(total_rho)) THEN
         total_rho = total_rho - charge
      END IF
      CALL pw_axpy(augmentation_g, rho_g)
      CALL pw_transfer(rho_g, rho_r)
      CALL augmentation_g%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_collocate

! **************************************************************************************************
!> \brief Integrate a Quickstep real local potential against USPP augmentation.
!>        This is the discrete adjoint of qs_uspp_collocate, including FFT normalization.
!> \param ks_env KS environment with the current projector overlaps.
!> \param potential Real potential; values include the volume element when weighted is true.
!> \param hamiltonian AO Hamiltonian for one spin, incremented in place.
!> \param weighted Whether potential already includes the real grid's volume element.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate(ks_env, potential, hamiltonian, weighted)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL integrate_grid(uspp, potential, hamiltonian, weighted)
   END SUBROUTINE qs_uspp_integrate

! **************************************************************************************************
!> \brief Real-grid adjoint shared by the ionic and density-dependent local potentials.
!> \param uspp Current projector state.
!> \param potential Real local potential.
!> \param hamiltonian AO Hamiltonian, incremented in place.
!> \param weighted Whether the real potential includes dV.
! **************************************************************************************************
   SUBROUTINE integrate_grid(uspp, potential, hamiltonian, weighted)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE                         :: multiplicity(:)
      TYPE(pw_c1d_gs_type)                               :: potential_g

      CALL timeset('qs_uspp_integrate', handle)
      CALL check_grid(uspp, potential%pw_grid)
      CALL potential_g%create(potential%pw_grid)
      CALL pw_transfer(potential, potential_g)
      IF (weighted) CALL pw_scale(potential_g, 1.0_dp/potential%pw_grid%dvol)
      ALLOCATE (multiplicity(SIZE(potential_g%array)))
      multiplicity = 1.0_dp
      IF (potential%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (potential%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      CALL uspp_projector_integrate(uspp%overlaps, uspp%augmentation, uspp%center_kinds, &
                                    uspp%position, potential%pw_grid%g, potential%pw_grid%vol, &
                                    potential_g%array, multiplicity, hamiltonian)
      CALL potential_g%release()
      CALL timestop(handle)
   END SUBROUTINE integrate_grid

! **************************************************************************************************
!> \brief Add B (D0 + integral Q V_sr) B^T to the core Hamiltonian.
!>        V_sr includes every non-ghost UPF ion, also NC ions without augmentation.
!>        The Gaussian-core Coulomb potential enters later through the Hartree grid.
!>        Keeping this density-independent term in H_core supplies its energy via Tr(P H_core).
!> \param uspp Current geometry-dependent USPP state.
!> \param qs_kind_set Kinds with original UPF local potentials and native projector ordering.
!> \param particle_set All ions, not only augmentation centers.
!> \param grid Same FFT grid as the physical electron density.
!> \param hamiltonian Core Hamiltonian, incremented in place.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_matrix(uspp, qs_kind_set, particle_set, grid, hamiltonian)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian

      COMPLEX(KIND=dp), ALLOCATABLE                      :: ionic_g(:)
      INTEGER                                            :: atom, center, handle, ikind
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: bare
      TYPE(pw_c1d_gs_type)                               :: potential_g
      TYPE(pw_r3d_rs_type)                               :: potential_r

      CALL timeset('qs_uspp_core_matrix', handle)
      CALL check_grid(uspp, grid)
      CALL uspp_projector_create_onsite(uspp%overlaps, bare)
      DO center = 1, SIZE(uspp%center_kinds)
         ikind = uspp%center_kinds(center)
         CALL dbcsr_get_block_p(bare, center, center, block, found)
         IF (found) block(:, :) = qs_kind_set(ikind)%sgp_potential%upf_projectors%d
      END DO
      CALL uspp_projector_add_operator(uspp%overlaps, bare, hamiltonian)
      CALL dbcsr_release(bare)

      CALL potential_g%create(grid)
      CALL potential_r%create(grid)
      ALLOCATE (ionic_g(SIZE(potential_g%array)))
      potential_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      DO atom = 1, SIZE(particle_set)
         ikind = particle_set(atom)%atomic_kind%kind_number
         IF (qs_kind_set(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(qs_kind_set(ikind)%sgp_potential))
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            CPASSERT(ALLOCATED(potential%upf_data))
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            CALL upf_local_fourier(potential%upf_local, grid%g, particle_set(atom)%r, grid%vol, ionic_g)
            potential_g%array = potential_g%array + ionic_g
         END ASSOCIATE
      END DO
      ! Explicit real projection also treats aliased Nyquist pairs consistently on even grids.
      CALL pw_transfer(potential_g, potential_r)
      CALL integrate_grid(uspp, potential_r, hamiltonian, weighted=.FALSE.)
      CALL potential_g%release()
      CALL potential_r%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_core_matrix

! **************************************************************************************************
!> \brief Require every retained G to occur once over the projector communicator.
!> \param uspp Current projector state.
!> \param grid Quickstep FFT grid; a reordered communicator with the same members is allowed.
! **************************************************************************************************
   SUBROUTINE check_grid(uspp, grid)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(pw_grid_type), INTENT(IN)                     :: grid

      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CPASSERT(group%compare(grid%para%group) /= mp_comm_unequal)
      CPASSERT(grid%grid_span == FULLSPACE .OR. grid%grid_span == HALFSPACE)
      IF (grid%grid_span == HALFSPACE) THEN
         CPASSERT(ALL(MODULO(grid%npts, 2) == 1))
      END IF
   END SUBROUTINE check_grid

END MODULE qs_uspp_methods
