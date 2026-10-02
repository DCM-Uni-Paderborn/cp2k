!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_grid_unittest
   USE cell_methods,                    ONLY: cell_create
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_cp2k_xc,                   ONLY: create_xc_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE particle_types,                  ONLY: particle_type
   USE skala_gpw_features,              ONLY: skala_grid_layout_type
   USE skala_gpw_functional,            ONLY: skala_grid_eval,&
                                              skala_point_eval

   IMPLICIT NONE

   TYPE(mp_comm_type), TARGET :: comm
   TYPE(skala_grid_layout_type) :: grid
   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(particle_type), POINTER :: atoms(:) => NULL()
   TYPE(section_type), POINTER :: definition => NULL()
   TYPE(section_vals_type), POINTER :: xc => NULL()
   CHARACTER(len=default_path_length) :: model
   CHARACTER(len=32) :: model_mode
   INTEGER :: ns, partition, i, s, first_z, last_z, nr, component, step, mode
   REAL(dp) :: hmat(3, 3), energy, reference, ep, em, derivative, h, error, maximum_error, &
               step_size, routing_tolerance, derivative_tolerance
   REAL(dp), ALLOCATABLE :: rho(:, :), grad(:, :, :), tau(:, :), drho(:, :), dgrad(:, :, :), dtau(:, :), &
                            vrho(:, :), vgrad(:, :, :), vtau(:, :), ref_vrho(:, :), ref_vgrad(:, :, :), ref_vtau(:, :)

   CALL get_command_argument(1, model)
   IF (LEN_TRIM(model) == 0) ERROR STOP "Pass a model exported by create_skala_grid_test_model.py"
   IF (LEN_TRIM(model) > default_string_length) &
      ERROR STOP "Use MODEL SKALA and set GAUXC_SKALA_MODEL for a long model path"
   CALL get_command_argument(2, model_mode)
   step_size = 1.0e-3_dp
   routing_tolerance = 1.0e-11_dp
   derivative_tolerance = 1.0e-8_dp
   SELECT CASE (TRIM(model_mode))
   CASE ("", "geometry")
      CONTINUE
   CASE ("skala-float32")
      step_size = 1.0e-2_dp
      routing_tolerance = 1.0e-5_dp
      derivative_tolerance = 1.0e-4_dp
   CASE DEFAULT
      ERROR STOP "Unknown model mode (expected geometry, skala-float32 or no mode for the analytic fixture)"
   END SELECT
   CALL cp2k_init()
   comm = mp_comm_world
   CALL comm%init()
   CALL create_xc_section(definition)
   CALL section_vals_create(xc, definition)
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%MODEL", c_val=TRIM(model))
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID", l_val=.TRUE.)
   hmat = RESHAPE([6.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 5.0_dp, 0.0_dp, 0.2_dp, 0.3_dp, 7.0_dp], [3, 3])
   CALL cell_create(cell, hmat=hmat, periodic=[1, 1, 1])
   ALLOCATE (atoms(2))
   atoms(1)%r = [0.8_dp, 0.6_dp, 1.1_dp]
   atoms(2)%r = [3.5_dp, 3.0_dp, 4.2_dp]
   grid%group => comm
   grid%npts = [4, 4, 3]
   IF (model_mode == "skala-float32") grid%npts = [8, 8, 7]
   grid%bounds(1, :) = 0
   grid%bounds(2, :) = grid%npts - 1
   first_z = grid%npts(3)*comm%mepos/comm%num_pe
   last_z = grid%npts(3)*(comm%mepos + 1)/comm%num_pe - 1
   grid%bounds_local = grid%bounds
   grid%bounds_local(:, 3) = [first_z, last_z]
   nr = PRODUCT(grid%npts(1:2))*(last_z - first_z + 1)
   grid%ngpts_local = nr
   grid%dvol = cell%deth/REAL(PRODUCT(grid%npts), dp)
   DO i = 1, 3
      grid%dh(:, i) = hmat(:, i)/REAL(grid%npts(i), dp)
   END DO
   maximum_error = 0.0_dp
   DO ns = 1, 2
      ALLOCATE (rho(nr, ns), grad(nr, 3, ns), tau(nr, ns), drho(nr, ns), dgrad(nr, 3, ns), dtau(nr, ns), &
                vrho(nr, ns), vgrad(nr, 3, ns), vtau(nr, ns), &
                ref_vrho(nr, ns), ref_vgrad(nr, 3, ns), ref_vtau(nr, ns))
      DO s = 1, ns
         DO i = 1, nr
            rho(i, s) = 0.1_dp + 0.03_dp*SIN(REAL(i + 16*first_z + s, dp))
            tau(i, s) = 0.04_dp + 0.01_dp*COS(REAL(i + 16*first_z + 2*s, dp))
            grad(i, :, s) = 0.01_dp*SIN(REAL(i + 16*first_z + s + [1, 2, 3], dp))
            drho(i, s) = 0.01_dp*COS(REAL(i + 16*first_z + s, dp))
            dtau(i, s) = 0.02_dp*SIN(REAL(i + 16*first_z + 2*s, dp))
            dgrad(i, :, s) = 0.01_dp*COS(REAL(i + 16*first_z + s + [1, 2, 3], dp))
         END DO
      END DO
      IF (model_mode == "skala-float32") CALL atomic_fields()
      DO partition = 1, 2
         CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_PARTITION", i_val=partition)
         DO mode = 0, 1
            CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_CHUNKS", l_val=(mode == 1))
            CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_CHUNK_MAX_ROWS", i_val=1)
            CALL evaluate(rho, grad, tau, energy)
            IF (model_mode == "skala-float32" .AND. energy >= -1.0e-6_dp) &
               ERROR STOP "The Skala field fixture must have a nontrivial XC energy"
            IF (mode == 0) THEN
               reference = energy
               ref_vrho = vrho
               ref_vgrad = vgrad
               ref_vtau = vtau
            ELSE
               error = ABS(energy - reference)
               IF (nr > 0) error = MAX(error, MAXVAL(ABS(vrho - ref_vrho)), &
                                       MAXVAL(ABS(vgrad - ref_vgrad)), MAXVAL(ABS(vtau - ref_vtau)))
               CALL comm%max(error)
               IF (error > routing_tolerance) ERROR STOP "Atom routing changes the energy or primitive derivatives"
            END IF
            DO component = 1, 3
               SELECT CASE (component)
               CASE (1)
                  derivative = SUM(ref_vrho*drho)*grid%dvol
               CASE (2)
                  derivative = SUM(ref_vgrad*dgrad)*grid%dvol
               CASE (3)
                  derivative = SUM(ref_vtau*dtau)*grid%dvol
               END SELECT
               CALL comm%sum(derivative)
               DO step = 1, 2
                  h = step_size/REAL(step, dp)
                  SELECT CASE (component)
                  CASE (1)
                     CALL evaluate(rho + h*drho, grad, tau, ep)
                     CALL evaluate(rho - h*drho, grad, tau, em)
                  CASE (2)
                     CALL evaluate(rho, grad + h*dgrad, tau, ep)
                     CALL evaluate(rho, grad - h*dgrad, tau, em)
                  CASE (3)
                     CALL evaluate(rho, grad, tau + h*dtau, ep)
                     CALL evaluate(rho, grad, tau - h*dtau, em)
                  END SELECT
                  error = ABS((ep - em)/(2.0_dp*h) - derivative)
                  maximum_error = MAX(maximum_error, error)
                  IF (error > derivative_tolerance) THEN
                     WRITE (*, '(A,4I3,3ES18.8)') "FD failure (spin,partition,mode,component,step,AD,FD): ", &
                        ns, partition, mode, component, h, derivative, (ep - em)/(2.0_dp*h)
                     ERROR STOP "Primitive-field derivative fails finite differences"
                  END IF
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (rho, grad, tau, drho, dgrad, dtau, vrho, vgrad, vtau, ref_vrho, ref_vgrad, ref_vtau)
   END DO
   CALL arbitrary_points()
   CALL comm%max(maximum_error)
   IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Maximum field/geometry derivative error: ", maximum_error
   DEALLOCATE (atoms)
   CALL cell_release(cell)
   CALL section_vals_release(xc)
   CALL section_release(definition)
   CALL cp2k_finalize()
CONTAINS
   SUBROUTINE arbitrary_points()
      INTEGER, PARAMETER                                 :: block_sizes(2) = [160, 192], nang = 16, &
                                                            total_rows = SUM(block_sizes)

      INTEGER                                            :: a, axis, b, c, distribution, first, &
                                                            global_row, ia, ih, ir, na_local, &
                                                            nr_local, nspin, row, source, spin, &
                                                            sx, sy, sz
      INTEGER, ALLOCATABLE                               :: ids(:), sizes(:)
      REAL(dp) :: ad(7), charge, descriptor_derivative, direction(3), displacement(3), dw, e, err, &
                  expected, mean, minus, occupation, phi, pi, plus, point(3), radius, ref_e, rho_sum, &
                  spin_scale, step_h, sum_descriptor, sum_weights, w, z, mean_point(3), &
                  energy_density, point_covector(3), centre_covector(3), weight_covector
      REAL(dp), ALLOCATABLE :: centres(:, :), descriptor_weights(:), dg(:, :, :), dr(:, :), &
                               dt(:, :), energy_weights(:), g(:, :, :), gathered(:, :, :), pg(:, :, :), points(:, :), &
                               pr(:, :), pt(:, :), r(:, :), reference_grad(:, :, :), t(:, :), &
                               dx(:, :), dc(:, :), de(:), dd(:), px(:, :), pc(:, :), pe(:), pd(:)

      pi = ACOS(-1.0_dp)
      DO nspin = 1, 2
         ALLOCATE (gathered(total_rows, 13, nspin), reference_grad(total_rows, 13, nspin))
         DO distribution = 0, 1
            na_local = 0
            nr_local = 0
            DO a = 1, 2
               IF (distribution == 0 .AND. comm%mepos /= 0) CYCLE
               IF (distribution == 1 .AND. MOD(a - 1, comm%num_pe) /= comm%mepos) CYCLE
               na_local = na_local + 1
               nr_local = nr_local + block_sizes(a)
            END DO
            ALLOCATE (sizes(na_local), centres(3, na_local), points(3, nr_local), ids(nr_local), &
                      energy_weights(nr_local), descriptor_weights(nr_local), &
                      r(nr_local, nspin), g(nr_local, 3, nspin), t(nr_local, nspin), &
                      dr(nr_local, nspin), dg(nr_local, 3, nspin), dt(nr_local, nspin), &
                      pr(nr_local, nspin), pg(nr_local, 3, nspin), pt(nr_local, nspin), &
                      dx(3, nr_local), dc(3, na_local), de(nr_local), dd(nr_local), &
                      px(3, nr_local), pc(3, na_local), pe(nr_local), pd(nr_local))
            b = 0
            row = 0
            DO a = 1, 2
               IF (distribution == 0 .AND. comm%mepos /= 0) CYCLE
               IF (distribution == 1 .AND. MOD(a - 1, comm%num_pe) /= comm%mepos) CYCLE
               b = b + 1
               sizes(b) = block_sizes(a)
               centres(:, b) = atoms(a)%r
               dc(:, b) = 0.03_dp*COS(REAL(a + [1, 2, 3], dp))
               DO ir = 1, block_sizes(a)/nang
                  radius = 0.12_dp + (REAL(ir, dp) - 0.5_dp)*0.24_dp
                  DO ia = 1, nang
                     z = 1.0_dp - 2.0_dp*(REAL(ia, dp) - 0.5_dp)/REAL(nang, dp)
                     phi = REAL(ia, dp)*pi*(3.0_dp - SQRT(5.0_dp))
                     direction = [SQRT(1.0_dp - z*z)*COS(phi), SQRT(1.0_dp - z*z)*SIN(phi), z]
                     row = row + 1
                     global_row = SUM(block_sizes(1:a - 1)) + (ir - 1)*nang + ia
                     ids(row) = global_row
                     point = atoms(a)%r + radius*direction
                     points(:, row) = point
                     descriptor_weights(row) = 4.0_dp*pi*radius**2*0.24_dp/REAL(nang, dp)
                     energy_weights(row) = descriptor_weights(row)*(0.6_dp + 0.1_dp*COS(phi))
                     ! This row still participates in the nonlocal descriptor.
                     IF (ir == 1 .AND. ia == 1) energy_weights(row) = 0.0_dp
                     dx(:, row) = 0.01_dp*SIN(REAL(global_row + [1, 2, 3], dp))
                     de(row) = 0.02_dp*descriptor_weights(row)*COS(REAL(global_row, dp))
                     dd(row) = 0.03_dp*descriptor_weights(row)*SIN(REAL(global_row, dp))
                     r(row, :) = 0.0_dp
                     g(row, :, :) = 0.0_dp
                     t(row, :) = 0.0_dp
                     DO source = 1, 2
                        DO sz = -1, 1
                           DO sy = -1, 1
                              DO sx = -1, 1
                                 displacement = point - atoms(source)%r - MATMUL(hmat, REAL([sx, sy, sz], dp))
                                 charge = EXP(-2.0_dp*SQRT(SUM(displacement**2)))/pi
                                 DO spin = 1, nspin
                                    occupation = 1.0_dp
                                    IF (nspin == 2) occupation = 0.5_dp + (3 - 2*spin)*(0.1_dp + 0.05_dp*source)
                                    r(row, spin) = r(row, spin) + occupation*charge
                                    t(row, spin) = t(row, spin) + 0.5_dp*occupation*charge
                                    g(row, :, spin) = g(row, :, spin) - &
                                                      2.0_dp*occupation*charge*displacement/SQRT(SUM(displacement**2))
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                     DO spin = 1, nspin
                        dr(row, spin) = 0.03_dp*r(row, spin)*COS(REAL(global_row + spin, dp))
                        dt(row, spin) = 0.03_dp*t(row, spin)*SIN(REAL(global_row + 2*spin, dp))
                        DO axis = 1, 3
                           dg(row, axis, spin) = 0.03_dp*r(row, spin)*COS(REAL(global_row + spin + axis, dp))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, pr, pg, pt, px, pc, pe, pd)
            IF (model_mode == "skala-float32" .AND. e >= -1.0e-6_dp) &
               ERROR STOP "The arbitrary-point fixture must have a nontrivial Skala energy"
            IF (model_mode /= "skala-float32") THEN
               expected = 0.0_dp
               err = 0.0_dp
               descriptor_derivative = 0.0_dp
               first = 1
               spin_scale = 1.0_dp
               IF (nspin == 1) spin_scale = 0.5_dp
               DO b = 1, na_local
                  mean = 0.0_dp
                  mean_point = 0.0_dp
                  centre_covector = 0.0_dp
                  sum_weights = SUM(energy_weights(first:first + sizes(b) - 1))
                  sum_descriptor = SUM(descriptor_weights(first:first + sizes(b) - 1))
                  DO row = first, first + sizes(b) - 1
                     mean = mean + SUM(r(row, :))*descriptor_weights(row)/sum_descriptor
                     IF (model_mode == "geometry") THEN
                        mean_point = mean_point + (points(:, row) - centres(:, b))*descriptor_weights(row)/sum_descriptor
                     END IF
                  END DO
                  DO row = first, first + sizes(b) - 1
                     w = energy_weights(row)
                     dw = descriptor_weights(row)
                     energy_density = spin_scale*(0.7_dp*SUM(r(row, :)**2) + &
                                                  0.3_dp*SUM(t(row, :)**2) + 0.2_dp*SUM(r(row, :)*t(row, :)) + &
                                                  0.1_dp*SUM(g(row, :, :)**2)) + 0.05_dp*mean**2
                     weight_covector = 0.1_dp*mean*sum_weights*(SUM(r(row, :)) - mean)/sum_descriptor
                     point_covector = 0.0_dp
                     IF (model_mode == "geometry") THEN
                        displacement = points(:, row) - centres(:, b)
                        energy_density = energy_density + 0.04_dp*SUM(displacement**2) + 0.03_dp*SUM(mean_point**2)
                        point_covector = 0.08_dp*w*displacement + 0.06_dp*sum_weights*dw*mean_point/sum_descriptor
                        weight_covector = weight_covector + &
                                          0.06_dp*sum_weights*SUM(mean_point*(displacement - mean_point))/sum_descriptor
                     END IF
                     expected = expected + w*energy_density
                     centre_covector = centre_covector - point_covector
                     err = MAX(err, ABS(pe(row) - energy_density), ABS(pd(row) - weight_covector), &
                               MAXVAL(ABS(px(:, row) - point_covector)))
                     rho_sum = 0.1_dp*mean*sum_weights*dw/sum_descriptor
                     DO spin = 1, nspin
                        err = MAX(err, ABS(pr(row, spin) - &
                                           (w*spin_scale*(1.4_dp*r(row, spin) + 0.2_dp*t(row, spin)) + rho_sum)))
                        err = MAX(err, ABS(pt(row, spin) - &
                                           w*spin_scale*(0.6_dp*t(row, spin) + 0.2_dp*r(row, spin))))
                        err = MAX(err, MAXVAL(ABS(pg(row, :, spin) - w*spin_scale*0.2_dp*g(row, :, spin))))
                     END DO
                  END DO
                  err = MAX(err, MAXVAL(ABS(pc(:, b) - centre_covector)))
                  descriptor_derivative = MAX(descriptor_derivative, MAXVAL(ABS(pr(first, :))))
                  first = first + sizes(b)
               END DO
               CALL comm%sum(expected)
               CALL comm%max(err)
               CALL comm%max(descriptor_derivative)
               IF (MAX(ABS(e - expected), err) > 1.0e-11_dp) &
                  ERROR STOP "Arbitrary-point analytic energy or raw covector is incorrect"
               IF (descriptor_derivative < 1.0e-10_dp) &
                  ERROR STOP "Zero-energy-weight rows lost their descriptor derivative"
            END IF
            gathered = 0.0_dp
            DO row = 1, nr_local
               gathered(ids(row), 1, :) = pr(row, :)
               gathered(ids(row), 2:4, :) = pg(row, :, :)
               gathered(ids(row), 5, :) = pt(row, :)
               gathered(ids(row), 6:8, 1) = px(:, row)
               gathered(ids(row), 9, 1) = pe(row)
               gathered(ids(row), 10, 1) = pd(row)
            END DO
            first = 1
            DO b = 1, na_local
               gathered(ids(first), 11:13, 1) = pc(:, b)
               first = first + sizes(b)
            END DO
            CALL comm%sum(gathered)
            IF (distribution == 0) THEN
               ref_e = e
               reference_grad = gathered
            ELSE
               IF (MAX(ABS(e - ref_e), MAXVAL(ABS(gathered - reference_grad))) > routing_tolerance) &
                  ERROR STOP "Redistributing complete atom blocks changed Skala energy or derivatives"
            END IF
            point = SUM(px, DIM=2) + SUM(pc, DIM=2)
            CALL comm%sum(point)
            IF (MAXVAL(ABS(point)) > derivative_tolerance) &
               ERROR STOP "Point and block-centre derivatives violate translation invariance"
            ad = [SUM(pr*dr), SUM(pg*dg), SUM(pt*dt), SUM(px*dx), SUM(pc*dc), SUM(pe*de), SUM(pd*dd)]
            CALL comm%sum(ad)
            DO c = 1, 7
               DO ih = 1, 2
                  step_h = step_size/REAL(ih, dp)
                  SELECT CASE (c)
                  CASE (1)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r + step_h*dr, g, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r - step_h*dr, g, t, xc, minus, pr, pg, pt)
                  CASE (2)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g + step_h*dg, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g - step_h*dg, t, xc, minus, pr, pg, pt)
                  CASE (3)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t + step_h*dt, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t - step_h*dt, xc, minus, pr, pg, pt)
                  CASE (4)
                     CALL skala_point_eval(comm, points + step_h*dx, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points - step_h*dx, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt)
                  CASE (5)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres + step_h*dc, &
                                           r, g, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres - step_h*dc, &
                                           r, g, t, xc, minus, pr, pg, pt)
                  CASE (6)
                     CALL skala_point_eval(comm, points, energy_weights + step_h*de, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights - step_h*de, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt)
                  CASE (7)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights + step_h*dd, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights - step_h*dd, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt)
                  END SELECT
                  err = ABS((plus - minus)/(2.0_dp*step_h) - ad(c))
                  maximum_error = MAX(maximum_error, err)
                  IF (err > derivative_tolerance) THEN
                     WRITE (*, '(A,4I3,3ES18.8)') "Point FD failure (spin,distribution,component,step,AD,FD): ", &
                        nspin, distribution, c, ih, step_h, ad(c), (plus - minus)/(2.0_dp*step_h)
                     ERROR STOP "Arbitrary-point derivative fails finite differences"
                  END IF
               END DO
            END DO
            energy_weights = 0.0_dp
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, pr, pg, pt, px, pc, pe, pd)
            err = ABS(e)
            ! The energy-weight derivative itself need not vanish at zero weights.
            IF (nr_local > 0) err = MAX(err, MAXVAL(ABS(pr)), MAXVAL(ABS(pg)), MAXVAL(ABS(pt)), &
                                        MAXVAL(ABS(px)), MAXVAL(ABS(pc)), MAXVAL(ABS(pd)))
            CALL comm%max(err)
            IF (err > 1.0e-12_dp) ERROR STOP "Zero energy weights must give zero energy, field and geometry covectors"
            DEALLOCATE (sizes, centres, points, ids, energy_weights, descriptor_weights, &
                        r, g, t, dr, dg, dt, pr, pg, pt, dx, dc, de, dd, px, pc, pe, pd)
         END DO
         DEALLOCATE (gathered, reference_grad)
      END DO
      IF (comm%mepos == 0) WRITE (*, '(A)') "Arbitrary-point atom-block energy and derivative checks passed."
   END SUBROUTINE arbitrary_points

   SUBROUTINE atomic_fields()
      INTEGER                                            :: a, ix, iy, iz, row, spin, sx, sy, sz
      REAL(dp)                                           :: charge, displacement(3), occupation, &
                                                            point(3), radius

      ! Occupied hydrogenic orbitals and periodic images give mutually consistent primitives.
      rho = 0.0_dp
      grad = 0.0_dp
      tau = 0.0_dp
      row = 0
      DO iz = first_z, last_z
         DO iy = 0, grid%npts(2) - 1
            DO ix = 0, grid%npts(1) - 1
               row = row + 1
               point = MATMUL(grid%dh, REAL([ix, iy, iz], dp))
               DO a = 1, SIZE(atoms)
                  DO sz = -1, 1
                     DO sy = -1, 1
                        DO sx = -1, 1
                           displacement = point - atoms(a)%r - MATMUL(hmat, REAL([sx, sy, sz], dp))
                           radius = SQRT(SUM(displacement**2))
                           charge = EXP(-2.0_dp*radius)/ACOS(-1.0_dp)
                           DO spin = 1, ns
                              occupation = 1.0_dp
                              IF (ns == 2) occupation = 0.5_dp + (3 - 2*spin)*(0.1_dp + 0.05_dp*a)
                              rho(row, spin) = rho(row, spin) + occupation*charge
                              tau(row, spin) = tau(row, spin) + 0.5_dp*occupation*charge
                              grad(row, :, spin) = grad(row, :, spin) - &
                                                   2.0_dp*occupation*charge*displacement/radius
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               DO spin = 1, ns
                  drho(row, spin) = 0.01_dp*rho(row, spin)*COS(REAL(row + 64*first_z + spin, dp))
                  dtau(row, spin) = 0.02_dp*tau(row, spin)*SIN(REAL(row + 64*first_z + spin, dp))
                  dgrad(row, :, spin) = 0.01_dp*grad(row, :, spin)* &
                                        COS(REAL(row + 64*first_z + spin + [1, 2, 3], dp))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE atomic_fields

   SUBROUTINE evaluate(r, g, t, e)
      REAL(dp), INTENT(IN)                               :: r(:, :), g(:, :, :), t(:, :)
      REAL(dp), INTENT(OUT)                              :: e

      CALL skala_grid_eval(grid, r, g, t, xc, atoms, cell, e, vrho, vgrad, vtau)
      IF (.NOT. ieee_is_finite(e) .OR. .NOT. ALL(ieee_is_finite(vrho)) .OR. &
          .NOT. ALL(ieee_is_finite(vgrad)) .OR. .NOT. ALL(ieee_is_finite(vtau))) &
         ERROR STOP "Nonfinite model energy or primitive derivative"
   END SUBROUTINE
END PROGRAM
