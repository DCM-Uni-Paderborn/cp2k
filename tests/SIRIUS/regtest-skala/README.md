# Skala through CP2K and SIRIUS

These manual integration fixtures require LibTorch and SIRIUS's external-XC
API. They are not registered in `TEST_FILES.toml`. Run them in a separate output
directory with `GAUXC_SKALA_MODEL` pointing to the CPU TorchScript model and
`CP2K_DATA_DIR` to CP2K's data directory.

The regular-grid callback supports norm-conserving pseudopotentials without
nonlinear core corrections. `PW_DFT/SKALA_GRID` selects reconstructed fields
on atom quadratures for PAW and FP-LAPW. PAW requires setup-consistent positive
core kinetic density. FP-LAPW requires a `SIRIUS_ATOM_FILE` in each `KIND`,
not a UPF file or Gaussian basis.

## Model-independent checks

From the source root, replacing `build` with the build directory:

```sh
python tools/regtesting/create_skala_grid_test_model.py build/grid-test.fun
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-test.fun
python tools/regtesting/create_skala_grid_test_model.py build/grid-geometry-test.fun --geometry
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-geometry-test.fun geometry
OMP_NUM_THREADS=1 build/bin/skala_atom_grid_unittest.psmp
OMP_NUM_THREADS=1 build/bin/skala_sirius_unittest.psmp build/grid-test.fun
python tools/regtesting/check_sirius_lapw_input.py build/bin/cp2k.psmp build/lapw-input-guards
python tools/regtesting/test_sirius_lapw_test_atom.py
```

The tests cover primitive-field and coordinate/weight derivatives, periodic
quadrature, reconstructed fields, orbital derivatives, SCF energy accounting
and input validation.
Repeat the three executables with `mpiexec -n 2`; four ranks additionally
exercise empty FFT slabs and atom blocks. The input checker accepts
`--launcher 'mpirun -n 2'`. Its output directory must not already exist.
The generated analytic models are test fixtures, not Skala. Coordinate and weight
derivatives hold the physical fields fixed and are not complete atomic forces.

## Actual-model inputs

| Input | Coverage | Required atom data |
| --- | --- | --- |
| `H2.inp` | Norm-conserving molecular SCF | CP2K data directory |
| `H-spin-kpoints.inp` | Collinear spin and shifted k points | CP2K data directory |
| `H-paw-atom-grid.inp` | Core-free PAW with atom quadrature | `H-paw-no-core.json` |
| `C-paw-core.inp` | Collinear PAW with two core and four valence electrons | `C-paw-core.json` |
| `He-lapw-atom-grid.inp` | All-valence FP-LAPW with atom quadrature | `He-lapw.json` |
| `Ne-lapw-atom-grid.inp` | FP-LAPW with an explicit 1s core | `Ne-lapw.json` |

Generate LAPW data with `create_sirius_lapw_test_atom.py He-lapw.json` or
`create_sirius_lapw_test_atom.py Ne-lapw.json --element Ne`.
The Ne fixture includes local 2s and 2p orbitals to resolve its valence shell
at the small test plane-wave cutoff.
For PAW hydrogen, start from SIRIUS's
`examples/pp-pw/C42H58ClNO2Si_420_atoms/H.pz-kjpaw_psl.0.1.UPF.json` and add
`paw_data.ae_core_kinetic_density` as an explicit zero array of the same
length as `ae_core_charge_density`, after checking that this setup has no core
electrons. Nonzero core charge requires independently supplied core tau.
The carbon fixture requires a four-valence-electron PAW setup whose
`paw_data.ae_core_kinetic_density` matches its two-electron AE core density.
For example, convert the matching carbon UPF and ATOMPAW XML files from the
[JTH PAW datasets](https://www.abinit.org/atomic_data/paw/) with SIRIUS:

```sh
upf_to_json C.UPF --paw-core-xml C.xml -o C-paw-core.json
```

The optional XML import verifies the element, valence charge and both core
charge profiles before adding the supplied AE and pseudo core kinetic densities.
Use both exports from the same atomic calculation. Resampling uses shared GSL.

Run an input with `OMP_NUM_THREADS=1 cp2k.psmp -i input.inp -o output.out`.
Check both SCF and band-solver convergence, and compare serial and MPI energies.
The low cutoffs and quadratures test integration, not physical convergence.
Point-quadrature forces and stress are supported for PW and PAW. Nonrelativistic
FP-LAPW forces require the joint Hessian and quadrature-direction callbacks,
with full variation for collinear states. FP-LAPW stress, GPU, SOC and
noncollinear support remain unavailable.

## Orbital derivatives

`skala_sirius_unittest.psmp SKALA LAPW 18 0.005` checks all-valence LAPW.
`skala_sirius_unittest.psmp SKALA H-paw-no-core.json 18 0.005` checks PAW.
The last argument enables finite differences at that step and half the step.
The supplied tolerance accounts for Skala's float32 arithmetic.

These variations hold occupations, radial functions and core states fixed.
PAW orbital rotations preserve the plane-wave norm but not necessarily the
PAW overlap norm. No density renormalization is applied. Primitive-field
checks use `skala_grid_unittest.psmp SKALA skala-float32`.

Total-force checks use `skala_sirius_unittest.psmp MODEL LAPW 18 "" atom-grid-forces`
with an analytic model, or append `skala-float32` when using `SKALA` as the model.

On macOS, pass any required `DYLD_LIBRARY_PATH` through `mpiexec ... env` so
all ranks load the same OpenMP runtime. Do not use `KMP_DUPLICATE_LIB_OK`.
