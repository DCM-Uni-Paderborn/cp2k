# Skala through CP2K and SIRIUS

These manual integration fixtures require LibTorch and the matching SIRIUS
external-XC APIs. They are not registered in `TEST_FILES.toml`. Run them in a
separate output directory with `SKALA_MODEL` pointing to the CPU TorchScript
model and `CP2K_DATA_DIR` to CP2K's data directory.

SIRIUS supplies densities, gradients and kinetic-energy densities and applies the
returned variational operators. CP2K evaluates Skala and its derivatives through
LibTorch. CPU/GPU orbitals and model inference are selected independently, with
host-staged fields and FP64 point evaluation. Parameter values, serialized models
and native CP2K evaluation defaults are unchanged.

The regular-grid energy callback supports norm-conserving pseudopotentials without
nonlinear core corrections. `PW_DFT/SKALA_GRID` selects atom quadrature for PAW and
FP-LAPW and is required for forces, stress and Skala CUDA inference in all three
methods. PAW and norm-conserving nonlinear core corrections require setup-consistent
positive core kinetic density. FP-LAPW requires a `SIRIUS_ATOM_FILE` in each `KIND`,
not a UPF file or Gaussian basis.

## Checks

The periodic atom-quadrature and adjoint checks run in CP2K's standard unit tests.
From the source root, replacing `build` with the build directory:

```sh
python tools/regtesting/create_skala_grid_test_model.py build/grid-test.fun
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-test.fun
python tools/regtesting/create_skala_grid_test_model.py build/grid-geometry-test.fun --geometry
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-geometry-test.fun geometry
OMP_NUM_THREADS=1 build/bin/skala_sirius_unittest.psmp build/grid-test.fun
python tools/regtesting/check_sirius_lapw_input.py build/bin/cp2k.psmp build/lapw-input-guards
python tools/regtesting/test_sirius_lapw_test_atom.py
```

These tests cover primitive fields, coordinate/weight derivatives, reconstruction,
orbital derivatives, SCF energy accounting and input validation. Repeat the
executables with `mpiexec -n 2`. The input checker accepts
`--launcher 'mpirun -n 2'` and requires a new output directory.
The generated analytic models are test fixtures, not Skala. Coordinate/weight
checks hold the fields fixed. The force and stress modes below check total derivatives.

## Actual-Model Inputs

| Input | Coverage | Required atom data |
| --- | --- | --- |
| `H2.inp` | Norm-conserving molecular SCF on the regular grid | CP2K data directory |
| `H2-spinor.inp` | Twelve noncollinear iterations with nonparallel initial spins and GPU orbitals/inference | CP2K data directory |
| `H-spin-kpoints.inp` | Collinear norm-conserving energy, forces and stress with shifted k points | CP2K data directory |
| `H-paw-atom-grid.inp` | Core-free PAW with atom quadrature | `H-paw-no-core.json` |
| `C-paw-core.inp` | Collinear PAW energy, forces and stress with two core and four valence electrons | `C-paw-core.json` |
| `He-lapw-atom-grid.inp` | All-valence FP-LAPW energy, forces and stress | `He-lapw.json` |
| `Ne-lapw-atom-grid.inp` | FP-LAPW energy, forces and stress with an explicit 1s core | `Ne-lapw.json` |

Generate LAPW data with `create_sirius_lapw_test_atom.py He-lapw.json` or
`create_sirius_lapw_test_atom.py Ne-lapw.json --element Ne`.
For PAW hydrogen, start from SIRIUS's
`examples/pp-pw/C42H58ClNO2Si_420_atoms/H.pz-kjpaw_psl.0.1.UPF.json` and add
an explicit zero `ae_core_kinetic_density` array inside `pseudo_potential.paw_data`,
matching its core-charge array length, after verifying that the setup is core-free.
The carbon fixture requires a four-valence-electron PAW setup with matching
two-electron AE core charge and positive tau. Convert matching UPF and ATOMPAW XML
exports from the same atomic calculation, for example the
[JTH PAW datasets](https://www.abinit.org/atomic_data/paw/):

```sh
upf_to_json C.UPF --paw-core-xml C.xml -o C-paw-core.json
```

The importer validates the element, valence charge and both core charge profiles.
Run inputs with `OMP_NUM_THREADS=1 cp2k.psmp -i input.inp -o output.out`.
The fixtures test implementation, not physical cutoff or quadrature convergence.

PW/PAW forces and stress require the point-geometry API. Nonrelativistic FP-LAPW
additionally requires joint Hessian and quadrature-direction APIs, with full
variation for scalar/collinear states and second variation for spinors.
The spinor API includes PW/PAW SOC and nonrelativistic FP-LAPW second-variation SOC.
Actual-model noncollinear SCF remains experimental. `H2-spinor.inp` exercises the
iteration path without asserting SCF convergence.
`KIND/SIRIUS_MAGNETIZATION_VECTOR` sets the initial Cartesian magnetization.
Without it, the scalar `MAGNETIZATION` remains along z.

Noncollinear Skala uses the global multicollinear projection of the complete
energy, following [Pu et al.](https://doi.org/10.1103/PhysRevResearch.5.013036).
`SKALA_GRID/N_SPIN_ANGULAR` selects a positive Lebedev spin-direction rule
(default 194, up to 974), independently of the real-space `N_ANGULAR` rule.
Energy, covectors and Hessian actions use the same spin rule.

## Derivative Modes

The `skala_sirius_unittest.psmp` arguments are
`MODEL SETUP SIDE STEP MODE PRESET ORBITALS INFERENCE FIELD_STEP BANDS`.
`SETUP` is `LAPW` for the built-in all-valence FP-LAPW fixture, a PAW JSON file
for a supplied setup, or empty for the synthetic PW/PAW fixtures.
`ORBITALS` and `INFERENCE` independently select `gpu`; omitted values select CPU.
The legacy `skala-float32` preset selects actual-model steps and solver thresholds,
not point-evaluation precision. Both presets retain the same derivative-error bounds.
`FIELD_STEP` defaults to `STEP`, and `BANDS` to four.

```sh
skala_sirius_unittest.psmp SKALA C-paw-core.json 4 0.00001 atom-grid-orbitals skala-float32 gpu gpu 0.00001 8
skala_sirius_unittest.psmp SKALA LAPW 18 "" atom-grid-forces skala-float32
skala_sirius_unittest.psmp SKALA LAPW 4 0.00001 spinor-lapw-soc-orbitals skala-float32 gpu gpu
```

`atom-grid-orbitals` checks scalar/collinear orbital derivatives. Spinor modes are
`spinor-orbitals`, `spinor-soc-orbitals`, `spinor-paw-orbitals`,
`spinor-paw-soc-orbitals`, `spinor-atom-grid-orbitals` and `spinor-lapw-soc-orbitals`.
These hold occupations, radial functions and core states fixed and do not require
SCF convergence. PAW rotations need not preserve its overlap norm.
For the synthetic PAW SOC mode, use `STEP=0.000002`.
`atom-grid-forces` and `atom-grid-stress` check total derivatives with geometry
updates and coupled LAPW response. `spinor-stress` checks fixed-orbital XC stress.
Each enabled derivative check uses two step sizes without density renormalization.

`skala_grid_unittest.psmp SKALA skala-float32 gpu` checks primitive-field derivatives
and matched CPU/GPU energies, covectors and Hessians. Use `spinor-float32` for spin
rotations, collinear limits and joint field/geometry Hessian actions.
Generate an analytic spin fixture with `--geometry --spin` and use `spinor`.
Adding `--float32` and selecting `spinor-precision` tests promotion of parameters,
floating-point buffers and casts without changing integer index buffers.

On macOS, pass any required `DYLD_LIBRARY_PATH` through `mpiexec ... env` so
all ranks load the same OpenMP runtime. Do not use `KMP_DUPLICATE_LIB_OK`.
