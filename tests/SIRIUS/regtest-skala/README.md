# Skala through CP2K and SIRIUS

These manual fixtures require LibTorch and matching SIRIUS external-XC APIs. They are not registered
in `TEST_FILES.toml`. Run in a separate output directory with `SKALA_MODEL` pointing to the CPU
TorchScript model and `CP2K_DATA_DIR` to CP2K's data directory.

SIRIUS supplies densities, gradients and kinetic-energy densities and applies the returned
variational operators. CP2K evaluates Skala and its derivatives through LibTorch. CPU/GPU orbital
processing and model inference are independent, with host-staged fields and FP64 point evaluation.

The regular-grid callback supports norm-conserving pseudopotentials without nonlinear core
corrections. `PW_DFT/SKALA_GRID` selects atom quadrature, required for PAW, FP-LAPW, forces, stress
and CUDA inference. Nonlinear core corrections require setup-consistent positive core kinetic
density. FP-LAPW requires `KIND/SIRIUS_ATOM_FILE`, not a UPF file or Gaussian basis.

## Unit Tests

Atom-quadrature and adjoint checks run in CP2K's standard unit tests. Additional analytic fixtures:

```sh
python tools/regtesting/create_skala_grid_test_model.py build/grid-test.fun
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-test.fun
python tools/regtesting/create_skala_grid_test_model.py build/grid-geometry-test.fun --geometry
OMP_NUM_THREADS=1 build/bin/skala_grid_unittest.psmp build/grid-geometry-test.fun geometry
OMP_NUM_THREADS=1 build/bin/skala_sirius_unittest.psmp build/grid-test.fun
python tools/regtesting/check_sirius_lapw_input.py build/bin/cp2k.psmp build/lapw-input-guards
python tools/regtesting/test_sirius_lapw_test_atom.py
```

Repeat executables with `mpiexec -n 2`. The input checker accepts `--launcher 'mpirun -n 2'` and
requires a new output directory. Generated models are analytic test fixtures, not Skala.

## Skala Inputs

| Input                   | Coverage                                                                                 | Required atom data   |
| ----------------------- | ---------------------------------------------------------------------------------------- | -------------------- |
| `H2.inp`                | Norm-conserving molecular SCF on the regular grid                                        | CP2K data directory  |
| `H2-spinor.inp`         | Twelve noncollinear iterations with nonparallel initial spins and GPU orbitals/inference | CP2K data directory  |
| `H-spin-kpoints.inp`    | Collinear norm-conserving energy, forces and stress with shifted k points                | CP2K data directory  |
| `H-paw-atom-grid.inp`   | Core-free PAW with atom quadrature                                                       | `H-paw-no-core.json` |
| `C-paw-core.inp`        | Collinear PAW energy, forces and stress with two core and four valence electrons         | `C-paw-core.json`    |
| `He-lapw-atom-grid.inp` | All-valence FP-LAPW energy, forces and stress                                            | `He-lapw.json`       |
| `Ne-lapw-atom-grid.inp` | FP-LAPW energy, forces and stress with an explicit 1s core                               | `Ne-lapw.json`       |

Generate LAPW data with `create_sirius_lapw_test_atom.py He-lapw.json` or select `--element Ne` for
an explicit 1s core (`H` and `Li` are also supported). For core-free PAW hydrogen, start from
SIRIUS's `examples/pp-pw/C42H58ClNO2Si_420_atoms/H.pz-kjpaw_psl.0.1.UPF.json` and add an explicit
zero `ae_core_kinetic_density` array inside `pseudo_potential.paw_data`, matching its core-charge
array length. Carbon requires a four-valence-electron PAW setup with matching two-electron AE core
charge and positive tau. Convert matching UPF and ATOMPAW XML exports from the same atomic
calculation, for example the [JTH PAW datasets](https://www.abinit.org/atomic_data/paw/):

```sh
upf_to_json C.UPF --paw-core-xml C.xml -o C-paw-core.json
```

Run with `OMP_NUM_THREADS=1 cp2k.psmp -i input.inp -o output.out`.

Forces/stress require the point-geometry API, with joint Hessian and quadrature-direction APIs for
FP-LAPW. Scalar/collinear states support scalar-relativistic full variation. Spinors support NR/ZORA
second variation and SOC; tau-dependent explicit core response supports NR/ZORA core states.
Noncollinear Skala SCF is experimental. `KIND/SIRIUS_MAGNETIZATION_VECTOR` sets the initial
Cartesian magnetization; scalar `MAGNETIZATION` initializes its z component.

Noncollinear Skala uses multicollinear projection of the complete energy, following
[Pu et al.](https://doi.org/10.1103/PhysRevResearch.5.013036). `SKALA_GRID/N_SPIN_ANGULAR` selects a
positive Lebedev spin-direction rule (default 194, up to 974), independently of the real-space
`N_ANGULAR` rule. Energy, covectors and Hessian actions use the same spin rule.

## Derivative Tests

The `skala_sirius_unittest.psmp` arguments are
`MODEL SETUP SIDE STEP MODE PRESET ORBITALS INFERENCE FIELD_STEP BANDS RELATIVITY ATOM_FILE CORE_RELATIVITY`.
`SETUP` is `LAPW`, a PAW JSON file, or empty for synthetic PW/PAW fixtures. `ORBITALS` and
`INFERENCE` select `gpu` independently (default CPU). `skala-float32` selects actual-model steps and
solver thresholds, not point precision or derivative-error bounds. Defaults are `FIELD_STEP=STEP`,
`BANDS=4`, `RELATIVITY=none` and `CORE_RELATIVITY=none` (also accepts `zora`). With `SETUP=LAPW`,
`ATOM_FILE` optionally replaces the all-valence fixture for explicit-core tests.

```sh
skala_sirius_unittest.psmp SKALA C-paw-core.json 4 0.00001 atom-grid-orbitals skala-float32 gpu gpu 0.00001 8
skala_sirius_unittest.psmp SKALA LAPW 18 "" atom-grid-forces skala-float32
skala_sirius_unittest.psmp SKALA LAPW 4 0.00001 spinor-lapw-soc-orbitals skala-float32 gpu gpu 0.00001 4 zora
```

`atom-grid-orbitals` checks scalar/collinear orbital derivatives. Spinor modes are
`spinor-orbitals`, `spinor-soc-orbitals`, `spinor-paw-orbitals`, `spinor-paw-soc-orbitals`,
`spinor-atom-grid-orbitals` and `spinor-lapw-soc-orbitals`. Orbital checks hold occupations, radial
functions and core states fixed, using a compact 12/26 real-space grid and 26 spin directions. PAW
rotations need not preserve its overlap norm. Synthetic PAW SOC uses `STEP=0.0000005`.
`atom-grid-forces`, `atom-grid-stress`, `spinor-lapw-soc-forces` and `spinor-lapw-soc-stress` check
total derivatives with coupled geometry response. `spinor-stress` checks fixed-orbital XC stress.

`skala_grid_unittest.psmp SKALA skala-float32 gpu` checks fields and CPU/GPU covectors/Hessians. Use
`spinor-float32` for spin rotations, collinear limits and joint field/geometry Hessians. Analytic
spin fixtures use generator flags `--geometry --spin` and test mode `spinor`; add `--float32` and
use `spinor-precision` to check model precision promotion.

On macOS, pass any required `DYLD_LIBRARY_PATH` through `mpiexec ... env` so all ranks load the same
OpenMP runtime. Do not use `KMP_DUPLICATE_LIB_OK`.
