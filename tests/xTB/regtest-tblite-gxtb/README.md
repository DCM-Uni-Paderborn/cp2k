# Native periodic g-xTB regression provenance

## Historical 9 September baseline

The prior energy references were refreshed on 2026-09-09 for the minimal
`DCM-Uni-Paderborn/save_tblite:lmseidler-integration` model at
`02730ac0ada8a56474fe53ffd6634e83e8cc2661`. They are not references for the
older performance branch's molecular parameterization.

Qualification used CP2K source `5a9779469bff065061fbf5974f8ef22245c9bf15`,
including the current-upstream k-point operator callback and the separation of
geometry-only and density-dependent repulsion energies. The tested executable's
SHA-256 was `c481326cca1789739f154863ca53e7e0634c5d1e54d0b4723c3aedda459a0628`.
All 34 inputs in this directory and the adjacent SPGLIB suite terminated
normally with this executable; energy, derivative, population, mixer-path and
special-k-point checks were evaluated from those outputs.

Independent consistency gates included direct CLI energy/force/stress checks,
implicit versus explicit Gamma, full versus K290/SPGLIB meshes, shifted silicon
meshes, spin-polarized time reversal, and an explicit H2 Gamma supercell versus
the equivalent 3x1x1 mesh. Analytic force/stress finite differences include an
explicit 2x2x2 SPGLIB-reduced methane calculation. Existing numeric tolerances
were retained; model changes must not be hidden by increasing them.
The 2x2x2 finite-difference case checks the initial special-k-point count:
the displaced geometries legitimately have lower symmetry.

Those small-system regression tests did not establish performance, dense-mesh
benchmark convergence, or qualification of optional acceleration combinations.

## Gamma matrix support regression

The two `NaCl_*_611_support` inputs are a deliberately distorted cubic salt,
not a rock-salt benchmark. They compare a Gamma-centered 6x1x1 primitive mesh
with its explicit 12-atom Gamma supercell. Ordinary Gaussian overlap screening
omitted Gamma atom blocks required by the separable ACP and exchange density.
The k-point path already retained the provider's convolution support.

The corrected pair agrees within 4.1e-11 Eh per primitive cell. The old Gamma
result differs by 5.65e-5 Eh for the supercell, so the regression distinguishes
the fix without relaxing energy tolerances. The even mesh must explicitly be
Gamma-centered; an unshifted MACDONALD input alone is not sufficient.
These references require the accompanying periodic Coulomb/auxiliary-charge
corrections in the provider and dependency.

## 25 September corrected periodic model

Current energy references include the provider's full periodic Hubbard
correction, retained nonzero-image self terms, the differentiable C2
Wigner-Seitz partition, and corrected multicharge real/reciprocal/image bounds.
The latter is distributed as a pinned source plus a tested correctness patch
in save_tblite, not as an unpublished dependency assumption.

The minimal and extended provider candidates were tested with the same
CP2K source. All 39 inputs in this directory and the SPGLIB suite terminate
normally with converged SCF. Their maximum cross-build energy difference is
2.85e-14 Eh. Executable SHA-256 identifiers:

- Minimal: `db56ac410a64236de86611fd0f4327646a8eea6135fc2bc331616c16290b92ea`.
- Extended: `8617dd64636f92d3e5cfa0e05ef8adbdf78ae0dabbf807b35e39b1b3612c42c1`.

References were reviewed against full/K290/SPGLIB, general/shifted mesh,
time-reversal, mixer/restart, CLI-Gamma and explicit-supercell comparisons.
The independent provider and multicharge image-sum/replication tests also pass.
The old numeric failures were retained in the validation archive. Changing
the references does not assert that the historical and corrected periodic
models are identical; no numeric tolerance was increased.

The C2 methane finite-difference step check was repeated for implicit Gamma,
explicit Gamma, SPGLIB Gamma and SPGLIB 2x2x2. Their virial difference sums
at DX=1e-4 are 1.216421e-6, 1.216958e-6, 1.218570e-6 and 1.237255e-6 Eh.
At DX=5e-5 they fall to 1.57812e-7, 1.59961e-7, 1.53515e-7 and 1.61995e-7;
at DX=2.5e-5 to 2.1391e-8, 2.1391e-8, 1.2797e-8 and 2.1277e-8 Eh.
Both providers pass the unchanged 1e-6 tolerance at the retained DX=5e-5.
The extended-provider corresponding maximum is 1.60921e-7 Eh.

## Smooth-image stress finite differences

For the cubic methane K290 case, use `DX 5e-5` with the existing `1e-6`
virial/force tolerances. With the current C2 Wigner-Seitz partition, the
virial absolute-difference sum at steps `1e-4`, `5e-5`, and `2.5e-5` is
respectively `1.238329e-6`, `1.60921e-7`, and `1.9129e-8` Eh. Tightening
`EPS_SCF` from `1e-10` to `1e-12` leaves the first value essentially
unchanged (`1.235106e-6`); a full, unreduced mesh also reproduces it
(`1.239403e-6`). This step-convergence check distinguishes finite-difference
discretization from an analytic-stress or symmetry-reduction discrepancy.
The reported quantity is the energy-valued `pv_virial` difference, not a
volume-normalized stress.

## Controlled NaCl timing check, 25 September

The corrected builds above were measured on Terok with the same full 6x6x6
input (SHA-256 `01acac45ebdf316d2c42580eedd6f14fb9090c61340586ef4f0e6cd3e8cf3885`).
Two interleaved repetitions per configuration ran without another native
calculation/compiler at launch. BLAS remained single-threaded, CPU affinity
was disjoint and reserved, and all energies were -622.6304237931206 Eh.

| Provider | MPI x OpenMP | Median elapsed compute / s | Sampled rank-summed peak RSS / GiB |
| --- | ---: | ---: | ---: |
| Minimal | 1 x 1 | 73.85 | 0.746-0.747 |
| Minimal | 1 x 4 | 65.28 | 0.748 |
| Minimal | 4 x 1 | 54.45 | 2.249-2.345 |
| Extended | 1 x 1 | 70.64 | 0.749-0.759 |
| Extended | 1 x 4 | 64.64 | 0.755-0.765 |
| Extended | 4 x 1 | 49.04 | 2.350 |

OpenMP speedup is modest (1.13x minimal, 1.09x extended); four MPI ranks
give 1.36x and 1.44x, respectively. The extended provider is not a universal
large speedup over the corrected minimal provider. Both share the CP2K-side
allocation/layout reuse, bounded overlap reuse and regular-grid Fock/overlap
FFTs. The earlier extended batched implementation redundantly repeated the
same overlap transform for each batch; avoiding that is not an approximation
to the exchange model. Global exchange coupling, density transforms and MPI
matrix communication still limit scaling, especially for a tiny AO basis.

These are two-repeat observations on a shared host, not confidence intervals
or a prediction for an 8x8x8 benchmark. RSS is sampled every 0.2 s, may miss
short peaks and counts shared pages once per rank; it is not PSS. Existing
historical timings use different code/physics and must not be treated as a
controlled same-build speedup series. Evidence: `current-controlled-k6-timings.json`
and its twelve raw run directories in the frozen validation archive.
