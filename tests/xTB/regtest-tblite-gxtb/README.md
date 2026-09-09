# Native periodic g-xTB regression provenance

The energy references were refreshed on 2026-09-09 for the minimal
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

These small-system regression tests do not establish performance, dense-mesh
benchmark convergence, or qualification of optional acceleration combinations.
