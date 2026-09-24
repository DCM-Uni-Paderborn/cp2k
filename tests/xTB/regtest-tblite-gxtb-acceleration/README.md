# Extended-provider g-xTB tests

This directory requires the `tblite_gxtb_acceleration` capability reported by
`cp2k --version`. CMake enables that capability only after compiling a probe
against the extended save_tblite API and checking its partial-exchange ABI.
The ordinary `tblite_gxtb` feature remains available with the minimal provider.

The tests cover explicit production and dense-oracle qualification modes for
distributed-image exchange, symmetry-fused traversal, Fourier transforms,
streamed ACP response, projector caching, and cross-mesh restarts. Restart
producer/consumer pairs remain in the same directory and in dependency order.

`CH4_gxtb_kp_symmetry_fused_production.inp` is deliberately separate from its
`QUALIFY` counterpart: its compact production path must not allocate the dense
qualification oracle as an unnoticed prerequisite.

During a provider/model port, energy references must be checked against the
same model and geometry on the independently qualified dense path. Old-model
references are not evidence of equivalence, and numerical gates must not be
relaxed merely to make a port pass.

## Historical energy references (9 September 2026)

The references were independently recalculated with the previously qualified
CP2K `93a083bdebc40c7c31c19678bbb6bbc80bb3356e` and the minimal
`lmseidler-integration` provider `c85459da430ffd9cf8992bcd929d5c435f64eb0b`.
The dense controls remove only `GXTB_ACCELERATION` and rename `PROJECT`;
the 4x4x4 control changes only the mesh of the 3x3x3 control. Geometries,
model parameters and SCF thresholds are unchanged.

| Independent control | Energy / Eh |
| --- | ---: |
| CH4, 2x2x2 | -40.467010243469048 |
| CH4, 3x3x3 | -40.466670333868578 |
| CH4, 4x4x4 | -40.466666199913966 |
| Si, shifted 2x2x2 | -579.001230548395938 |

The acceleration/restart combination is CP2K
`531ff1324aadf766f75e30327264e1988fda9169` with performance provider
`9f1c5f83463e933da9a32cc42360147d25602a0f`. All thirteen explicit serial
acceleration/restart controls agree with the independently qualified dense
controls (H2O restart is checked against its qualified density-mixer baseline).
Six controls each also pass with two and four MPI ranks; actual inter-k-point
groups and distinct singleton CPU bindings were verified. These small-system
checks establish numerical equivalence, not a large-system speedup factor.

Frozen executable SHA-256 values:

- Dense reference: `499782fc7464b033cf410790c2ccee93fc18b31fafb54a790085e5318147243b`.
- Performance combination: `a5c73193439aaac3b77351140d2b488716395ad2195782e7e63322b676b2bb8a`.

The complete dynamically loaded runtime closure was verified before and after
each run, not just the executable. The nine old-model energy references
in this directory and the two CH4 restart tests were refreshed only after
these independent controls passed. No matcher tolerance was changed.

## Corrected periodic-model references (25 September 2026)

The current references include the periodic Hubbard/self-image, smooth-WS and
multicharge image-sum corrections described in the adjacent native suite.
All nine inputs were repeated with the extended provider, separately at
1 MPI x 1 OpenMP, 4 x 1 and 1 x 4. Independent controls with the minimal
provider remove only `GXTB_ACCELERATION`; coordinates, mesh, SCF thresholds,
printing and restart producer/consumer relationships are unchanged.

| Independent dense control | Energy / Eh |
| --- | ---: |
| CH4, 2x2x2 | -40.467008824367582 |
| CH4, 3x3x3 | -40.466668927226940 |
| CH4, 4x4x4 restart | -40.466664793424440 |
| Si, shifted 2x2x2 | -579.001250005684600 |

Across all 27 extended comparisons to the nine dense controls, the maximum
energy difference is 1.43e-14 Eh, force-component difference 2.88e-16 Eh/bohr,
and stress-component difference 5.08e-10 bar (derivatives at printed precision).
All non-energy qualification/production/restart matchers pass. The old energy
reference failures are preserved in the audit, and every matcher tolerance
is unchanged. In particular, this is not a change to the accelerated model
to fit the old energies, nor an acceptance of stale-model energy references.

Frozen executable SHA-256 values:

- Dense minimal: `db56ac410a64236de86611fd0f4327646a8eea6135fc2bc331616c16290b92ea`.
- Extended: `8617dd64636f92d3e5cfa0e05ef8adbdf78ae0dabbf807b35e39b1b3612c42c1`.

Input/output and complete loaded-runtime hashes, CPU reservations and pre-exec
affinity proofs accompany every result. Evidence is in
`current-acceleration-reference-review.json` and the four source reports.
These tests qualify this combination on the listed systems; performance is
measured separately, not inferred from passing numerical checks.
