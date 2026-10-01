# Native scalar USPP implementation

## State and ordering

The NC implementation is submitted separately in CP2K PR #6143. This branch adds a reader for the
angular-resolved augmentation in SSSP v2.0 and a direct reciprocal augmentation reference. These
building blocks do **not** yet enable a USPP Quickstep calculation. The ordinary reader continues to
reject unsupported calculation paths unless augmentation import is explicitly requested by a caller.

The implementation order remains USPP GPW, USPP GAPW_XC, USPP GAPW, then PAW GPW, PAW GAPW_XC, and
PAW GAPW. Each step needs an energy functional, its Hamiltonian, its metric, and derivatives of that
same discretized functional. An import audit alone does not establish any of those calculation
modes.

## Dataset convention

For scalar projectors write `i=(p,l,m)`, with a radial projector index `p` and CP2K's real,
unit-normalized spherical harmonics. Preserve the tabulated beta amplitudes; the existing NC reader
removes the radial factor in PP_BETA and converts PP_DIJ from Ry to Hartree. Neither conversion
applies to augmentation.

PP_QIJL stores `r^2 Q_pq^L(r)`. The radial pair index is `q*(q-1)/2+p`, with `p<=q`. The angular
expansion is

```text
Q_ij(r) = sum_LM Gaunt(l_i m_i, l_j m_j, LM) Q_pq^L(r) Y_LM(rhat).
```

Allowed L range from `abs(l_i-l_j)` to `l_i+l_j` in steps of two. Every allowed channel must be
present, or be explicitly declared null. Duplicate channels, incorrect composite indices, invalid
angular selection rules and inconsistent array lengths are errors. The parser preserves PP_Q
independently of PP_QIJL.

The independent audit found 86 distinct USPP files and 3228 QIJL channels in the four SSSP v2.0
collections. All use `q_with_l=true, nqf=0`; both `US` and `USPP` occur in the header. The largest
difference between a stored PP_Q coefficient and a Simpson integral of its L=0 QIJL function is
2.1190930433e-5 (Mg GBRV PBE). This is a dataset/quadrature consistency observation, not a
total-energy error.

The compiled Debug reader passed the independent XML comparison for all 86 files and 3228 channels.
Five malformed-input checks reject a wrong composite index, wrong angular momentum, inconsistent tag
suffix, missing channel and duplicate channel. The fixture also covers attribute-free PP_R, PP_RAB,
PP_DIJ and PP_Q tags, null channels, and the 62-electron core of 16-valence-electron Pt. All five NC
units still pass on one and two MPI ranks after these changes.

Quantum ESPRESSO also obtains its runtime overlap coefficients from G=0 of the augmentation Fourier
transform: see its pinned
[initialization](https://github.com/QEF/q-e/blob/61569eb480231b649c47f631df9c5c83537df461/upflib/init_us_1.f90),
[radial transform](https://github.com/QEF/q-e/blob/61569eb480231b649c47f631df9c5c83537df461/upflib/qrad_mod.f90),
and
[angular reconstruction](https://github.com/QEF/q-e/blob/61569eb480231b649c47f631df9c5c83537df461/upflib/qvan2.f90).
The independent Python audit integrates the whole supplied radial mesh; QE's radial transform uses
its effective beta cutoff. A separate XML/SciPy audit found that every QIJL array in all 86 files
vanishes beyond the maximum PP_BETA cutoff_radius_index. Using the same quadrature on the full and
truncated meshes gives exactly zero difference at G=0, 1 and 5 bohr^-1 for all 3228 channels. Two
files have an even cutoff index. This establishes cutoff insensitivity for those sampled radial
transforms; it does not replace comparison of complete implementations or their SCF results.

## GPW reference: one augmentation representation

Let `B_mu,i = <chi_mu|beta_i>` and let P be the Gaussian density matrix. Per atom,
`D_ij = (B^dagger P B)_ij`, with spin components retained where required. Then

```text
n(r) = n_G(r) + sum_aij D^a_ij Q^a_ij(r)
S = S_G + sum_a B^a q^a B^{a,dagger}
q_ij = integral Q_ij(r) dr
Tr(P S) = integral n(r) dr.
```

The final equality must hold for the implemented representation, including the normalization of
CP2K's Fourier coefficients. It cannot rely on accidental agreement of PP_Q with a separately
discretized augmentation density.

A reciprocal-space reference avoids dividing QIJL by r^2 at the origin. Define

```text
Qbar_pq^L(G) = integral dr j_L(G*r) [r^2 Q_pq^L(r)]
Q_ij(G) = (4*pi/Omega) exp(-i G.R_a)
          * sum_LM (-i)^L Gaunt(i,j,LM) Y_LM(Ghat) Qbar_pq^L(G).
```

The `G=0` limit is evaluated analytically, and `q=Omega*Q(G=0)`. Collocation of `sum D Q(G)` and
integration of `V(G)` against Q must be adjoint operations with the same Fourier normalization,
radial quadrature, angular convention and cutoffs. Required tests include charge/metric equality,
Hermiticity, rotational covariance, a constant-potential shift and a random density/potential
adjoint test.

`src/uspp_augmentation.F` now implements this direct Fourier reference, including full-mesh
quadrature, analytic G=0, magnetic-projector expansion, Gaunt contractions, translation phase,
density collocation and its adjoint potential integration. Its preparation reuses native
`upf_projector_type` data with Cartesian degree at least `max(lbeta)`. The reciprocal coefficients
use integral/volume normalization; integration takes explicit conjugate-pair multiplicities.
Self-conjugate modes have unit multiplicity. The current collocation API takes the real symmetric
projector density matrix. No global Quickstep overlap or density array has yet been changed.

`uspp_fourier_unittest` compares 15 complete 19-by-19 matrices against independent Cartesian
Gaussian Fourier integrals, including s/p/d/f projectors, two distinct radial p projectors,
multipoles through L=6, odd/even uniform meshes, a shifted logarithmic mesh,
zero/small/axial/general reciprocal vectors and a displaced center. The maximum absolute discrepancy
in volume times Q(G) is 4.7532e-15 in the Terok GNU Debug build. The test also verifies
charge/metric equality, reciprocal conjugation, projector symmetry, adjointness, full/reduced
reciprocal-list agreement and a constant-potential shift. Both one- and two-rank MPI runs pass.
These checks establish the kernel; they are not end-to-end USPP energies or a performance benchmark.

Next integration steps are the generalized atomic initial guess, projector density and metric
assembly, augmentation contributions to the local effective potential and total energy, and their
force/stress derivatives. Direct quadrature is retained as a correctness reference when introducing
reciprocal interpolation and distributed grid operations.

## Hamiltonian and local-potential accounting

For the ordinary USPP functional the augmentation contribution to the screened projector matrix
contains the **total** local effective potential. It includes the local ionic potential as well as
Hartree and XC. Including only Hartree/XC would omit an energy term and change the pseudopotential.

CP2K splits the local ionic potential into a Gaussian-core Coulomb part and a short-range
complement. Its existing native AO short-range integral remains in the core Hamiltonian. The
augmented density additionally requires the integral of that short-range potential against Q;
interactions with short-range potentials on other atoms must also be retained. The Gaussian-core
part is already included when Hartree electrostatics is evaluated on the augmented valence density
plus the signed Gaussian core density. Its self-energy and ion-ion corrections must follow the
existing NC bookkeeping.

The implementation must retain the ordinary Gaussian overlap S_G separately from the generalized
metric S. Audit every consumer of the overlap rather than replacing an array globally. The SCF
eigensolver, occupations, orbital orthonormalization and electron-count checks need S; AO integrals
and the unaugmented density still refer to the original Gaussian basis.

Atomic initial guesses also require attention: the present SGP conversion returns no USPP operators.
Enabling the parser without replacing that path would produce an invalid atomic reference. A radial
generalized atomic solver and its augmentation density are part of the native implementation.

## Derivatives and subsequent decompositions

Forces must include derivatives of B, the augmentation translation phase and the generalized-overlap
Pulay term. The latter involves the energy-weighted density matrix. Stress additionally
differentiates reciprocal vectors, volume normalization, radial transforms and the
Gaussian/projector integrals. Check forces and stress against the discretized total energy before
accepting a mode.

GAPW_XC introduces a specified hard/soft decomposition in XC while keeping GPW electrostatics. GAPW
subsequently extends that decomposition to electrostatics. Specify independently whether the split
acts on the Gaussian orbital density, the USPP augmentation, or both. Nonoverlap and
environmental-density assumptions cannot be inferred from the word "one-center". The PRB manuscript
records the Coulomb cross terms and the XC decomposition residual that must be evaluated or bounded.
An exact reference for overlapping corrections uses a partition of unity and the full environmental
density; geometry derivatives of that partition must then be included too.

PAW adds AE/pseudo partial waves, projector duality, compensation multipoles, frozen-core
conventions and AE-minus-pseudo one-center energies. Its regular PAW one-center corrections are
required already in PAW GPW. Additional GAPW density splitting is a separate implementation step.
