# SSSP downloads and UPF data preparation

SSSP is a curated selection from several pseudopotential families, not a pseudopotential type. Its
collections contain norm-conserving, ultrasoft and PAW datasets in UPF format.

Keep downloaded collections outside the source and installed data trees. The supplemented H USPP
under `tests/QS/regtest-upf/` is an offline regression input originating from PSLibrary 1.0.0 and
selected by SSSP v2.0 PBE Precision. Its adjacent `SOURCE` file records the preparation and license.
It is not installed in the CP2K data directory.

## Install a selection with prepared USPPs

Choose a writable directory outside the source or installed data tree, for example:

```bash
python3 tools/pseudopotentials/fetch_sssp.py \
  "$HOME/.cache/cp2k/upf" --install --library pbe-prec
```

This installs unchanged NC and PAW files from SSSP. Each selected USPP is replaced by its
supplemented counterpart from the fixed
[inverse-uspp v0.3.0 prerelease](https://github.com/DCM-Uni-Paderborn/inverse-uspp/releases/tag/v0.3.0).
`inverse_uspp_release.json` pins its catalog and archive checksums. Matching the original UPF hash
preserves the SSSP functional, valence partition and Efficiency/Precision selection.

Original NC/PAW files and unchanged SSSP cutoff metadata reside under
`DIRECTORY/{pbe,pbesol}/{efficiency,precision}/`. Prepared USPPs are stored once under
`DIRECTORY/prepared-uspp/v0.3.0/`, together with their catalog, source information and license. Each
selection directory contains `selection.json`, which maps the original filenames to actual file
locations relative to that directory. In CP2K input, set `POTENTIAL UPF` and give the resolved path
in `POTENTIAL_FILE_NAME` in the same `&KIND` section. This keyword accepts paths longer than the
80-character limit of `POTENTIAL UPF <filename>`. Choose the matching basis separately.

Set `METHOD` in `&QS` to `GPW`, `GAPW_XC`, or `GAPW` to choose the density evaluation.
For USPP and PAW datasets with overlapping reconstruction or projector supports, set
`USPP_METHOD INVERSE`. The default `INVERSE_1C` requires disjoint supports. `USPP_METHOD DIRECT`
expands the soft orbitals for either family. In USPP LDA/GGA it uses the original projectors and
augmentation without full reference waves. PAW uses its valence partial waves and one-center
corrections in all three orbital representations.

The prepared collection stores the SSSP selection under `data/{pbe,pbesol}/sssp/`. Add
`--additional-uspp` to install non-SSSP PSLibrary and GBRV USPPs for the selected functionals. They
reside under `data/{pbe,pbesol}/additional/{pslibrary,gbrv}/` and do not replace the SSSP selection.
Other functionals in the external release are outside this SSSP installer's scope.

The prerelease distinguishes compatible generator replays (route 3) from explicitly auxiliary
reconstructions (route 4). Original potential, projector, augmentation and core fields are retained.
All 860 release files passed development-reader and atomic inverse checks. This does not establish
all system calculations or the reduced one-center identities. Auxiliary reconstructions require the
original functional's density, metric, kinetic-energy and geometry corrections.

The 76 SSSP USPPs with nonzero NLCC also supply an explicitly prescribed orbital-rescaling core
kinetic density. This auxiliary meta-GGA prescription retains the original potential and core
charge. It is not recovered generator core tau. The other ten SSSP USPPs need no such addition.
Their selected occupations and reference sources are recorded in the catalog. External core
references are installed under `inputs/core-kinetic/` for the selected files.

Add `--paw-core` to install compatible complete frozen-core reference waves for the selected La, Pt,
Pa and U PAW files from the fixed
[PAW core prerelease](https://github.com/DCM-Uni-Paderborn/inverse-uspp/releases/tag/paw-core-v0.1.0).
Its archives and catalog are pinned in `paw_core_release.json`. The selection index points matching
original hashes to the supplemented files under `DIRECTORY/prepared-paw-core/v0.1.0/`. Other PAW
files retain their original contents. Only the core-orbital table is supplemented. Original
potentials, augmentation, valence waves and core densities remain unchanged. This supplies the core
reference data required by tau-dependent XC without establishing every SCF or geometry combination.

Add `--nc-core` after `DIRECTORY` to install auxiliary NC core kinetic fields from the pinned
[SSSP v2.0 NC core prerelease](https://github.com/DCM-Uni-Paderborn/inverse-uspp/releases/tag/nc-core-v0.1.0).
It covers the 40 distinct NC files with nonzero NLCC. Supply `--nc-core CATALOG` to use a local
catalog instead. The installer matches each field to the SHA-256 of its original UPF and records the
catalog checksum in `selection.json`. The prepared NC files reside under
`DIRECTORY/prepared-nc-core/`. The positive-amplitude model preserves the original signed core
charge and operator tables. The installer retains the exact catalog as
`DIRECTORY/prepared-nc-core/catalog.json.gz` with the added fields and their provenance.

Omit `--library` and combine both core options to install all four SSSP v2 selections with
prepared USPPs and the PAW and NC core supplements.

```bash
python3 tools/pseudopotentials/fetch_sssp.py \
  "$HOME/.cache/cp2k/upf-complete" --install --paw-core --nc-core
```

## Download original collections

From the CP2K repository root, choose a writable directory outside the source or installed data
tree, for example:

```bash
python3 tools/pseudopotentials/fetch_sssp.py \
  "$HOME/.cache/cp2k/sssp/2.0" --library pbe-prec
```

This verifies the archive against `sssp_v2_archives.json`, retains it as a cache, and extracts the
original archive contents below `DIRECTORY/pbe-prec/`, retaining the upstream directory layout. The
PBE Precision UPFs are in `pbe-prec/mix-sssp-prec-pbe-lib-v2/library/`. Select `pbe-eff`,
`pbe-prec`, `pbesol-eff` or `pbesol-prec`; repeat `--library` for several collections or omit it for
all four.

Use `--install --uspp-source original` to install all unchanged original UPFs under
`DIRECTORY/{pbe,pbesol}/{efficiency,precision}/`, together with the original cutoff metadata and a
selection index. Use a separate destination when changing the source mode or prepared release.
Installation removes newly downloaded archives and temporary files. Existing files with different
contents are rejected. The prepared mode still downloads complete source archives before selecting
their NC/PAW files, so filtering the installation does not reduce that archive transfer.

Neither mode runs automatically during a calculation, configuration, build or installation. Do not
redirect `CP2K_DATA_DIR` to a cache containing only UPFs: the standard basis and potential files
must remain available. Use explicit paths when referring to external original UPFs.

## Prepare missing USPP partial waves

`add_upf_full_waves.py` handles scalar ultrasoft UPF 2.0.1 files without spin-orbit or PAW flags. It
is a verifier and supplementer, **not an automatic converter for every SSSP potential**.

1. Obtain the original potential's recipe and compatible generator version.
1. Replay the generator while requesting projector-indexed AE and PS partial waves. For the supplied
   H example, use Quantum ESPRESSO `qe-6.3` and `lsave_wfc=.true.` in `&inputp`.
1. Compare the replay with the original and add the verified full-wave block:

```bash
python3 tools/pseudopotentials/add_upf_full_waves.py \
  "$HOME/.cache/cp2k/sssp/2.0/pbe-prec/mix-sssp-prec-pbe-lib-v2/library/H.us.pbe.z_1.ld1.psl.v1.0.0-high.upf" \
  /path/to/generated-H.upf \
  /path/to/prepared/H.us.pbe.z_1.ld1.psl.v1.0.0-high-full-waves.upf
```

The prepared output's parent directory must already exist and the output file must not exist. The
tool does not run the external generator. It requires agreement of the original mesh, quadrature,
potential, projectors and augmentation with the replay; only the atomic guess waves and density
allow roundoff differences. Original bytes are retained apart from the `has_wfc` flag, and
`PP_FULL_WFC` is added. Missing core kinetic data are not inferred.

For small numerical differences from a generator replay, `--replay-tolerance VALUE` permits a
relative maximum difference per numeric table and per reference-energy or cutoff attribute. The
radial mesh, quadrature, functional, electron count and projector assignment must still agree
exactly. All original arrays remain unchanged. Check the reconstructed calculation's internal
consistency before using such supplemented data. The default still requires exact operator tables.

Keep generated data in a separate working directory. Preserve the original collections and record
the preparation, checksums and licensing alongside supplemented files.

## Supplement missing PAW core orbitals

`add_upf_paw_core.py` replaces an incomplete scalar PAW GIPAW core-orbital table with verified
generator core orbitals. The generator must export the actual selected core states, which need not
be the first entries of its atomic orbital list.

```bash
python3 tools/pseudopotentials/add_upf_paw_core.py \
  /path/to/original-PAW.upf /path/to/generated-complete-core.upf \
  /path/to/prepared/PAW-complete-core.upf --core-tolerance 1e-9
```

The radial mesh and quadrature must agree exactly. The explicit tolerance applies to the AE core
density, shared signed core orbitals up to an overall phase, and reconstruction of the original AE
core density from normalized closed shells. The core occupations must sum to the original frozen
electron count and exclude occupied valence states. Only the core-orbital table is replaced. All
potential, projector, augmentation and valence-wave bytes remain unchanged. This checks the core
reference independently of regeneration of the valence pseudopotential. It does not infer missing
USPP model-core kinetic data or establish complete SCF and geometry support.

## Prepare an explicit USPP model core kinetic density

For tau-dependent XC with NLCC, `add_upf_core_tau.py` can define an auxiliary core model from
original GIPAW core orbitals. Select every core-shell occupation explicitly in the original GIPAW
index order. For the Al-q3 and P-q5 examples the frozen Ne core gives `2 2 6 0`, excluding the
valence 3s orbital also present in the GIPAW table.

```bash
python3 tools/pseudopotentials/add_upf_core_tau.py \
  /path/to/prepared/Al-full-waves.upf /path/to/prepared/Al-full-waves-core-tau.upf \
  --model orbital-rescale --core-occupations 2 2 6 0
```

For GBRV data, `--vanderbilt-ae /path/to/atom.ae` instead reads the retained Vanderbilt all-electron
state. Supply one occupation for every orbital in that file, assigning zero to each valence orbital.
The reader supports nonrelativistic and scalar-relativistic states with four-byte Fortran record
markers and eight-byte real values in either byte order. It checks the nuclear charge, closed-shell
occupations and normalization before and after natural cubic interpolation onto the unchanged UPF
radial mesh. The atomic-state hash is recorded in the output. A new generator run is unnecessary
when this state is already available.

An incomplete GIPAW table can instead be completed with `--qe-ae-wfc /path/to/ld1.wfc`. The radial
nodes must agree to the QE text format's rounding. Shared orbitals must agree up to phase within a
relative radial L2 difference of `1e-6`. Keep occupations in GIPAW index order, then append those
for missing QE columns in header order, including zeros for unoccupied states.

The common rescaling of the selected core orbitals reproduces the supplied nonnegative NLCC density.
It defines `tau_model = tau_W[n_model] + n_model * (tau_AE - tau_W[n_AE]) / n_AE`, where
`tau_W[n] = |grad(n)|^2/(8*n)`. The Pauli remainder is evaluated as positive orbital squares. This
is an explicitly chosen auxiliary model, not recovery of an original generator kinetic field. The
core orbitals must be normalized closed shells with occupations summing to the frozen-core charge.
Occupied valence states, negative NLCC densities, missing data and existing core tau are rejected.
All original sections are retained apart from the kinetic-data availability flag. The added
`PP_TAUMOD` declares the canonical spin-summed `hartree_per_volume` convention and records the model
and source provenance. Valence kinetic augmentation still requires full partial waves.

For scalar ONCV norm-conserving UPFs with signed NLCC samples, select the positive-amplitude model
explicitly. It derives the kinetic field from the positive amplitude of the NLCC table and preserves
the original signed charge. The tool checks the closed-shell occupations against the embedded atomic
input and checks normalization on the native reference mesh.

```bash
python3 tools/pseudopotentials/add_upf_core_tau.py original-NC.upf prepared-NC.upf \
  --model positive-amplitude --oncv-core core-waves.dat
```

The export begins with `ONCV_AE_CORE_V1 element Z shell_count radial_count xc_code`. Each occupied
core shell begins with `SHELL n l occupation`, followed by `radial_count` rows of radius in bohr and
signed reduced radial wave `u(r)`. Supply compatible reference waves from the atomic generator. The
converter records the export and original file hashes.

## Use the supplied H example

For a CP2K version with the native inverse-USPP implementation, select the regression dataset using
an explicit path to the source tree:

```text
&KIND H
  POTENTIAL UPF
  POTENTIAL_FILE_NAME /path/to/cp2k/tests/QS/regtest-upf/H.us.pbe.z_1.ld1.psl.v1.0.0-high-full-waves.upf
&END KIND
```

Calculation support depends on the CP2K build and Gaussian basis. Inverse USPP and tau-dependent
XC with any USPP orbital representation require complete projector-indexed AE/PS waves. Inverse
USPP retains the original density, metric and kinetic-energy corrections and their geometry
derivatives. For tau-dependent XC, NC and USPP with NLCC also require a compatible model core
kinetic density. PAW uses its AE/PS valence waves and, where present, the occupied frozen-core
waves required to build the core kinetic density.

Choose and converge a compatible Gaussian basis, grid and response settings for the intended
calculation. SSSP plane-wave cutoffs are not GPW multigrid or Gaussian-basis recommendations. The
deliberately small H2 regression checks implementation consistency, not production convergence or
the accuracy of the complete SSSP collections.
