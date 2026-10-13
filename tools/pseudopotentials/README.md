# SSSP downloads and UPF data preparation

Keep original downloaded collections separate from the prepared data shipped with CP2K. `data/UPF/`
currently contains one supplemented H USPP example derived from SSSP v2.0 PBE Precision, **not a
complete or unchanged SSSP collection**. Its `SOURCE` file records the original and supplemented
checksums, generator replay, attribution and license.

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

The optional `--install` mode instead places unchanged original files under
`DIRECTORY/{pbe,pbesol}/{efficiency,precision}/` and removes newly downloaded archives and temporary
extraction files. Its destination is also explicitly chosen by the user and must be writable. This
mode does **not** prepare partial waves or turn an unsupported UPF into a supported one. Existing
destination files with different contents are rejected.

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

Keep new generated data in a separate working directory until their provenance and CP2K calculation
checks have been reviewed. Only qualified prepared data should be added under `data/UPF/`, with
distinct filenames and their preparation recorded in `SOURCE`. Original collections remain upstream
inputs and must not overwrite supplemented data.

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

The common rescaling of the selected core orbitals reproduces the supplied nonnegative NLCC density.
It defines `tau_model = tau_W[n_model] + n_model * (tau_AE - tau_W[n_AE]) / n_AE`, where
`tau_W[n] = |grad(n)|^2/(8*n)`. The Pauli remainder is evaluated as positive orbital squares. This
is an explicitly chosen auxiliary model, not recovery of an original generator kinetic field. The
core orbitals must be normalized closed shells with occupations summing to the frozen-core charge.
Occupied valence states, negative NLCC densities, missing data and existing core tau are rejected.
All original sections are retained apart from the kinetic-data availability flag. The added
`PP_TAUMOD` declares the canonical spin-summed `hartree_per_volume` convention and records the model
and source provenance. Valence kinetic augmentation still requires full partial waves.

## Use the supplied H example

For a CP2K version with the native inverse-USPP implementation, the supplied PBE H dataset can be
selected using the standard CP2K data-directory lookup:

```text
&KIND H
  POTENTIAL UPF UPF/H.us.pbe.z_1.ld1.psl.v1.0.0-high-full-waves.upf
&END KIND
```

The initial inverse-USPP implementation supports fully periodic GPW calculations at Gamma and
requires matched full partial waves. Installing the dataset alone does not enable that
implementation or its unsupported methods and properties. NC, USPP and PAW datasets must not be
treated as interchangeable.

Choose and converge a compatible Gaussian basis, grid and response settings for the intended
calculation. SSSP plane-wave cutoffs are not GPW multigrid or Gaussian-basis recommendations. The
deliberately small H2 regression checks implementation consistency, not production convergence or
the accuracy of the complete SSSP collections.
