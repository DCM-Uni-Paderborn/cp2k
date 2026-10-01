# Tacho sparse Pfaffians

`CP2K_USE_TACHO=ON` enables `SPECTRAL_LOCALIZER / SOLVER TACHO` for two-dimensional class-AII/Z2
analysis. It also requires `CP2K_USE_MUMPS=ON` for independent metric and gap checks. Both options
are disabled by default.

The tested Tacho source is
[iyamazaki/Trilinos, tacho-sk](https://github.com/iyamazaki/Trilinos/tree/tacho-sk) at revision
`02ef047cb4e1e16959714ce50462bb09511d0e9a`. Build with `Trilinos_ENABLE_ShyLU_NodeTacho=ON`,
`Kokkos_ENABLE_SERIAL=ON`, `Tacho_ENABLE_INT_INT=ON`, BLAS/LAPACK and MATLAB/tests/examples
disabled. Provide `ShyLU_NodeTacho_DIR` and Boost's CMake configuration. Match CP2K's C++ ABI and
BLAS integer width. Tacho is BSD-2-Clause licensed; Boost Graph is header-only and uses the Boost
Software License. No MATLAB is needed.

The CP2K adapter supplies skew-specific delayed-pivot numerical factorization using Tacho's symbolic
analysis and BLAS interfaces. It retains permutation signs and uses no pivot perturbation. This is
not a generic standard-Trilinos Pfaffian interface, and other versions have not been validated.

DBCSR assembly and MUMPS checks remain distributed, but sparse Pfaffian data and numerical factors
reside on one source rank. `PFAFFIAN_MEMORY` limits numerical factor/front/update/solve buffers in
MiB, not whole-process memory. It excludes symbolic ordering, graph matching, sparse inputs and
MUMPS. Unresolved gaps, failed residual checks and exhausted budgets produce no index.
Three-dimensional AII calculations continue to require `SOLVER DENSE`.
