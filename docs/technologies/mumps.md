# MUMPS

Enable the optional MPI backend with `CP2K_USE_MUMPS=ON` and `MUMPS_ROOT`, or provide
`MUMPS_INCLUDE_DIR` and `MUMPS_LIBRARIES`. Static linkage must include ordering libraries. Match
CP2K's MPI, Fortran ABI and integer/BLAS conventions. The tested version is MUMPS 5.9.1, distributed
under CeCILL-C.

`SPECTRAL_LOCALIZER / SOLVER MUMPS` supplies class-A inertia and generalized gap brackets, including
complex Hamiltonians through realification. It does not provide a Pfaffian sign.
`SPECTRAL_LOCALIZER / SOLVER DBCSR` and `QUADRATIC_PSEUDOSPECTRUM / SOLVER ITERATIVE` use native
pivoted LDL factors and do not require MUMPS.

Assembly and factorizations are distributed. Quadratic-response vectors and metric right-hand sides
use replicated/centralized storage. Fill-in and root fronts still limit scaling. Unresolved gaps or
invalid metrics produce no index. Dense reference methods also remain available without this
dependency.
