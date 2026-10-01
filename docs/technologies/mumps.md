# MUMPS

Enable the optional MPI backend with `CP2K_USE_MUMPS=ON` and `MUMPS_ROOT`, or provide
`MUMPS_INCLUDE_DIR` and `MUMPS_LIBRARIES`. Static linkage must include ordering libraries. Match
CP2K's MPI, Fortran ABI and integer/BLAS conventions. The tested version is MUMPS 5.9.1, distributed
under CeCILL-C.

`SPECTRAL_LOCALIZER / SOLVER MUMPS` supplies class-A inertia and generalized gap brackets, including
complex Hamiltonians through realification. It does not provide a Pfaffian sign.
`QUADRATIC_PSEUDOSPECTRUM / SOLVER ITERATIVE` reuses metric factors with DBCSR operator products and
a bounded Ritz space. `MAX_SUBSPACE` and `MAX_ITER` control that iteration.

Assembly and factorizations are distributed. Quadratic-response vectors and metric right-hand sides
use replicated/centralized storage. Fill-in and root fronts still limit scaling. Unresolved gaps or
invalid metrics produce no index. Dense reference paths remain available without this dependency.
