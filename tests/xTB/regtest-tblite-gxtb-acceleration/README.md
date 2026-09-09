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
