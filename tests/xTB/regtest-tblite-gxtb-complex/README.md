# Anisotropic finite-mesh symmetry and complex derivatives

The methane geometry is tetrahedral in a cubic primitive cell, but its 3x1x1
and 3x3x1 Monkhorst-Pack meshes do not preserve every crystal rotation.
These tests retain only operations preserving the complete finite mesh,
including the additional operations fixing Gamma. This is required for the
coupled g-xTB exchange: the BvK supercell has the symmetry of the mesh, not
automatically the full symmetry of the primitive crystal.

Full-grid, K290 and SPGLIB 3x1x1 calculations agree in central energy. The
3x3x1 full/SPGLIB pair additionally exercises a genuine spatial reduction
from nine to three k points, rather than just time reversal. Every case
checks analytic forces and stress against finite differences. Reduction
counts refer to the initial geometry; FD distortions lower its symmetry.

Before the whole-mesh subgroup correction, the SPGLIB 3x1x1 input aborted at
the first SCF Fock validation with residual 2.6931e-5. Do not replace this
input with an isotropic mesh or relax the covariance/FD tolerances.

The full-grid and reduced central energies should agree within 1e-9 Eh.
The DEBUG matchers check the derivative residuals against the same 1e-6
tolerances used by the other g-xTB regression directories.
