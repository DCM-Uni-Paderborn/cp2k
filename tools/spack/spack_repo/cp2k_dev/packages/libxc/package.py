# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.libxc.package import (
    AutotoolsBuilder as BuiltinAutotoolsBuilder,
    CMakeBuilder as BuiltinCMakeBuilder,
    Libxc as BuiltinLibxc,
)

from spack.package import *


class Libxc(BuiltinLibxc):
    # Share the toolchain's derivative corrections without copying the patch.
    patch(
        "https://raw.githubusercontent.com/DCM-Uni-Paderborn/cp2k/4eb521b0bde9980bd2e298389b59f1bf79919eb7/tools/toolchain/scripts/stage3/libxc-7.1.2-derivatives.patch",
        sha256="5d1797565de91bb9c20312cdbf492c5bbbdc701be7daa29b5061f933cd2109b0",
        when="@7.1.2",
    )

    depends_on("python", type="build", when="@7.1.2")
    depends_on("py-sympy@1.13.1", type="build", when="@7.1.2")
    depends_on("py-mpmath@1.3.0", type="build", when="@7.1.2")

    def regenerate_r2scan(self):
        if not self.spec.satisfies("@7.1.2"):
            return
        maxorder = 4 if "+lxc" in self.spec else 3 if "+kxc" in self.spec else 2
        # Preserve the generator's functional-specific differentiability limits.
        with working_dir(self.stage.source_path):
            for functional in ("mgga_x_r2scan", "mgga_c_r2scan"):
                self.spec["python"].command(
                    "-s", "scripts/sympy2c/compile.py", functional, "--maxorder", str(maxorder)
                )


class CMakeBuilder(BuiltinCMakeBuilder):
    @run_before("cmake")
    def regenerate_r2scan(self):
        self.pkg.regenerate_r2scan()


class AutotoolsBuilder(BuiltinAutotoolsBuilder):
    @run_before("configure")
    def regenerate_r2scan(self):
        self.pkg.regenerate_r2scan()
