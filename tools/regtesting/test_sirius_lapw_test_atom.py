#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check the small LAPW setup generators without CP2K or a model runtime."""

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from create_sirius_lapw_test_atom import (
    helium_setup,
    hydrogen_setup,
    lithium_setup,
    neon_setup,
)


class LapwTestAtom(unittest.TestCase):
    def test_electron_integrals(self):
        for factory in (hydrogen_setup, helium_setup, lithium_setup, neon_setup):
            atom = factory()
            with self.subTest(element=atom["symbol"]):
                radius = atom["free_atom"]["radial_grid"]
                density = atom["free_atom"]["density"]
                self.assertEqual(len(radius), len(density))
                self.assertTrue(all(b > a > 0 for a, b in zip(radius, radius[1:])))
                self.assertTrue(all(math.isfinite(d) and d >= 0 for d in density))
                radial = [4 * math.pi * r * r * d for r, d in zip(radius, density)]
                electrons = math.fsum(
                    (b - a) * (u + v) / 2
                    for a, b, u, v in zip(radius, radius[1:], radial, radial[1:])
                )
                self.assertAlmostEqual(electrons, atom["number"], delta=1e-5)

    def test_core_and_valence(self):
        for factory, number, angular in (
            (lithium_setup, 3, {0}),
            (neon_setup, 10, {0, 1}),
        ):
            atom = factory()
            with self.subTest(element=atom["symbol"]):
                self.assertEqual(atom["core"], "1s")
                self.assertEqual(atom["number"], number)
                self.assertTrue(all(c["n"] >= 2 for c in atom["valence"] if "n" in c))
                self.assertEqual({orbital["l"] for orbital in atom["lo"]}, angular)
                self.assertTrue(
                    all(c["n"] == 2 for orbital in atom["lo"] for c in orbital["basis"])
                )

    def test_independent_setups(self):
        before = helium_setup()
        neon = neon_setup()
        neon["valence"][0]["basis"][0]["enu"] = 0
        self.assertEqual(helium_setup(), before)
        self.assertEqual(neon_setup()["valence"][0]["basis"][0]["enu"], -0.3)

    def test_cli(self):
        script = Path(__file__).with_name("create_sirius_lapw_test_atom.py")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "atom.json"
            for args, factory in (
                ([], helium_setup),
                (["--element", "H"], hydrogen_setup),
                (["--element", "Li"], lithium_setup),
                (["--element", "Ne"], neon_setup),
            ):
                with self.subTest(args=args):
                    subprocess.run(
                        [sys.executable, str(script), str(output), *args], check=True
                    )
                    self.assertEqual(json.loads(output.read_text()), factory())


if __name__ == "__main__":
    unittest.main()
