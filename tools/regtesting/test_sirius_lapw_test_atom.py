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

from create_sirius_lapw_test_atom import helium_setup, neon_setup


class LapwTestAtom(unittest.TestCase):
    def test_electron_integral(self):
        for factory, core in ((helium_setup, ""), (neon_setup, "1s")):
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
                self.assertEqual(atom["core"], core)

    def test_cli(self):
        script = Path(__file__).with_name("create_sirius_lapw_test_atom.py")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "atom.json"
            for args, factory in (
                ([], helium_setup),
                (["--element", "Ne"], neon_setup),
            ):
                with self.subTest(args=args):
                    subprocess.run(
                        [sys.executable, str(script), str(output), *args], check=True
                    )
                    self.assertEqual(json.loads(output.read_text()), factory())


if __name__ == "__main__":
    unittest.main()
