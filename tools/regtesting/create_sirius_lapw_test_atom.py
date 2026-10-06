#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Write small LAPW integration fixtures, not production atom setups."""

import argparse
import json
import math
from pathlib import Path


def helium_setup():
    radius = [1e-7 + (20.0 - 1e-7) * (i / 2999) ** 3 for i in range(3000)]
    basis = [{"enu": -0.3, "dme": dme, "auto": 0} for dme in (0, 1)]
    return {
        "name": "Helium integration fixture",
        "symbol": "He",
        "number": 2,
        "mass": 4.002602,
        "rmin": 1e-6,
        "rmt": 1.8,
        "nrmt": 801,
        "core": "",
        "valence": [{"basis": basis}]
        + [{"n": l + 1, "l": l, "basis": basis} for l in range(4)],
        "lo": [
            {
                "l": 0,
                "basis": [
                    {"n": 1, "enu": -0.7, "dme": dme, "auto": 0} for dme in (0, 1)
                ],
            }
        ],
        "free_atom": {
            "radial_grid": radius,
            "density": [16.0 / math.pi * math.exp(-4.0 * r) for r in radius],
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    atom = helium_setup()
    args.output.write_text(json.dumps(atom, indent=2) + "\n")
