#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check SIRIUS Skala input rejection through the CP2K executable, without SCF."""

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess

from create_sirius_lapw_test_atom import helium_setup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cp2k", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--launcher", default="", help="For example, 'mpirun -n 2'.")
    args = parser.parse_args()
    executable = args.cp2k.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    source = (root / "tests/SIRIUS/regtest-skala/He-lapw-atom-grid.inp").read_text()
    cases = [
        (
            "missing-atom-file",
            source.replace("      SIRIUS_ATOM_FILE He-lapw.json\n", ""),
            "requires SIRIUS_ATOM_FILE",
        ),
        (
            "pseudopotential-method",
            source.replace("full_potential_lapwlo", "pseudopotential"),
            "SIRIUS_ATOM_FILE requires ELECTRONIC_STRUCTURE_METHOD",
        ),
        (
            "pseudopotential-kind",
            source.replace("POTENTIAL ALL", "POTENTIAL GTH-PBE-q2"),
            "not a pseudopotential",
        ),
        (
            "gaussian-basis",
            source.replace(
                "POTENTIAL ALL", "POTENTIAL ALL\n      BASIS_SET DZVP-MOLOPT-SR-GTH"
            ),
            "not BASIS_SET or BASIS",
        ),
        (
            "missing-atom-grid",
            re.sub(
                r"    &SKALA_GRID\n.*?    &END SKALA_GRID\n", "", source, flags=re.S
            ),
            "requires PW_DFT/SKALA_GRID",
        ),
        ("wrong-element", source, "does not match the requested element"),
        (
            "invalid-spin-angular",
            source.replace(
                "      N_ANGULAR 110", "      N_ANGULAR 110\n      N_SPIN_ANGULAR 8"
            ),
            "N_SPIN_ANGULAR requires an available positive Lebedev rule",
        ),
        (
            "vector-without-spinors",
            source.replace(
                "      POTENTIAL ALL",
                "      POTENTIAL ALL\n      SIRIUS_MAGNETIZATION_VECTOR 0.3 -0.2 0.1",
            ),
            "SIRIUS_MAGNETIZATION_VECTOR requires NUM_MAG_DIMS 3",
        ),
        (
            "vector-without-sirius",
            (root / "tests/SIRIUS/regtest-skala/H2.inp")
            .read_text()
            .replace("METHOD SIRIUS", "METHOD Quickstep")
            .replace(
                "      POTENTIAL GTH-PBE-q1",
                "      POTENTIAL GTH-PBE-q1\n      SIRIUS_MAGNETIZATION_VECTOR 0.3 -0.2 0.1",
            ),
            "SIRIUS_MAGNETIZATION_VECTOR requires METHOD SIRIUS",
        ),
        (
            "missing-force-grid",
            (root / "tests/SIRIUS/regtest-skala/H2.inp")
            .read_text()
            .replace("RUN_TYPE ENERGY", "RUN_TYPE ENERGY_FORCE"),
            "Native Skala forces and stress with SIRIUS require PW_DFT/SKALA_GRID",
        ),
    ]
    for name, text, expected in cases:
        directory = output / name
        directory.mkdir()
        (directory / "input.inp").write_text(text)
        atom = helium_setup()
        if name == "wrong-element":
            atom["symbol"] = "Ne"
            atom["number"] = 10
        (directory / "He-lapw.json").write_text(json.dumps(atom))
        command = shlex.split(args.launcher) + [str(executable), "-i", "input.inp"]
        with (directory / "output.log").open("w") as log:
            process = subprocess.Popen(
                command,
                cwd=directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                result = process.wait(timeout=90)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise RuntimeError(
                    f"{name}: input rejection hung; see {directory}"
                ) from None
        # CP2K wraps abort messages inside an ASCII-art frame.
        message = " ".join(
            " ".join(line.strip(" *|").split())
            for line in (directory / "output.log").read_text().splitlines()
        )
        if result == 0 or expected not in message:
            raise RuntimeError(
                f"{name}: expected rejection was not observed; see {directory}"
            )
        if "[find] iteration :" in message:
            raise RuntimeError(f"{name}: rejection must precede SCF; see {directory}")
        print(f"{name}: rejected as expected", flush=True)


if __name__ == "__main__":
    main()
