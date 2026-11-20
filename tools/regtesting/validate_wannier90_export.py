#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run and validate CP2K-generated Wannier90 exports.

Run external Wannier90, compare final spreads/centres, or compare exported
``.eig``, ``.mmn`` and optional rectangular ``.amn`` matrices. Only the export
comparison requires NumPy. It does not modify projections or CP2K inputs.

Raw matrices and AMN column-subspace comparisons assume the same MO band gauge.
MMN singular values are invariant under band rotations, but matching them alone
does not establish a consistent gauge over the mesh.
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

WF_RE = re.compile(
    r"WF centre and spread\s+(\d+)\s+\(\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?),\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?),\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*\)\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
OMEGA_RE = re.compile(
    r"Omega\s+(I|D|OD|Total)\s*=\s*" r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
ERROR_RE = re.compile(r"(^|\s)(ERROR:|Error in|Exiting\.\.\.\.\.\.\.)")


@dataclass(frozen=True)
class WannierFunction:
    index: int
    center: tuple[float, float, float]
    spread: float


@dataclass(frozen=True)
class WoutSummary:
    seed: str
    path: Path
    functions: tuple[WannierFunction, ...]
    omega_i: float
    omega_d: float
    omega_od: float
    omega_total: float
    all_done: bool


def parse_wout(seed: str, workdir: Path) -> WoutSummary:
    path = workdir / f"{seed}.wout"
    text = path.read_text(encoding="utf-8", errors="replace")
    if ERROR_RE.search(text):
        raise ValueError(f"{path}: Wannier90 reported an error")
    final_pos = text.rfind("Final State")
    if final_pos < 0:
        raise ValueError(f"{path}: missing final Wannier90 state")
    final_text = text[final_pos:]

    functions: list[WannierFunction] = []
    for match in WF_RE.finditer(final_text):
        functions.append(
            WannierFunction(
                index=int(match.group(1)),
                center=(
                    float(match.group(2)),
                    float(match.group(3)),
                    float(match.group(4)),
                ),
                spread=float(match.group(5)),
            )
        )
    if not functions:
        raise ValueError(f"{path}: missing final Wannier-function centres")

    omegas: dict[str, float] = {}
    for match in OMEGA_RE.finditer(final_text):
        omegas[match.group(1)] = float(match.group(2))
    required = {"I", "D", "OD", "Total"}
    missing = required - omegas.keys()
    if missing:
        raise ValueError(
            f"{path}: missing omega field(s): {', '.join(sorted(missing))}"
        )

    return WoutSummary(
        seed=seed,
        path=path,
        functions=tuple(functions),
        omega_i=omegas["I"],
        omega_d=omegas["D"],
        omega_od=omegas["OD"],
        omega_total=omegas["Total"],
        all_done="All done: wannier90 exiting" in text,
    )


def run_wannier90(executable: str, seed: str, workdir: Path, timeout: float) -> None:
    result = subprocess.run(
        [executable, seed],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )
    (workdir / f"{seed}.wannier90.stdout").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"{executable} {seed} failed with exit code {result.returncode}"
        )


def print_summary(summary: WoutSummary) -> None:
    done = "yes" if summary.all_done else "no"
    print(f"seed = {summary.seed}")
    print(f"wout = {summary.path}")
    print(f"all_done = {done}")
    print(f"num_wann = {len(summary.functions)}")
    print(f"omega_i = {summary.omega_i:.12e}")
    print(f"omega_d = {summary.omega_d:.12e}")
    print(f"omega_od = {summary.omega_od:.12e}")
    print(f"omega_total = {summary.omega_total:.12e}")
    for function in summary.functions:
        cx, cy, cz = function.center
        print(
            f"wf {function.index:5d}: center = "
            f"({cx:.12e}, {cy:.12e}, {cz:.12e}), "
            f"spread = {function.spread:.12e}"
        )


def max_abs_delta(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError(f"length mismatch: {len(left)} != {len(right)}")
    return max((abs(x - y) for x, y in zip(left, right)), default=0.0)


def compare_summaries(
    left: WoutSummary, right: WoutSummary
) -> tuple[float, float, float]:
    if len(left.functions) != len(right.functions):
        raise ValueError(
            f"number of Wannier functions differs: "
            f"{len(left.functions)} != {len(right.functions)}"
        )
    left_spreads = [function.spread for function in left.functions]
    right_spreads = [function.spread for function in right.functions]
    spread_delta = max_abs_delta(left_spreads, right_spreads)

    left_centers = [value for function in left.functions for value in function.center]
    right_centers = [value for function in right.functions for value in function.center]
    center_delta = max_abs_delta(left_centers, right_centers)
    omega_delta = abs(left.omega_total - right.omega_total)
    return omega_delta, spread_delta, center_delta


def read_eig(path: Path, nbands: int, nkpts: int) -> list[float]:
    values: dict[tuple[int, int], float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{path}: malformed eigenvalue entry")
        band, kpoint = map(int, fields[:2])
        value = float(fields[2])
        if not (1 <= band <= nbands and 1 <= kpoint <= nkpts):
            raise ValueError(f"{path}: eigenvalue index out of range")
        if (kpoint, band) in values or not math.isfinite(value):
            raise ValueError(f"{path}: duplicate or nonfinite eigenvalue")
        values[kpoint, band] = value
    if len(values) != nbands * nkpts:
        raise ValueError(f"{path}: missing eigenvalues")
    return [values[index] for index in sorted(values)]


def read_amn(path: Path, nbands: int, nkpts: int) -> list[NDArray[np.complex128]]:
    import numpy as np

    with path.open(encoding="utf-8") as stream:
        stream.readline()  # Wannier90 comment line.
        dimensions = tuple(map(int, stream.readline().split()))
        if len(dimensions) != 3 or dimensions[:2] != (nbands, nkpts):
            raise ValueError(f"{path}: AMN dimensions disagree with MMN dimensions")
        nproj = dimensions[2]
        if not 0 < nproj <= nbands:
            raise ValueError(f"{path}: invalid AMN projection count")
        matrices = np.full((nkpts, nbands, nproj), np.nan, dtype=np.complex128)
        for line in stream:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"{path}: malformed AMN entry")
            band, projection, kpoint = map(int, fields[:3])
            if not (
                1 <= band <= nbands
                and 1 <= projection <= nproj
                and 1 <= kpoint <= nkpts
            ):
                raise ValueError(f"{path}: AMN index out of range")
            value = complex(float(fields[3]), float(fields[4]))
            if not np.isfinite(value) or not np.isnan(
                matrices[kpoint - 1, band - 1, projection - 1]
            ):
                raise ValueError(f"{path}: duplicate or nonfinite AMN entry")
            matrices[kpoint - 1, band - 1, projection - 1] = value
    if not np.all(np.isfinite(matrices)):
        raise ValueError(f"{path}: missing AMN entries")
    return [matrices[kpoint] for kpoint in range(nkpts)]


def compare_amn(
    left: Sequence[NDArray[np.complex128]],
    right: Sequence[NDArray[np.complex128]],
) -> dict[str, float]:
    """Compare projection columns in a common MO band gauge."""
    import numpy as np

    if not left or len(left) != len(right):
        raise ValueError("incompatible AMN k-point counts")
    raw_max = singular_max = subspace_max = 0.0
    for kpoint, (a, b) in enumerate(zip(left, right), start=1):
        if a.ndim != 2 or a.shape != b.shape or not 0 < a.shape[1] <= a.shape[0]:
            raise ValueError(f"incompatible AMN matrices at k-point {kpoint}")
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise ValueError(f"nonfinite AMN matrix at k-point {kpoint}")
        ua, sa, _ = np.linalg.svd(a, full_matrices=False)
        ub, sb, _ = np.linalg.svd(b, full_matrices=False)
        if not np.all(np.isfinite(sa)) or not np.all(np.isfinite(sb)):
            raise ValueError(f"nonfinite AMN singular values at k-point {kpoint}")
        tolerance = np.finfo(float).eps * max(a.shape)
        if sa[-1] <= tolerance * sa[0] or sb[-1] <= tolerance * sb[0]:
            raise ValueError(f"rank-deficient AMN projections at k-point {kpoint}")
        overlap = ua.conj().T @ ub
        raw_max = max(raw_max, float(np.max(np.abs(a - b))))
        singular_max = max(singular_max, float(np.max(np.abs(sa - sb))))
        # Projector residuals stay accurate for nearly identical column spaces.
        subspace_max = max(
            subspace_max,
            float(np.linalg.norm(ub - ua @ overlap, ord=2)),
            float(np.linalg.norm(ua - ub @ overlap.conj().T, ord=2)),
        )
    return {
        "amn_max_abs_delta": raw_max,
        "amn_max_singular_delta": singular_max,
        "amn_subspace_deviation": subspace_max,
    }


def command_compare_export(args: argparse.Namespace) -> int:
    try:
        from compare_wannier90_mmn import compare_mmn, read_mmn
    except ImportError as exc:
        raise RuntimeError("compare-export requires NumPy") from exc

    tolerances = {
        "eig_max_abs_delta": args.eig_tol,
        "mmn_max_abs_delta": args.mmn_tol,
        "mmn_max_singular_delta": args.mmn_singular_tol,
        "amn_max_abs_delta": args.amn_tol,
        "amn_max_singular_delta": args.amn_singular_tol,
        "amn_subspace_deviation": args.amn_subspace_tol,
    }
    for name, tolerance in tolerances.items():
        if tolerance is not None and (not math.isfinite(tolerance) or tolerance < 0):
            raise ValueError(f"{name}: tolerance must be finite and nonnegative")
        if name.startswith("amn_") and tolerance is not None and not args.compare_amn:
            raise ValueError("AMN tolerances require --compare-amn")
    left = Path(args.left_workdir) / args.left_seed
    right = Path(args.right_workdir) / args.right_seed
    left_mmn = read_mmn(str(left) + ".mmn")
    right_mmn = read_mmn(str(right) + ".mmn")
    raw, singular, _, _ = compare_mmn(left_mmn, right_mmn)
    nbands, nkpts = left_mmn[:2]
    left_eig = read_eig(Path(str(left) + ".eig"), nbands, nkpts)
    right_eig = read_eig(Path(str(right) + ".eig"), nbands, nkpts)
    metrics = {
        "eig_max_abs_delta": max_abs_delta(left_eig, right_eig),
        "mmn_max_abs_delta": raw,
        "mmn_max_singular_delta": singular,
    }
    if args.compare_amn:
        metrics.update(
            compare_amn(
                read_amn(Path(str(left) + ".amn"), nbands, nkpts),
                read_amn(Path(str(right) + ".amn"), nbands, nkpts),
            )
        )
    failed = False
    for name, value in metrics.items():
        if not math.isfinite(value):
            raise ValueError(f"{name}: nonfinite comparison result")
        print(f"{name} = {value:.16e}")
        tolerance = tolerances[name]
        failed |= tolerance is not None and value > tolerance
    return int(failed)


def command_run(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir)
    if args.run:
        run_wannier90(args.executable, args.seed, workdir, args.timeout)
    summary = parse_wout(args.seed, workdir)
    print_summary(summary)
    if args.require_all_done and not summary.all_done:
        print("Wannier90 did not report normal completion.", file=sys.stderr)
        return 1
    if args.max_omega_total is not None and summary.omega_total > args.max_omega_total:
        print(
            f"Omega total {summary.omega_total:.12e} exceeds "
            f"{args.max_omega_total:.12e}.",
            file=sys.stderr,
        )
        return 1
    return 0


def command_compare(args: argparse.Namespace) -> int:
    left = parse_wout(args.left_seed, Path(args.left_workdir))
    right = parse_wout(args.right_seed, Path(args.right_workdir))
    omega_delta, spread_delta, center_delta = compare_summaries(left, right)
    print(f"omega_total_delta = {omega_delta:.12e}")
    print(f"spread_delta = {spread_delta:.12e}")
    print(f"center_delta = {center_delta:.12e}")

    failed = False
    if args.omega_total_tol is not None and omega_delta > args.omega_total_tol:
        failed = True
    if args.spread_tol is not None and spread_delta > args.spread_tol:
        failed = True
    if args.center_tol is not None and center_delta > args.center_tol:
        failed = True
    return 1 if failed else 0


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run", help="run and/or parse one Wannier90 seed"
    )
    run_parser.add_argument("seed", help="Wannier90 seed name")
    run_parser.add_argument(
        "--workdir", default=".", help="directory containing the seed files"
    )
    run_parser.add_argument(
        "--run", action="store_true", help="run wannier90.x before parsing"
    )
    run_parser.add_argument(
        "--executable", default="wannier90.x", help="Wannier90 executable"
    )
    run_parser.add_argument(
        "--timeout", type=float, default=300.0, help="execution timeout in seconds"
    )
    run_parser.add_argument("--max-omega-total", type=float, default=None)
    run_parser.add_argument("--require-all-done", action="store_true")
    run_parser.set_defaults(func=command_run)

    compare_parser = subparsers.add_parser(
        "compare", help="compare two parsed .wout files"
    )
    compare_parser.add_argument("left_seed")
    compare_parser.add_argument("right_seed")
    compare_parser.add_argument("--left-workdir", default=".")
    compare_parser.add_argument("--right-workdir", default=".")
    compare_parser.add_argument("--omega-total-tol", type=float, default=None)
    compare_parser.add_argument("--spread-tol", type=float, default=None)
    compare_parser.add_argument("--center-tol", type=float, default=None)
    compare_parser.set_defaults(func=command_compare)

    export_parser = subparsers.add_parser(
        "compare-export",
        help="compare EIG/MMN and optional AMN data (requires NumPy)",
        description=(
            "Compare exported matrices without modifying them. Raw matrices and "
            "AMN column subspaces require a common MO band gauge. Matching MMN "
            "singular values alone does not prove mesh-wide gauge equivalence."
        ),
    )
    export_parser.add_argument("left_seed")
    export_parser.add_argument("right_seed")
    export_parser.add_argument("--left-workdir", default=".")
    export_parser.add_argument("--right-workdir", default=".")
    export_parser.add_argument("--compare-amn", action="store_true")
    export_parser.add_argument(
        "--eig-tol",
        type=float,
        default=None,
        help="absolute eigenvalue tolerance in eV",
    )
    export_parser.add_argument("--mmn-tol", type=float, default=None)
    export_parser.add_argument("--mmn-singular-tol", type=float, default=None)
    export_parser.add_argument("--amn-tol", type=float, default=None)
    export_parser.add_argument("--amn-singular-tol", type=float, default=None)
    export_parser.add_argument(
        "--amn-subspace-tol",
        type=float,
        default=None,
        help="maximum column-projector residual in a common MO band gauge",
    )
    export_parser.set_defaults(func=command_compare_export)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, subprocess.TimeoutExpired, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
