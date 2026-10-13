#!/usr/bin/env python3
"""Add an explicitly chosen orbital-rescaling NLCC kinetic model to a USPP.

The original GIPAW radial core orbitals and user-specified closed-shell
occupations define the reference core. A common radial rescaling reproduces
the supplied nonnegative NLCC charge. This defines auxiliary model orbitals,
not recovered generator core tau or normalized atomic core eigenstates.
Original sections are retained, apart from the kinetic-data availability flag.
"""

import argparse
from hashlib import sha256
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from add_upf_full_waves import flag, real, top_level_spans, values

ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba "
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb "
    "Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh "
    "Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def slopes(r, y, monotone=False):
    """Node derivatives of a natural cubic spline or shape-preserving Hermite.

    Signed reference orbitals use the natural cubic spline. The nonnegative
    model density uses shape-preserving Hermite slopes, including zero slopes
    at internal zero plateaus. No density or kinetic value is clipped.
    """
    n = len(r)
    h = [r[i + 1] - r[i] for i in range(n - 1)]
    d = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    if monotone:
        result = [0.0] * n
        for i in range(1, n - 1):
            if d[i - 1] != 0 and d[i] != 0 and (d[i - 1] > 0) == (d[i] > 0):
                a, b = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
                result[i] = (a + b) / (a / d[i - 1] + b / d[i])
        for index, first, second in [(0, 0, 1), (-1, -1, -2)]:
            v = ((2 * h[first] + h[second]) * d[first] - h[first] * d[second]) / (
                h[first] + h[second]
            )
            if v == 0 or d[first] == 0 or (v > 0) != (d[first] > 0):
                v = 0.0
            elif (d[second] == 0 or (d[first] > 0) != (d[second] > 0)) and abs(
                v
            ) > 3 * abs(d[first]):
                v = 3 * d[first]
            result[index] = v
        if r[0] == 0:
            result[0] = 0.0
        return result
    # Thomas solve for the natural spline's second derivatives.
    diagonal, rhs, curvature = [1.0] * n, [0.0] * n, [0.0] * n
    for i in range(1, n - 1):
        factor = h[i - 1] / diagonal[i - 1] if i > 1 else 0.0
        diagonal[i] = 2 * (h[i - 1] + h[i]) - factor * h[i - 1]
        rhs[i] = 6 * (d[i] - d[i - 1]) - factor * rhs[i - 1]
    for i in range(n - 2, 0, -1):
        curvature[i] = (rhs[i] - h[i] * curvature[i + 1]) / diagonal[i]
    result = [
        d[i] - h[i] * (2 * curvature[i] + curvature[i + 1]) / 6 for i in range(n - 1)
    ]
    result.append(d[-1] + h[-1] * (curvature[-2] + 2 * curvature[-1]) / 6)
    return result


def weights(rab):
    """The original UPF composite Simpson rule, with a final trapezoid if needed."""
    n = len(rab)
    last = n - 1 if n % 2 else n - 2
    w = [(2 if i % 2 == 0 else 4) / 3 if i <= last else 0.0 for i in range(n)]
    w[0] = w[last] = 1 / 3
    if last != n - 1:
        w[-2] += 0.5
        w[-1] += 0.5
    return [a * b for a, b in zip(w, rab)]


def integer(value):
    v = real(value)
    if v != int(v):
        raise ValueError("Noninteger orbital quantum number or core electron count")
    return int(v)


def reference_orbitals(tree, r, rab, occupations):
    header = tree.find("PP_HEADER")
    table = tree.find("PP_GIPAW/PP_GIPAW_CORE_ORBITALS")
    if table is None:
        raise ValueError("Missing original GIPAW core orbital table")
    count = int(table.get("number_of_core_orbitals", "-1"))
    if len(occupations) != count or len(table) != count:
        raise ValueError("Supply one core occupation for every original GIPAW index")
    symbol = header.get("element", "").strip()
    if symbol not in ELEMENTS:
        raise ValueError("Unknown element")
    ncore = ELEMENTS.index(symbol) + 1 - integer(header.get("z_valence", "nan"))
    if ncore <= 0 or sum(occupations) != ncore:
        raise ValueError(
            "Core occupations must sum to Z minus the explicit valence charge"
        )
    valence = set()
    for wave in tree.findall("PP_PSWFC/*"):
        if real(wave.get("occupation", "0")) <= 0:
            continue
        label = re.fullmatch(
            r"([0-9]+)([SPDFGH])", wave.get("label", "").strip().upper()
        )
        if label is None:
            raise ValueError("Cannot identify the occupied atomic valence states")
        valence.add((int(label[1]), "SPDFGH".index(label[2])))
    if not valence:
        raise ValueError("Missing occupied atomic valence-state labels")
    w = weights(rab)
    selected, seen = [], set()
    for i, occupation in enumerate(occupations, 1):
        node = table.find(f"PP_GIPAW_CORE_ORBITAL.{i}")
        if node is None or integer(node.get("index", "-1")) != i:
            raise ValueError("Invalid GIPAW core indexing")
        l, n = integer(node.get("l", "-1")), integer(node.get("n", "-1"))
        if l < 0 or n <= l or (n, l) in seen:
            raise ValueError("Invalid or repeated GIPAW orbital quantum numbers")
        seen.add((n, l))
        if occupation not in (0, 2 * (2 * l + 1)):
            raise ValueError(
                "This model requires explicitly selected closed core shells"
            )
        if not occupation:
            continue
        if (n, l) in valence:
            raise ValueError(
                "An occupied valence state cannot be assigned to the frozen core"
            )
        u = values(node)
        if (
            len(u) != len(r)
            or abs(math.fsum(a * b * b for a, b in zip(w, u)) - 1) > 1e-4
        ):
            raise ValueError("Original core orbital size or normalization check failed")
        radial = [a / b for a, b in zip(u, r) if b > 0]
        if r[0] == 0:
            if u[0] != 0:
                raise ValueError("A regular core orbital must have u(0)=0")
            origin = 0.0
            if l == 0:
                for j in range(3):
                    coefficient = math.prod(
                        -r[k + 1] / (r[j + 1] - r[k + 1]) for k in range(3) if j != k
                    )
                    origin += coefficient * radial[j]
            radial.insert(0, origin)
        derivative = slopes(r, radial)
        if r[0] == 0 and l > 1:
            derivative[0] = 0.0
        selected.append((occupation, l, radial, derivative))
    return selected


def model_kinetic(r, density, orbitals):
    """tau_M = tau_W[n_M] + n_M (tau_A - tau_W[n_A]) / n_A.

    The reference Pauli term is evaluated as positive orbital squares.
    Scaling the reference waves at each node avoids underflow in core tails.
    """
    derivative, kinetic = slopes(r, density, monotone=True), []
    for i, nm in enumerate(density):
        if nm == 0:
            if derivative[i] != 0:
                raise ValueError("A zero NLCC density has a nonzero model derivative")
            kinetic.append(0.0)
            continue
        scale = max(abs(wave[i]) for _, _, wave, _ in orbitals)
        if scale == 0:
            raise ValueError("Nonzero NLCC charge outside the reference core support")
        a = [wave[i] / scale for _, _, wave, _ in orbitals]
        b = [slope[i] / scale for _, _, _, slope in orbitals]
        norm = math.fsum(f * x * x for (f, _, _, _), x in zip(orbitals, a))
        eta = math.fsum(f * x * y for (f, _, _, _), x, y in zip(orbitals, a, b)) / norm
        residual = []
        for (f, l, _, _), x, y in zip(orbitals, a, b):
            angular = (
                l * (l + 1) * (x / r[i]) ** 2
                if r[i] > 0
                else (2 * y * y if l == 1 else 0.0)
            )
            residual.append(f * ((y - eta * x) ** 2 + angular))
        tau = (
            0.125 * (derivative[i] / math.sqrt(nm)) ** 2
            + nm * 0.5 * math.fsum(residual) / norm
        )
        if not math.isfinite(tau):
            raise ValueError("Nonfinite model core kinetic density")
        kinetic.append(tau)
    return kinetic


def enrich(data, occupations):
    tree = ET.fromstring(data)
    header = tree.find("PP_HEADER")
    if tree.tag != "UPF" or tree.get("version") != "2.0.1" or header is None:
        raise ValueError("Requires UPF version 2.0.1")
    if (
        not flag(header, "is_ultrasoft")
        or flag(header, "is_paw")
        or flag(header, "has_so")
    ):
        raise ValueError(
            "Requires scalar ultrasoft data without PAW or spin-orbit flags"
        )
    if tree.find("PP_TAUMOD") is not None:
        raise ValueError("Original already supplies model core kinetic data")
    nodes = [tree.find(tag) for tag in ("PP_MESH/PP_R", "PP_MESH/PP_RAB", "PP_NLCC")]
    if any(node is None for node in nodes):
        raise ValueError("Missing radial mesh, quadrature or NLCC density")
    r, rab, density = (values(node) for node in nodes)
    if len(r) < 4 or len(rab) != len(r) or len(density) != len(r):
        raise ValueError("Incomplete radial mesh or NLCC density")
    if (
        len(r) != int(header.get("mesh_size", "-1"))
        or r[0] < 0
        or any(b <= a for a, b in zip(r, r[1:]))
        or any(v <= 0 for v in rab)
    ):
        raise ValueError("Invalid original radial mesh or quadrature")
    if (
        not flag(header, "core_correction")
        or any(v < 0 for v in density)
        or not any(density)
    ):
        raise ValueError(
            "This orbital-rescaling model requires a nonzero nonnegative NLCC density"
        )
    kinetic = model_kinetic(r, density, reference_orbitals(tree, r, rab, occupations))
    spans, end = top_level_spans(data)
    first, last = spans["PP_HEADER"]
    old = data[first:last]
    pattern = rb"(?i)(with_metagga_info\s*=\s*)([\"'])(.*?)(\2)"
    if re.search(pattern, old):
        updated = re.sub(pattern, lambda m: m[1] + m[2] + b"T" + m[4], old, count=1)
    else:
        marker = old.rfind(b"/>")
        if marker < 0:
            raise ValueError("Requires a self-closing UPF header")
        updated = old[:marker] + b' with_metagga_info="T"' + old[marker:]
    provenance = (
        "\n<!-- CP2K auxiliary NLCC core kinetic model ORBITAL_RESCALE\n"
        f"Source SHA256 {sha256(data).hexdigest()}\n"
        f"Original GIPAW index occupations {' '.join(map(str, occupations))}\n"
        "Natural cubic reference-orbital slopes and shape-preserving NLCC slopes.\n"
        "Common orbital rescaling preserves the supplied NLCC charge.\n"
        "This field is a specified model, not recovered original generator tau. -->\n"
    )
    block = f'<PP_TAUMOD type="real" size="{len(r)}" columns="3" convention="hartree_per_volume">\n'
    block += "\n".join(
        " ".join(f"{v:.17e}" for v in kinetic[i : i + 3]) for i in range(0, len(r), 3)
    )
    block += "\n</PP_TAUMOD>\n"
    return (
        data[:first]
        + updated
        + data[last:end]
        + (provenance + block).encode()
        + data[end:]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--model", choices=["orbital-rescale"], required=True)
    parser.add_argument("--core-occupations", nargs="+", type=int, required=True)
    args = parser.parse_args()
    try:
        original = args.original.read_bytes()
        output = enrich(original, args.core_occupations)
        with args.output.open("xb") as stream:
            stream.write(output)
    except (OSError, ValueError, AttributeError, ET.ParseError) as error:
        parser.exit(1, f"Core kinetic preparation failed: {error}\n")
    print(f"Prepared auxiliary model core tau: {args.output}")
    print(f"Original SHA256: {sha256(original).hexdigest()}")
    print(f"Prepared SHA256: {sha256(output).hexdigest()}")


if __name__ == "__main__":
    main()
