#!/usr/bin/env python3
"""Add an explicitly chosen NLCC kinetic model to a scalar USPP or ONCV UPF.

GIPAW radial core orbitals, optionally completed by matching QE atomic waves,
or a retained Vanderbilt atom.ae file and closed-shell occupations define the
reference core. A common radial rescaling reproduces
the supplied nonnegative NLCC charge. This defines auxiliary model orbitals,
not recovered generator core tau or normalized atomic core eigenstates.
For ONCV data, a separate positive-amplitude model uses closed-core reference
waves from an ONCV export or a matching independent scalar atomic state.
It retains the original signed NLCC charge.
Original sections are retained, apart from the kinetic-data availability flag.
"""

import argparse
from bisect import bisect_right
from hashlib import sha256
import io
import math
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET

from add_upf_full_waves import (
    flag,
    integer,
    mask_upf_info,
    parse_upf,
    real,
    top_level_spans,
    values,
)

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


def resample(r, u, target):
    """Natural cubic interpolation, extending only a terminated zero tail by zero."""
    derivative = slopes(r, u)
    result = []
    for x in target:
        if x > r[-1]:
            if u[-1] != 0 or u[-2] != 0:
                raise ValueError("Target mesh exceeds a nonzero reference orbital tail")
            result.append(0.0)
            continue
        i = min(bisect_right(r, x) - 1, len(r) - 2)
        if i < 0:
            raise ValueError("Target mesh precedes the reference orbital mesh")
        h = r[i + 1] - r[i]
        t = (x - r[i]) / h
        result.append(
            (1 + 2 * t) * (1 - t) ** 2 * u[i]
            + t * (1 - t) ** 2 * h * derivative[i]
            + t * t * (3 - 2 * t) * u[i + 1]
            + t * t * (t - 1) * h * derivative[i + 1]
        )
    return result


def read_vanderbilt_state(data, z):
    """Read rwae's sequential records with four-byte markers and eight-byte reals."""
    if len(data) < 8:
        raise ValueError("Truncated Vanderbilt atomic state")
    endian = "<" if struct.unpack_from("<i", data)[0] == 56 else ">"
    offset = 0

    def record(fmt=None):
        nonlocal offset
        if offset + 8 > len(data):
            raise ValueError("Truncated Vanderbilt record")
        size = struct.unpack_from(endian + "i", data, offset)[0]
        if not 0 < size <= len(data) - offset - 8:
            raise ValueError("Invalid Vanderbilt record size")
        end = offset + 4 + size
        if struct.unpack_from(endian + "i", data, end)[0] != size:
            raise ValueError("Inconsistent Vanderbilt record markers")
        raw = data[offset + 4 : end]
        offset = end + 4
        if fmt is None:
            return raw
        if size != struct.calcsize(endian + fmt):
            raise ValueError("Unsupported Vanderbilt record layout")
        return struct.unpack(endian + fmt, raw)

    charge, ion, xc, mesh, relativity, aasf, bbsf, rmax = record("dddiiddd")
    if (
        charge != z
        or relativity not in (0, 2)
        or mesh < 5
        or mesh % 2 != 1
        or bbsf <= 0
        or rmax <= 0
        or not all(math.isfinite(x) for x in (charge, ion, xc, aasf, bbsf, rmax))
    ):
        raise ValueError("Invalid or incompatible scalar Vanderbilt atomic state")
    (count,) = record("i")
    if count <= 0 or count > (len(data) - offset) // 28:
        raise ValueError("Invalid Vanderbilt orbital count")
    shells = [record("idd") for _ in range(count)]
    if len(record()) != 8 * mesh:
        raise ValueError("Invalid Vanderbilt radial potential size")
    a, b = math.exp(-aasf) / z, 1 / bbsf
    r = [a * (math.exp(b * i) - 1) for i in range(mesh)]
    if (
        not all(math.isfinite(x) for x in r)
        or any(y <= x for x, y in zip(r, r[1:]))
        or r[-1] < rmax
    ):
        raise ValueError("Invalid Vanderbilt radial mesh")
    quadrature = weights([b * (x + a) for x in r])
    orbitals = []
    for label, occupied, energy in shells:
        raw = record()
        if len(raw) < 4:
            raise ValueError("Truncated Vanderbilt orbital")
        (size,) = struct.unpack_from(endian + "i", raw)
        if not 0 < size <= mesh or len(raw) != 4 + 8 * size:
            raise ValueError("Invalid Vanderbilt orbital size")
        u = list(struct.unpack_from(endian + str(size) + "d", raw, 4))
        u.extend([0.0] * (mesh - size))
        if not all(math.isfinite(x) for x in (*u, occupied, energy)) or occupied < 0:
            raise ValueError("Invalid Vanderbilt orbital values")
        orbitals.append((label // 100, label % 100 // 10, occupied, u))
    if offset != len(data):
        raise ValueError("Unexpected trailing Vanderbilt records")
    return ion, xc, relativity, r, quadrature, orbitals


def select_vanderbilt_core(orbitals, occupations, quadrature):
    if len(orbitals) != len(occupations):
        raise ValueError("Supply one core occupation for every Vanderbilt orbital index")
    selected_orbitals = []
    for (n, l, occupied, u), selected in zip(orbitals, occupations):
        if selected and (
            selected != occupied
            or abs(math.fsum(w * v * v for w, v in zip(quadrature, u)) - 1) > 1e-4
        ):
            raise ValueError("Vanderbilt core occupation or normalization check failed")
        selected_orbitals.append((n, l, u if selected else []))
    return selected_orbitals


def vanderbilt_orbitals(data, z, occupations, target):
    _, _, _, r, quadrature, orbitals = read_vanderbilt_state(data, z)
    return [
        (n, l, resample(r, u, target) if u else [])
        for n, l, u in select_vanderbilt_core(orbitals, occupations, quadrature)
    ]


def qe_atomic_orbitals(data, target):
    """Read scalar ld1.wfc columns on the original grid, printed to 12 decimals."""
    lines = data.decode("utf8").splitlines()
    if not lines or lines[0].split()[:2] != ["#", "r"]:
        raise ValueError("Requires the QE atomic wavefunction column header")
    labels = lines[0].split()[2:]
    shells = []
    for label in labels:
        match = re.fullmatch(r"([0-9]+)([SPDFGH])", label.upper())
        if match is None:
            raise ValueError("Invalid QE atomic orbital label")
        n, l = int(match[1]), "SPDFGH".index(match[2])
        if n <= l or (n, l) in shells:
            raise ValueError("Invalid or repeated QE atomic orbital")
        shells.append((n, l))
    rows = [[real(x) for x in line.split()] for line in lines[1:] if line.strip()]
    if (
        not shells
        or len(rows) != len(target)
        or any(
            len(row) != len(shells) + 1 or not all(map(math.isfinite, row))
            for row in rows
        )
    ):
        raise ValueError("Incomplete or nonfinite QE atomic wavefunction table")
    # Retain the exact UPF nodes. Only the radial column's decimal rounding is allowed.
    if any(
        abs(row[0] - radius) > 5.01e-13 + 4 * math.ulp(radius)
        for row, radius in zip(rows, target)
    ):
        raise ValueError("QE atomic wavefunctions use a different radial mesh")
    return [(n, l, [row[i + 1] for row in rows]) for i, (n, l) in enumerate(shells)]


def gipaw_core_orbitals(tree, r, rab, qe_wfc=None):
    table = tree.find("PP_GIPAW/PP_GIPAW_CORE_ORBITALS")
    if table is None:
        raise ValueError("Missing original GIPAW core orbital table")
    count = int(table.get("number_of_core_orbitals", "-1"))
    if len(table) != count:
        raise ValueError("Incomplete original GIPAW orbital table")
    orbitals = []
    for i in range(1, count + 1):
        node = table.find(f"PP_GIPAW_CORE_ORBITAL.{i}")
        if node is None or integer(node.get("index", "-1")) != i:
            raise ValueError("Invalid GIPAW core indexing")
        orbitals.append(
            (
                integer(node.get("n", "-1")),
                integer(node.get("l", "-1")),
                values(node),
            )
        )
    if qe_wfc is not None:
        original = {(n, l): u for n, l, u in orbitals}
        common = 0
        w = weights(rab)
        for n, l, u in qe_atomic_orbitals(qe_wfc, r):
            if (n, l) not in original:
                orbitals.append((n, l, u))
                continue
            a = original[n, l]
            if len(a) != len(r):
                raise ValueError("Incomplete original GIPAW orbital")
            norm = math.fsum(weight * x * x for weight, x in zip(w, a))
            error = min(
                math.fsum(
                    weight * (x - sign * y) ** 2 for weight, x, y in zip(w, a, u)
                )
                for sign in (-1, 1)
            )
            if not math.isfinite(norm) or not math.isfinite(error):
                raise ValueError("Nonfinite original GIPAW orbital")
            if norm == 0 and error == 0:
                continue
            if norm <= 0 or error > 1e-12 * norm:
                raise ValueError("QE and original GIPAW orbitals disagree")
            common += 1
        if not common:
            raise ValueError("QE atomic waves have no shared GIPAW reference")
    return orbitals


def reference_orbitals(tree, r, rab, occupations, vanderbilt=None, qe_wfc=None):
    header = tree.find("PP_HEADER")
    symbol = header.get("element", "").strip()
    if symbol not in ELEMENTS:
        raise ValueError("Unknown element")
    ncore = ELEMENTS.index(symbol) + 1 - integer(header.get("z_valence", "nan"))
    if ncore <= 0 or sum(occupations) != ncore:
        raise ValueError(
            "Core occupations must sum to Z minus the explicit valence charge"
        )
    if vanderbilt is not None and qe_wfc is not None:
        raise ValueError("Select one external atomic reference format")
    if vanderbilt is not None:
        orbitals = vanderbilt_orbitals(
            vanderbilt, ELEMENTS.index(symbol) + 1, occupations, r
        )
    else:
        orbitals = gipaw_core_orbitals(tree, r, rab, qe_wfc)
        if len(occupations) != len(orbitals):
            raise ValueError(
                "Supply occupations in GIPAW order, then missing QE header order"
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
    for (n, l, u), occupation in zip(orbitals, occupations):
        if l < 0 or n <= l or (n, l) in seen:
            raise ValueError("Invalid or repeated reference orbital quantum numbers")
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


def pauli_kinetic(r, density, orbitals):
    """Evaluate n_model (tau_A - tau_W[n_A]) / n_A as scaled squares."""
    scale = max(abs(wave) for _, _, wave, _ in orbitals)
    if scale == 0:
        raise ValueError("Nonzero NLCC charge outside the reference core support")
    a = [wave / scale for _, _, wave, _ in orbitals]
    b = [slope / scale for _, _, _, slope in orbitals]
    norm = math.fsum(f * x * x for (f, _, _, _), x in zip(orbitals, a))
    eta = math.fsum(f * x * y for (f, _, _, _), x, y in zip(orbitals, a, b)) / norm
    residual = []
    for (f, l, _, _), x, y in zip(orbitals, a, b):
        angular = (
            l * (l + 1) * (x / r) ** 2 if r > 0 else (2 * y * y if l == 1 else 0.0)
        )
        residual.append(f * ((y - eta * x) ** 2 + angular))
    return density * 0.5 * math.fsum(residual) / norm


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
        samples = [(f, l, wave[i], slope[i]) for f, l, wave, slope in orbitals]
        tau = 0.125 * (derivative[i] / math.sqrt(nm)) ** 2 + pauli_kinetic(
            r[i], nm, samples
        )
        if not math.isfinite(tau):
            raise ValueError("Nonfinite model core kinetic density")
        kinetic.append(tau)
    return kinetic


def oncv_document(data):
    """Return the UPF byte interval and mask unescaped generator notes for parsing."""
    match = re.search(rb"<UPF\b[^>]*>.*?</UPF>", data, re.S)
    if match is None:
        raise ValueError("Missing a complete ONCV UPF block")
    start, end = match.span()
    prefix = tuple(data[:start].strip().splitlines())
    suffix = data[end:].strip()
    if (prefix, suffix) not in (
        ((), b""),
        ((b"Begin PSP_UPF",), b"END_PSP"),
        ((b"END_PSP", b"Begin PSP_UPF"), b"END_PSP"),
    ):
        raise ValueError("Unrecognized framing around the ONCV UPF block")
    document = match[0]
    return start, document, mask_upf_info(document)


def oncv_configuration(document, header):
    """Read the occupied reference configuration and frozen-shell assignment."""
    match = re.search(rb"<PP_INPUTFILE\b[^>]*>(.*?)</PP_INPUTFILE>", document, re.S)
    if match is None:
        raise ValueError("Missing embedded ONCVPSP input")
    records = [
        line.split()
        for line in match[1].decode().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    element, z, nc, nv, xc, *_ = next(iter(records), ())
    z, nc, nv, xc = map(integer, (z, nc, nv, xc))
    q = real(header.get("z_valence", "nan"))
    if (
        element not in ELEMENTS
        or element != header.get("element", "").strip()
        or z != ELEMENTS.index(element) + 1
        or nc <= 0
        or nv <= 0
        or len(records) < 1 + nc + nv
    ):
        raise ValueError("Invalid embedded ONCV atomic configuration")
    orbitals = [
        (integer(n), integer(l), real(f)) for n, l, f, *_ in records[1 : 1 + nc + nv]
    ]
    core = orbitals[:nc]
    if (
        len({(n, l) for n, l, _ in orbitals}) != nc + nv
        or z - sum(f for _, _, f in core) != q
    ):
        raise ValueError("ONCV orbital assignment and frozen-core charge disagree")
    if any(l < 0 or n <= l or f != 2 * (2 * l + 1) for n, l, f in core):
        raise ValueError("The ONCV reference core must contain only closed shells")
    return element, z, nc, xc, orbitals


def oncv_reference(document, header, data):
    """Match the exported shells to the embedded atomic input before interpolation."""
    element, z, nc, xc, orbitals = oncv_configuration(document, header)
    stream = io.StringIO(data.decode("ascii"))
    schema, symbol, charge, count, nodes, functional = stream.readline().split()
    if (schema, symbol, integer(charge), integer(count), integer(functional)) != (
        "ONCV_AE_CORE_V1",
        element,
        z,
        nc,
        xc,
    ):
        raise ValueError("ONCV core export and embedded atomic input disagree")
    nodes = integer(nodes)
    if nodes < 4:
        raise ValueError("Incomplete ONCV reference mesh")
    native_r, reference = [], []
    for n, l, f in orbitals[:nc]:
        marker, nn, ll, ff = stream.readline().split()
        if (marker, integer(nn), integer(ll), real(ff)) != ("SHELL", n, l, f):
            raise ValueError("ONCV core export has a different occupied shell order")
        radii, reduced = [], []
        for _ in range(nodes):
            r, u = map(real, stream.readline().split())
            radii.append(r)
            reduced.append(u)
        if not native_r:
            if radii[0] <= 0 or any(b <= a for a, b in zip(radii, radii[1:])):
                raise ValueError("Invalid ONCV reference mesh")
            native_r = radii
        elif radii != native_r:
            raise ValueError("ONCV reference orbitals must share the same radial mesh")
        norm = math.fsum(
            (b - a) * (u * u + v * v) / 2 for a, b, u, v in zip(radii, radii[1:], reduced, reduced[1:])
        )
        if abs(norm - 1) > 5e-4:
            raise ValueError("ONCV reference orbital is not normalized on its native mesh")
        radial = [u / r for r, u in zip(radii, reduced)]
        reference.append((f, l, radial, slopes(radii, radial)))
    if stream.read().strip():
        raise ValueError("Unexpected trailing records in the ONCV core export")
    return native_r, reference


def oncv_vanderbilt_reference(document, header, data):
    """Match an independent neutral PBE atom and select the ONCV frozen shells."""
    _, z, nc, xc, configuration = oncv_configuration(document, header)
    ion, functional, relativity, r, quadrature, orbitals = read_vanderbilt_state(data, z)
    # ONCVPSP uses 4 for PBE and Vanderbilt uses 5. The archived neutral
    # scalar PBE states are independently generated reference orbitals.
    if (xc, ion, functional, relativity) != (4, 0, 5, 2):
        raise ValueError(
            "Independent Vanderbilt references require a matching neutral scalar PBE atom"
        )
    occupied = sorted((n, l, f) for n, l, f, _ in orbitals if f > 0)
    if occupied != sorted((n, l, f) for n, l, f in configuration if f > 0):
        raise ValueError("Vanderbilt and ONCV occupied configurations disagree")
    frozen = set(configuration[:nc])
    occupations = [f if (n, l, f) in frozen else 0 for n, l, f, _ in orbitals]
    reference = []
    selected = select_vanderbilt_core(orbitals, occupations, quadrature)
    for (_, l, u), f in zip(selected, occupations):
        if not f:
            continue
        if abs(u[0]) >= 1e-14:
            raise ValueError("Vanderbilt reduced radial orbital is nonzero at the origin")
        radial = [value / x for x, value in zip(r[1:], u[1:])]
        origin = 0.0
        if l == 0:
            # The regular s radial function can have a finite value at r=0.
            origin = math.fsum(
                radial[j]
                * math.prod(-r[k + 1] / (r[j + 1] - r[k + 1]) for k in range(3) if k != j)
                for j in range(3)
            )
        radial.insert(0, origin)
        reference.append((f, l, radial, slopes(r, radial)))
    return r, reference


def oncv_qe_reference(document, header, donor, log, qe_wfc):
    """Match an independent scalar QE atom and extract normalized frozen shells."""
    element, z, nc, xc, configuration = oncv_configuration(document, header)
    functional = {4: "PBE", -116133: "PBESOL"}.get(xc)
    tree, _ = parse_upf(donor)
    source = tree.find("PP_HEADER")
    if (
        functional is None
        or source is None
        or source.get("element", "").strip() != element
        or source.get("functional", "").strip().upper() != functional
        or source.get("relativistic", "").strip().lower() != "scalar"
        or flag(source, "has_so")
    ):
        raise ValueError("Independent QE references require matching scalar PBE or PBEsol data")
    report = re.search(rb" All-electron run (.*?) End of All-electron run ", log, re.S)
    if report is None:
        raise ValueError("Missing complete QE all-electron run")
    ae = report[1].decode("utf8")
    charge = re.search(r"atomic number is\s+([\d.]+)", ae)
    dft = re.search(r"\bdft\s*=\s*(\w+)", ae)
    error = re.search(r"final scf error:\s*(\S+)", ae)
    if (
        charge is None
        or real(charge[1]) != z
        or dft is None
        or dft[1].upper() != functional
        or error is None
        or not 0 <= real(error[1]) <= 1e-10
        or "scalar relativistic calculation" not in ae
    ):
        raise ValueError("QE log must identify a matching converged scalar atom")
    records = re.findall(r"(?m)^\s*(\d+)\s+(\d+)\s+\S+\s+\d+\(\s*([\d.]+)\s*\)", ae)
    occupied = [(int(n), int(l), real(f)) for n, l, f in records if real(f) > 0]
    if sorted(occupied) != sorted(row for row in configuration if row[2] > 0):
        raise ValueError("QE and ONCV occupied configurations disagree")

    r, rab = (values(tree.find(tag)) for tag in ("PP_MESH/PP_R", "PP_MESH/PP_RAB"))
    if (
        len(r) < 5
        or len(r) != int(source.get("mesh_size", "-1"))
        or len(rab) != len(r)
        or r[0] <= 0
        or any(b <= a for a, b in zip(r, r[1:]))
        or any(x <= 0 for x in rab)
    ):
        raise ValueError("Invalid QE reference radial mesh or quadrature")
    orbitals = gipaw_core_orbitals(tree, r, rab, qe_wfc)
    waves = {(n, l): u for n, l, u in orbitals}
    if len(waves) != len(orbitals):
        raise ValueError("Repeated QE reference orbital")
    quadrature = weights(rab)
    reference = []
    for n, l, f in configuration[:nc]:
        u = waves.get((n, l), [])
        norm = math.fsum(w * v * v for w, v in zip(quadrature, u))
        if len(u) != len(r) or not math.isfinite(norm) or abs(norm - 1) > 1e-4:
            raise ValueError("Missing or unnormalized QE closed-core wave")
        radial = [v / x for x, v in zip(r, u)]
        reference.append((f, l, radial, slopes(r, radial)))
    return r, reference


def positive_amplitude_kinetic(r, density, native_r, reference):
    """Build the auxiliary NC kinetic field while retaining the signed NLCC table."""
    if min(density) < -1e-8 * max(density):
        raise ValueError("Negative NLCC samples exceed 1e-8 of the positive peak")
    amplitude = [math.sqrt(max(value, 0.0)) for value in density]
    derivative = slopes(r, amplitude, monotone=True)
    kinetic = [0.0] * len(r)
    for i, radius in enumerate(r):
        if amplitude[i] == 0:
            if derivative[i] != 0:
                raise ValueError("A zero model amplitude has a nonzero derivative")
            continue
        if radius == 0:
            continue
        if not native_r[0] <= radius <= native_r[-1]:
            raise ValueError("Nonzero NLCC amplitude outside the ONCV reference mesh")
        j = min(bisect_right(native_r, radius) - 1, len(native_r) - 2)
        h = native_r[j + 1] - native_r[j]
        t = (radius - native_r[j]) / h
        samples = []
        for f, l, y, dy in reference:
            c0, c1 = y[j], h * dy[j]
            c2 = 3 * (y[j + 1] - c0) - h * (2 * dy[j] + dy[j + 1])
            c3 = 2 * (c0 - y[j + 1]) + h * (dy[j] + dy[j + 1])
            value = c0 + t * (c1 + t * (c2 + t * c3))
            slope = (c1 + t * (2 * c2 + t * 3 * c3)) / h
            samples.append((f, l, value, slope))
        kinetic[i] = 0.5 * derivative[i] ** 2 + pauli_kinetic(radius, amplitude[i] ** 2, samples)
    if r[0] == 0 and amplitude[0] > 0:
        # The origin carries zero radial volume weight.
        kinetic[0] = kinetic[1]
    if not all(map(math.isfinite, kinetic)):
        raise ValueError("Nonfinite ONCV model core kinetic density")
    return kinetic


def enrich(
    data,
    occupations=None,
    vanderbilt=None,
    qe_wfc=None,
    oncv_core=None,
    *,
    oncv_vanderbilt=None,
    oncv_qe=None,
):
    nc_inputs = (oncv_core, oncv_vanderbilt, oncv_qe)
    nc_model = any(value is not None for value in nc_inputs)
    if nc_model:
        start, document, masked = oncv_document(data)
        tree = ET.fromstring(masked)
    else:
        start, document = 0, data
        tree, masked = parse_upf(data)
    header = tree.find("PP_HEADER")
    if tree.tag != "UPF" or tree.get("version") != "2.0.1" or header is None:
        raise ValueError("Requires UPF version 2.0.1")
    if flag(header, "is_paw") or flag(header, "has_so"):
        raise ValueError("Requires scalar UPF data without PAW or spin-orbit flags")
    if not nc_model:
        if not flag(header, "is_ultrasoft") or occupations is None:
            raise ValueError(
                "Orbital rescaling requires ultrasoft data and explicit core occupations"
            )
    elif (
        header.get("pseudo_type", "").strip().upper() != "NC"
        or flag(header, "is_ultrasoft")
        or sum(value is not None for value in nc_inputs) != 1
        or any(value is not None for value in (occupations, vanderbilt, qe_wfc))
    ):
        raise ValueError(
            "The ONCV model requires scalar NC data and its own embedded shell occupations"
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
    if not flag(header, "core_correction") or max(density) <= 0:
        raise ValueError("Requires an enabled, nonzero NLCC density")
    if not nc_model:
        if min(density) < 0:
            raise ValueError("The USPP orbital-rescaling model requires a nonnegative NLCC density")
        kinetic = model_kinetic(
            r, density, reference_orbitals(tree, r, rab, occupations, vanderbilt, qe_wfc)
        )
    else:
        if oncv_core is not None:
            reference = oncv_reference(document, header, oncv_core)
        elif oncv_vanderbilt is not None:
            reference = oncv_vanderbilt_reference(document, header, oncv_vanderbilt)
        else:
            reference = oncv_qe_reference(document, header, *oncv_qe)
        kinetic = positive_amplitude_kinetic(r, density, *reference)
    if not nc_model:
        source = (
            "Original GIPAW"
            if vanderbilt is None
            else f"Vanderbilt atom.ae SHA256 {sha256(vanderbilt).hexdigest()}\nVanderbilt"
        )
        if qe_wfc is not None:
            source = f"QE ld1.wfc SHA256 {sha256(qe_wfc).hexdigest()}\nGIPAW plus QE"
        provenance = (
            "\n<!-- CP2K auxiliary NLCC core kinetic model ORBITAL_RESCALE\n"
            f"Source SHA256 {sha256(data).hexdigest()}\n"
            f"{source} index occupations {' '.join(map(str, occupations))}\n"
            "Natural cubic reference-orbital slopes and shape-preserving NLCC slopes.\n"
            "Common orbital rescaling preserves the supplied NLCC charge.\n"
            "This field is a specified model, not recovered original generator tau. -->\n"
        )
    else:
        if oncv_qe is not None:
            donor, log, waves = oncv_qe
            source = (
                f"QE reference UPF SHA256 {sha256(donor).hexdigest()}\n"
                f"QE all-electron log SHA256 {sha256(log).hexdigest()}\n"
            )
            if waves is not None:
                source += f"QE ld1.wfc SHA256 {sha256(waves).hexdigest()}\n"
            source += "Independent scalar QE reference matched to the ONCV atomic configuration.\n"
        elif oncv_vanderbilt is not None:
            source = (
                f"Vanderbilt atom.ae SHA256 {sha256(oncv_vanderbilt).hexdigest()}\n"
                "Independent neutral scalar PBE reference matched to the ONCV atomic configuration.\n"
            )
        else:
            source = (
                f"ONCV core export SHA256 {sha256(oncv_core).hexdigest()}\n"
                "Closed-core references follow the embedded ONCV atomic input.\n"
            )
        provenance = (
            "\n<!-- CP2K auxiliary NLCC core kinetic model POSITIVE_AMPLITUDE\n"
            f"Source SHA256 {sha256(data).hexdigest()}\n"
            f"{source}"
            "The original signed NLCC charge is retained. -->\n"
        )
    return insert_core_tau(data, start, document, masked, kinetic, provenance)


def insert_core_tau(data, start, document, masked, kinetic, provenance):
    spans, end = top_level_spans(masked)
    first, last = spans["PP_HEADER"]
    old = document[first:last]
    pattern = rb"(?i)(with_metagga_info\s*=\s*)([\"'])(.*?)(\2)"
    if re.search(pattern, old):
        updated = re.sub(pattern, lambda m: m[1] + m[2] + b"T" + m[4], old, count=1)
    else:
        marker = old.rfind(b"/>")
        if marker < 0:
            raise ValueError("Requires a self-closing UPF header")
        updated = old[:marker] + b' with_metagga_info="T"' + old[marker:]
    block = f'<PP_TAUMOD type="real" size="{len(kinetic)}" columns="3" convention="hartree_per_volume">\n'
    block += "\n".join(
        " ".join(f"{v:.17e}" for v in kinetic[i : i + 3])
        for i in range(0, len(kinetic), 3)
    )
    block += "\n</PP_TAUMOD>\n"
    return (
        data[: start + first]
        + updated
        + document[last:end]
        + (provenance + block).encode()
        + data[start + end:]
    )


def supplement_nc_core(data, kinetic, catalog_sha256):
    """Attach a tabulated positive-amplitude model to its matching NC original."""
    start, document, masked = oncv_document(data)
    tree = ET.fromstring(masked)
    header = tree.find("PP_HEADER")
    if (
        header is None
        or header.get("pseudo_type", "").strip().upper() != "NC"
        or any(flag(header, key) for key in ("is_ultrasoft", "is_paw", "has_so"))
        or not flag(header, "core_correction")
        or tree.find("PP_TAUMOD") is not None
    ):
        raise ValueError("Requires a scalar NC UPF with an unsupplemented core density")
    density = values(tree.find("PP_NLCC"))
    if (
        not isinstance(kinetic, list)
        or len(kinetic) != len(density)
        or len(kinetic) != int(header.get("mesh_size", "0"))
        or any(
            type(t) not in (int, float)
            or not math.isfinite(t)
            or t < 0
            or (n <= 0 and t != 0)
            for n, t in zip(density, kinetic)
        )
        or not any(kinetic)
    ):
        raise ValueError("Invalid tabulated NC core kinetic field")
    provenance = (
        "\n<!-- CP2K auxiliary NLCC core kinetic model POSITIVE_AMPLITUDE\n"
        f"Source SHA256 {sha256(data).hexdigest()}\n"
        f"SSSP NC core catalog SHA256 {catalog_sha256}\n"
        "The original signed NLCC charge is retained. -->\n"
    )
    return insert_core_tau(data, start, document, masked, kinetic, provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--model", choices=["orbital-rescale", "positive-amplitude"], required=True
    )
    parser.add_argument(
        "--core-occupations", nargs="+", type=int, help="Explicit USPP shell occupations"
    )
    reference = parser.add_mutually_exclusive_group()
    reference.add_argument(
        "--vanderbilt-ae",
        type=Path,
        help="Use a retained Vanderbilt atom.ae file for USPP or the positive-amplitude PBE NC model",
    )
    parser.add_argument(
        "--qe-ae-wfc",
        type=Path,
        help="Complete GIPAW with missing scalar ld1.wfc columns on the original mesh; "
        "for orbital-rescale, append occupations in missing QE header order",
    )
    reference.add_argument(
        "--oncv-core",
        type=Path,
        help="ONCV_AE_CORE_V1 closed-core export for the positive-amplitude NC model",
    )
    reference.add_argument(
        "--qe-reference-upf", type=Path, help="Matching independent QE UPF for the NC model"
    )
    parser.add_argument(
        "--qe-reference-log", type=Path, help="All-electron log from the same QE atomic run"
    )
    args = parser.parse_args()
    has_qe_reference = args.qe_reference_upf is not None
    if has_qe_reference != (args.qe_reference_log is not None):
        parser.error("--qe-reference-upf and --qe-reference-log must be supplied together")
    if args.model == "positive-amplitude" and (
        (args.oncv_core is None and args.vanderbilt_ae is None and not has_qe_reference)
        or args.core_occupations is not None
        or (args.qe_ae_wfc is not None and not has_qe_reference)
    ):
        parser.error(
            "positive-amplitude uses ONCV, Vanderbilt or QE references and the embedded ONCV occupations"
        )
    if args.model == "orbital-rescale" and (
        args.oncv_core is not None
        or has_qe_reference
        or args.core_occupations is None
        or (args.vanderbilt_ae is not None and args.qe_ae_wfc is not None)
    ):
        parser.error(
            "orbital-rescale requires --core-occupations and uses USPP reference data"
        )
    try:
        original = args.original.read_bytes()
        vanderbilt = args.vanderbilt_ae.read_bytes() if args.vanderbilt_ae else None
        qe_wfc = args.qe_ae_wfc.read_bytes() if args.qe_ae_wfc else None
        oncv_core = args.oncv_core.read_bytes() if args.oncv_core else None
        if args.model == "positive-amplitude":
            qe_reference = None
            if has_qe_reference:
                qe_reference = (
                    args.qe_reference_upf.read_bytes(),
                    args.qe_reference_log.read_bytes(),
                    qe_wfc,
                )
            output = enrich(
                original, oncv_core=oncv_core, oncv_vanderbilt=vanderbilt, oncv_qe=qe_reference
            )
        else:
            output = enrich(original, args.core_occupations, vanderbilt, qe_wfc)
        with args.output.open("xb") as stream:
            stream.write(output)
    except (OSError, ValueError, OverflowError, AttributeError, ET.ParseError) as error:
        parser.exit(1, f"Core kinetic preparation failed: {error}\n")
    print(f"Prepared auxiliary model core tau: {args.output}")
    print(f"Original SHA256: {sha256(original).hexdigest()}")
    print(f"Prepared SHA256: {sha256(output).hexdigest()}")


if __name__ == "__main__":
    main()
