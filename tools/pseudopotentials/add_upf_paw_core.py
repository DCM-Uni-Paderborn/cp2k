#!/usr/bin/env python3
"""Supplement scalar PAW core orbitals while retaining the original potential.

Use a generator export containing the actual closed core, not an orbital-list
prefix. Require the original radial mesh, AE core density and shared signed
core waves to agree within the explicit tolerance. This checks core data only,
not reproduction of the regenerated valence pseudopotential. Only the GIPAW
core-orbital table is replaced. OUTPUT must not already exist.
"""

import argparse
import copy
from hashlib import sha256
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from xml.parsers import expat

from add_upf_core_tau import integer, reference_orbitals, weights
from add_upf_full_waves import compare, flag, normalized, values


def core_span(data):
    parser = expat.ParserCreate()
    stack, spans = [], []
    start = None

    def opening(name, attributes):
        nonlocal start
        stack.append(name)
        if stack == ["UPF", "PP_GIPAW", "PP_GIPAW_CORE_ORBITALS"]:
            start = parser.CurrentByteIndex

    def closing(name):
        if stack == ["UPF", "PP_GIPAW", "PP_GIPAW_CORE_ORBITALS"]:
            end = parser.CurrentByteIndex
            if data[end : end + 2] == b"</":
                end = data.index(b">", end) + 1
            spans.append((start, end))
        stack.pop()

    parser.StartElementHandler = opening
    parser.EndElementHandler = closing
    parser.Parse(data, True)
    if len(spans) != 1:
        raise ValueError("Requires one GIPAW core-orbital table")
    return spans[0]


def core_waves(tree, mesh):
    table = tree.find("PP_GIPAW/PP_GIPAW_CORE_ORBITALS")
    if table is None or len(table) != int(table.get("number_of_core_orbitals", "-1")):
        raise ValueError("Invalid GIPAW core-orbital count")
    result = {}
    for i, node in enumerate(table, 1):
        key = integer(node.get("n", "-1")), integer(node.get("l", "-1"))
        if key in result or integer(node.get("index", "-1")) != i:
            raise ValueError("Repeated or invalid core-orbital index")
        wave = values(node)
        if len(wave) != mesh:
            raise ValueError("Incomplete core orbital")
        result[key] = wave
    return table, result


def enrich(original_data, generated_data, tolerance=0.0):
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Core tolerance must be finite and nonnegative")
    original, generated = (
        ET.fromstring(data) for data in (original_data, generated_data)
    )
    headers = []
    for tree in (original, generated):
        header = tree.find("PP_HEADER")
        if tree.tag != "UPF" or tree.get("version") != "2.0.1" or header is None:
            raise ValueError("Requires UPF version 2.0.1")
        if not flag(header, "is_paw") or flag(header, "has_so"):
            raise ValueError("Requires scalar PAW data")
        headers.append(header)
    for name in ("element", "functional", "relativistic", "z_valence", "mesh_size"):
        if normalized(headers[0].get(name, "")) != normalized(headers[1].get(name, "")):
            raise ValueError(f"Different core reference attribute {name}")
    for name in ("PP_MESH/PP_R", "PP_MESH/PP_RAB", "PP_PAW/PP_AE_NLCC"):
        a, b = original.find(name), generated.find(name)
        if a is None or b is None:
            raise ValueError(f"Missing required core data {name}")
        compare(a, b, tolerance if name.endswith("PP_AE_NLCC") else 0.0)
    r, rab = (values(original.find("PP_MESH/" + name)) for name in ("PP_R", "PP_RAB"))
    if (
        len(r) < 5
        or len(r) != integer(headers[0].get("mesh_size", "-1"))
        or len(rab) != len(r)
        or r[0] < 0
        or any(w <= 0 for w in rab)
        or any(y <= x for x, y in zip(r, r[1:]))
    ):
        raise ValueError("Invalid core radial mesh")
    _, a = core_waves(original, len(r))
    generated_table, b = core_waves(generated, len(r))
    common = a.keys() & b.keys()
    if not common:
        raise ValueError("No shared original core orbitals")
    for key in common:
        scale = max(max(map(abs, a[key])), max(map(abs, b[key])))
        difference = min(
            max(abs(x - sign * y) for x, y in zip(a[key], b[key])) for sign in (-1, 1)
        )
        if difference > tolerance * scale:
            raise ValueError(f"Different original core orbital {key}")
    combined = copy.deepcopy(original)
    parent = combined.find("PP_GIPAW")
    parent.remove(parent.find("PP_GIPAW_CORE_ORBITALS"))
    parent.append(copy.deepcopy(generated_table))
    occupations = [
        2 * (2 * integer(node.get("l", "-1")) + 1) for node in generated_table
    ]
    orbitals = reference_orbitals(combined, r, rab, occupations)
    density = values(original.find("PP_PAW/PP_AE_NLCC"))
    if len(density) != len(r):
        raise ValueError("Incomplete AE core density")
    reconstructed = [
        math.fsum(f * wave[i] ** 2 for f, _, wave, _ in orbitals) for i in range(len(r))
    ]
    target = [4 * math.pi * value for value in density]
    scale = max(map(abs, target))
    w = weights(rab)
    norm = math.fsum(weight * radius**2 * y for weight, radius, y in zip(w, r, target))
    if scale == 0 or norm <= 0:
        raise ValueError("Nonpositive AE core charge")
    maximum = max(abs(x - y) for x, y in zip(reconstructed, target)) / scale
    integral = math.fsum(
        weight * radius**2 * abs(x - y)
        for weight, radius, x, y in zip(w, r, reconstructed, target)
    )
    integral /= norm
    if max(maximum, integral) > tolerance:
        raise ValueError(
            "Complete core orbitals do not reproduce the original AE core density"
        )
    first, last = core_span(original_data)
    start, end = core_span(generated_data)
    return original_data[:first] + generated_data[start:end] + original_data[last:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("generated", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--core-tolerance", type=float, default=0.0)
    args = parser.parse_args()
    try:
        original, generated = args.original.read_bytes(), args.generated.read_bytes()
        result = enrich(original, generated, args.core_tolerance)
        with args.output.open("xb") as output:
            output.write(result)
    except (OSError, ValueError, KeyError, ET.ParseError, expat.ExpatError) as error:
        parser.error(str(error))
    for name, data in (
        ("Original", original),
        ("Generator", generated),
        ("Output", result),
    ):
        print(f"{name} SHA-256: {sha256(data).hexdigest()}")
    print(
        f"Preserved original potential and supplemented PAW core orbitals: {args.output}"
    )


if __name__ == "__main__":
    main()
