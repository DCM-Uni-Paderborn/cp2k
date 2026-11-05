#!/usr/bin/env python3
"""Download checksum-pinned SSSP v2.0 and prepared scalar USPP data.

Without --install, retain and extract the complete original SSSP archives.
With --install, use original NC/PAW files and prepared USPPs by default.
--uspp-source original installs all original SSSP files instead. Additional
PSLibrary and GBRV USPPs are available with --additional-uspp. All destinations are
chosen explicitly. Downloads do not enable unsupported calculation methods.
--paw-core supplies compatible complete frozen-core waves for selected PAW files.
"""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile
import urllib.request
import xml.etree.ElementTree as ET


def download(source, root, temporary):
    path = root / source["file"]
    if not path.exists():
        path = temporary / source["file"]
        if not path.exists():
            with urllib.request.urlopen(source["url"], timeout=120) as response:
                with path.open("wb") as output:
                    shutil.copyfileobj(response, output)
    checksum = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    if checksum.hexdigest() != source["sha256"]:
        raise RuntimeError(f"Checksum mismatch: {source['file']}")
    return path


def stage_file(data, target, staging, pending):
    if target.exists():
        if target.read_bytes() != data:
            raise RuntimeError(f"Existing file differs: {target}")
    elif target in pending:
        if pending[target].read_bytes() != data:
            raise RuntimeError(f"Conflicting installation data: {target}")
    else:
        staged = staging / f"upf-install-{len(pending):04d}"
        staged.write_bytes(data)
        pending[target] = staged


def install_staged(pending):
    for target, staged in pending.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        staged.replace(target)


class PreparedCollection:
    def __init__(self, root, temporary, core=False):
        self.root, self.temporary = root, temporary
        self.core = core
        self.release = json.loads(
            Path(__file__)
            .with_name("paw_core_release.json" if core else "inverse_uspp_release.json")
            .read_text()
        )
        self.destination = (
            root
            / ("prepared-paw-core" if core else "prepared-uspp")
            / self.release["version"]
        )
        self.catalog_data = download(
            self.release["catalog"], root, temporary
        ).read_bytes()
        self.entries = json.loads(self.catalog_data)["entries"]
        self.originals = {}
        for row in self.entries:
            path = PurePosixPath(row["file"])
            prefix = (
                ("data", "paw-core", row["functional"])
                if core
                else ("data", row["functional"])
            )
            if (
                path.is_absolute()
                or ".." in path.parts
                or path.parts[: len(prefix)] != prefix
            ):
                raise RuntimeError("Invalid prepared UPF path")
            if core or row["sssp"]:
                key = (row["functional"], row["original_sha256"])
                if key in self.originals:
                    raise RuntimeError("Duplicate original UPF in release catalog")
                self.originals[key] = row

    def stage(self, rows, staging, pending):
        stage_file(
            self.catalog_data, self.destination / "catalog.json", staging, pending
        )
        for functional in sorted({row["functional"] for row in rows}):
            wanted = {
                row["file"]: row for row in rows if row["functional"] == functional
            }
            source = self.release["archives"][functional]
            archive = download(source, self.root, self.temporary)
            metadata = {"LICENSE", "README.md"}
            if self.core:
                metadata.add("inputs/paw-core/qe-6.3-core-export.patch")
            else:
                metadata.add("sources.json")
            for row in rows:
                if row["functional"] != functional:
                    continue
                if self.core:
                    recipe = PurePosixPath(row["generator_input"])
                    if (
                        recipe.parts[:3] != ("inputs", "paw-core", functional)
                        or ".." in recipe.parts
                    ):
                        raise RuntimeError("Invalid PAW core generator-input path")
                    metadata.add(recipe.as_posix())
                elif row["source_library"] == "GBRV":
                    recipe = row["generator_input"]
                    if (
                        not recipe.startswith("inputs/gbrv/")
                        or ".." in PurePosixPath(recipe).parts
                    ):
                        raise RuntimeError("Invalid generator-input path")
                    metadata.update(f"{recipe}/{kind}.adat" for kind in ("ae", "ps"))
                reference = row.get("model_core_tau", {}).get("reference_file")
                if reference:
                    path = PurePosixPath(reference)
                    if (
                        path.parts[:2] != ("inputs", "core-kinetic")
                        or ".." in path.parts
                    ):
                        raise RuntimeError("Invalid core-reference path")
                    metadata.add(reference)
            with tarfile.open(archive) as contents:
                for member in contents:
                    try:
                        name = (
                            PurePosixPath(member.name)
                            .relative_to(self.release.get("archive_prefix", "."))
                            .as_posix()
                        )
                    except ValueError:
                        continue
                    if name not in wanted and name not in metadata:
                        continue
                    if not member.isfile():
                        raise RuntimeError(
                            f"Expected a regular release file: {member.name}"
                        )
                    with contents.extractfile(member) as input_file:
                        data = input_file.read()
                    if name in wanted:
                        row = wanted.pop(name)
                        if sha256(data).hexdigest() != row["sha256"]:
                            raise RuntimeError(
                                f"Prepared UPF checksum mismatch: {member.name}"
                            )
                    else:
                        metadata.remove(name)
                    stage_file(data, self.destination / name, staging, pending)
            if wanted or metadata:
                raise RuntimeError("Incomplete prepared UPF archive")


def install_library(root, name, source, prepared=None, paw_core=None):
    functional, table = name.split("-")
    destination = root / functional / {"eff": "efficiency", "prec": "precision"}[table]
    with tempfile.TemporaryDirectory(prefix="sssp-download-", dir=root) as temporary:
        staging = Path(temporary)
        archive = download(source, root, staging)
        names, pending, selection, replacements = set(), {}, {}, []
        core_replacements = []
        with tarfile.open(archive) as contents:
            for member in contents:
                if member.isdir():
                    continue
                filename = Path(member.name).name
                if (
                    not member.isfile()
                    or filename in names
                    or not (filename.endswith(".upf") or filename == "cutoffs.json")
                ):
                    raise RuntimeError(f"Unexpected archive member: {member.name}")
                names.add(filename)
                with contents.extractfile(member) as input_file:
                    data = input_file.read()
                target = destination / filename
                if filename.endswith(".upf"):
                    # Some original UPFs contain non-XML prose in PP_INFO.
                    match = re.search(
                        rb"""<PP_HEADER\b(?:[^>"']|"[^"]*"|'[^']*')*>""", data
                    )
                    if match is None:
                        raise RuntimeError(f"Missing UPF header: {filename}")
                    opening = match.group()
                    header = ET.fromstring(
                        opening if opening.endswith(b"/>") else opening[:-1] + b"/>"
                    )
                    pseudo_type = header.get("pseudo_type", "").strip().upper()
                    if pseudo_type not in ("NC", "PAW", "US", "USPP"):
                        raise RuntimeError(f"Unknown pseudopotential type: {filename}")
                    original_hash = sha256(data).hexdigest()
                    entry = {"type": pseudo_type, "original_sha256": original_hash}
                    if prepared is not None and pseudo_type in ("US", "USPP"):
                        row = prepared.originals.get((functional, original_hash))
                        if row is None or name not in row["sssp_collections"]:
                            raise RuntimeError(f"No matching prepared USPP: {filename}")
                        replacements.append(row)
                        target = prepared.destination / row["file"]
                        entry.update(sha256=row["sha256"], route=row["route"])
                    elif (
                        paw_core is not None
                        and pseudo_type == "PAW"
                        and (row := paw_core.originals.get((functional, original_hash)))
                        is not None
                    ):
                        core_replacements.append(row)
                        target = paw_core.destination / row["file"]
                        entry.update(
                            sha256=row["sha256"],
                            core_reference_tolerance=row["core_reference_tolerance"],
                        )
                    else:
                        stage_file(data, target, staging, pending)
                        entry["sha256"] = original_hash
                    entry["file"] = Path(
                        os.path.relpath(target, destination)
                    ).as_posix()
                    selection[filename] = entry
                else:
                    stage_file(data, target, staging, pending)
        if replacements:
            prepared.stage(replacements, staging, pending)
        if core_replacements:
            paw_core.stage(core_replacements, staging, pending)
        index = {
            "sssp_library": name,
            "sssp_archive_sha256": source["sha256"],
            "prepared_release": prepared.release["version"] if prepared else None,
            "files": selection,
        }
        if paw_core is not None:
            index["paw_core_release"] = paw_core.release["version"]
        stage_file(
            (json.dumps(index, indent=2, sort_keys=True) + "\n").encode(),
            destination / "selection.json",
            staging,
            pending,
        )
        install_staged(pending)
    print(
        f"{name}: installed {len(selection)} UPFs, {len(replacements)} prepared USPPs to {destination}",
        flush=True,
    )
    if core_replacements:
        print(
            f"{name}: supplied complete core waves for {len(core_replacements)} PAW datasets",
            flush=True,
        )


def main():
    manifest = json.loads(Path(__file__).with_name("sssp_v2_archives.json").read_text())
    libraries = {
        item["archive"]
        .removeprefix("SSSP-lib-")
        .removesuffix("-v2.tar.gz"): dict(item, file=item["archive"])
        for item in manifest["archives"]
    }
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="Explicit external data destination"
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="Install original NC/PAW and prepared USPP data, without retaining new archives",
    )
    parser.add_argument(
        "--uspp-source",
        choices=["prepared", "original"],
        help="USPP source for --install (default: prepared)",
    )
    parser.add_argument(
        "--additional-uspp",
        action="store_true",
        help="Also install non-SSSP PSLibrary and GBRV USPPs for the selected functionals",
    )
    parser.add_argument(
        "--paw-core",
        action="store_true",
        help="Supplement matching PAW files with complete frozen-core reference waves",
    )
    parser.add_argument(
        "--library",
        action="append",
        choices=sorted(libraries),
        help="Select a library; repeat to fetch more than one (default: all four)",
    )
    args = parser.parse_args()
    if args.uspp_source and not args.install:
        parser.error("--uspp-source requires --install")
    if args.additional_uspp and (not args.install or args.uspp_source == "original"):
        parser.error("--additional-uspp requires installation with prepared USPPs")
    if args.paw_core and not args.install:
        parser.error("--paw-core requires --install")
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    selected = list(dict.fromkeys(args.library or libraries))
    with tempfile.TemporaryDirectory(prefix="upf-download-", dir=root) as temporary:
        staging = Path(temporary)
        prepared = (
            PreparedCollection(root, staging)
            if args.install and args.uspp_source != "original"
            else None
        )
        paw_core = (
            PreparedCollection(root, staging, core=True) if args.paw_core else None
        )
        for name in selected:
            source = libraries[name]
            if args.install:
                install_library(root, name, source, prepared, paw_core)
                continue
            archive = download(source, root, staging)
            if archive.parent != root:
                archive = archive.replace(root / source["file"])
            destination = root / name
            destination.mkdir(exist_ok=True)
            with tarfile.open(archive) as contents:
                contents.extractall(destination, filter="data")
            print(f"{name}: verified and extracted to {destination}", flush=True)
        if args.additional_uspp:
            functionals = {name.split("-")[0] for name in selected}
            rows = [
                r
                for r in prepared.entries
                if not r["sssp"] and r["functional"] in functionals
            ]
            pending = {}
            prepared.stage(rows, staging, pending)
            install_staged(pending)
            print(f"Installed {len(rows)} additional, non-SSSP USPPs", flush=True)


if __name__ == "__main__":
    main()
