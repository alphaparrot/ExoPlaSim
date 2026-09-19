#!/usr/bin/env python
"""Remove reproducible work files from completed reference variants.

Native MOST fields, diagnostics, plots, and provenance remain in place.
Restart files are retained only for cases that explicitly test restart behavior.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from reference_cases.native_output import native_files


CASES = Path(__file__).resolve().parent
EXECUTABLE = re.compile(r"(?:most_plasim_t\d+_l\d+_p\d+|burn\d+)\.x\Z")
DIAGNOSTICS = ("ice_output", "ocean_output")
NATIVE_ARCHIVE = re.compile(r"MOST(?:_REST)?\.\d{5}\Z")
RESTART_CASES = frozenset({"D4", "I1", "I5"})
SAVED_RESTART = re.compile(r"MOST_REST\.\d{5}(?:\.gz)?\Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_file(path: Path, expected_hash: str | None = None) -> int:
    """Write, verify, and atomically install a lossless compressed copy."""
    target = path.with_name(path.name + ".gz")
    if target.exists():
        raise FileExistsError(target)
    original_hash = sha256(path)
    if expected_hash is not None and original_hash != expected_hash:
        raise RuntimeError(f"native output differs from recorded hash: {path}")
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.name + ".",
                                     suffix=".tmp", delete=False) as temporary:
        temporary_path = Path(temporary.name)
    try:
        with path.open("rb") as source, temporary_path.open("wb") as sink:
            with gzip.GzipFile(fileobj=sink, mode="wb", filename="", mtime=0,
                               compresslevel=6) as packed:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    packed.write(block)
        check = hashlib.sha256()
        with gzip.open(temporary_path, "rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                check.update(block)
        if check.hexdigest() != original_hash:
            raise RuntimeError(f"compressed data failed verification: {path}")
        packed_bytes = temporary_path.stat().st_size
        os.replace(temporary_path, target)
        path.unlink()
        return packed_bytes
    finally:
        temporary_path.unlink(missing_ok=True)


def compact_variant(workdir: Path, record: dict, dry_run: bool = False) -> dict:
    """Compact a successful reference workdir; safe to rerun."""
    if record.get("status") != "solver output captured":
        raise ValueError(f"not a completed reference variant: {workdir}")
    if not native_files(workdir):
        raise ValueError(f"native model output missing: {workdir}")
    before = sum(path.stat().st_size for path in workdir.rglob("*") if path.is_file())
    actions = []

    def remove(path: Path, reason: str) -> None:
        size = path.stat().st_size
        actions.append({"file": path.name, "action": "remove", "reason": reason,
                        "original_bytes": size})
        if not dry_run:
            path.unlink()

    for path in sorted(workdir.iterdir()):
        if not path.is_file() or path.is_symlink():
            continue
        if EXECUTABLE.fullmatch(path.name):
            remove(path, "copy of build artifact")
        elif path.suffix == ".sra" and path.name in record.get("inputs", {}):
            expected = record["inputs"][path.name]
            candidates = (workdir.parents[1] / "inputs").rglob(path.name)
            if sha256(path) != expected or not any(
                    source.is_file() and sha256(source) == expected
                    for source in candidates):
                raise RuntimeError(f"input differs from recorded hash: {path}")
            remove(path, "case input retained in inputs/")

    restart = workdir / "plasim_restart"
    saved_restarts = sorted(workdir.glob("MOST_REST.[0-9][0-9][0-9][0-9][0-9]"))
    if restart.is_file() and saved_restarts and sha256(restart) == sha256(saved_restarts[-1]):
        remove(restart, "identical to final MOST_REST segment")

    if record["case"] not in RESTART_CASES:
        for path in sorted(workdir.iterdir()):
            if path.is_file() and not path.is_symlink() and SAVED_RESTART.fullmatch(path.name):
                remove(path, "restart validated during run; case does not test restart behavior")

    for name in DIAGNOSTICS:
        path = workdir / name
        if not path.is_file():
            continue
        size = path.stat().st_size
        if dry_run:
            actions.append({"file": name, "action": "gzip", "original_bytes": size})
        else:
            packed_bytes = gzip_file(path)
            actions.append({"file": name, "action": "gzip", "original_bytes": size,
                            "stored_bytes": packed_bytes})

    for path in sorted(workdir.iterdir()):
        if not path.is_file() or not NATIVE_ARCHIVE.fullmatch(path.name):
            continue
        size = path.stat().st_size
        if dry_run:
            actions.append({"file": path.name, "action": "gzip", "original_bytes": size})
        else:
            expected = record.get("outputs", {}).get(path.name, {}).get("sha256")
            packed_bytes = gzip_file(path, expected)
            actions.append({"file": path.name, "action": "gzip", "original_bytes": size,
                            "stored_bytes": packed_bytes})

    existing = workdir / "storage_manifest.json"
    if not actions and existing.is_file():
        return json.loads(existing.read_text()) | {"actions": []}

    result = {"version": 1, "workdir": workdir.name, "original_bytes": before,
              "actions": actions}
    if not dry_run:
        result["stored_bytes"] = sum(path.stat().st_size for path in workdir.rglob("*")
                                     if path.is_file())
        existing.write_text(json.dumps(result, indent=2) + "\n")
    return result


def compact_crash(workdir: Path, dry_run: bool = False) -> int:
    """Discard only copied executables from an unsuccessful workdir."""
    saved = 0
    for path in workdir.iterdir():
        if path.is_file() and not path.is_symlink() and EXECUTABLE.fullmatch(path.name):
            saved += path.stat().st_size
            if not dry_run:
                path.unlink()
    return saved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", nargs="?", help="Comma-separated case IDs")
    parser.add_argument("--all", action="store_true", help="All cases with references")
    parser.add_argument("--dry-run", action="store_true", help="Preview without changing files")
    args = parser.parse_args()
    if (args.selection is None) == (not args.all):
        parser.error("provide a comma-separated case selection or --all")
    ids = [path.name for path in CASES.iterdir() if re.fullmatch(r"[A-Z]\d+", path.name)]
    selected = ids if args.all else [part.strip().upper() for part in args.selection.split(",")]
    if not selected or len(selected) != len(set(selected)) or any(case_id not in ids for case_id in selected):
        parser.error("unknown or duplicate case ID")
    saved = 0
    crash_saved = 0
    count = 0
    for case_id in sorted(selected):
        root = CASES / case_id / "reference"
        if not root.is_dir() or root.is_symlink():
            continue
        for record_path in sorted(root.glob("*.json")):
            if record_path.name == "provenance.json":
                continue
            record = json.loads(record_path.read_text())
            workdir = root / record_path.stem
            if not workdir.is_dir() or workdir.is_symlink():
                continue
            result = compact_variant(workdir, record, args.dry_run)
            removed = (sum(action["original_bytes"] - action.get("stored_bytes", 0)
                           for action in result["actions"])
                       if not args.dry_run else 0)
            if not args.dry_run:
                saved += removed
            count += 1
            print(f"{case_id}/{workdir.name}: {len(result['actions'])} actions"
                  + (" planned" if args.dry_run else f", saved {removed / 1048576:.1f} MiB"))
        for crash in sorted(root.glob("ref_*_crashed")):
            if crash.is_dir() and not crash.is_symlink():
                amount = compact_crash(crash, args.dry_run)
                crash_saved += amount
                if amount:
                    print(f"{case_id}/{crash.name}: copied executables"
                          f" {'would save' if args.dry_run else 'saved'} {amount / 1048576:.1f} MiB")
    print(f"{count} completed variant(s)" +
          (f" previewed; crash files would save {crash_saved / 1048576:.1f} MiB"
           if args.dry_run else
           f" compacted; saved {(saved + crash_saved) / 1073741824:.2f} GiB"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
