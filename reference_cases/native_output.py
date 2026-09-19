"""Read native model output from either raw or losslessly compressed files."""

from __future__ import annotations

import gzip
import re
from pathlib import Path


NATIVE_NAME = re.compile(r"MOST\.\d{5}(?:\.gz)?\Z")


def logical_name(path: Path) -> str:
    return path.name[:-3] if path.name.endswith(".gz") else path.name


def native_files(workdir: Path) -> list[Path]:
    found = {}
    for path in workdir.iterdir():
        if not path.is_file() or not NATIVE_NAME.fullmatch(path.name):
            continue
        name = logical_name(path)
        if name in found:
            raise RuntimeError(f"both compressed and raw native output exist: {workdir / name}")
        found[name] = path
    return [found[name] for name in sorted(found)]


def read_native(path: Path):
    from exoplasim.pyburn import readallvariables

    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as source:
            return readallvariables(source.read())
    return readallvariables(path.read_bytes())
