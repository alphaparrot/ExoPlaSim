#!/usr/bin/env python
"""Create T63 Earth-like boundary files by nearest-neighbour remapping T42 files.

The source is the bundled T42 Earth forcing copied to ``_assets/T42``. This is a
resolution-matched fixture for exercising T63; it adds no new spatial information.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


ASSETS = Path(__file__).resolve().parent / "_assets"
SOURCE = ASSETS / "T42"
DESTINATION = ASSETS / "T63_earth"
SOURCE_SHAPE = (64, 128)
TARGET_SHAPE = (96, 192)


def remap(source: Path, target: Path) -> None:
    lines = source.read_text().splitlines()
    values_per_line = len(lines[1].split())
    rows = SOURCE_SHAPE[0] * SOURCE_SHAPE[1] // values_per_line
    if values_per_line not in (4, 8) or len(lines) % (rows + 1):
        raise ValueError(f"unexpected SRA record layout in {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    y = np.floor((np.arange(TARGET_SHAPE[0]) + 0.5) *
                 SOURCE_SHAPE[0] / TARGET_SHAPE[0]).astype(int)
    x = np.floor((np.arange(TARGET_SHAPE[1]) + 0.5) *
                 SOURCE_SHAPE[1] / TARGET_SHAPE[1]).astype(int)
    with target.open("w") as output:
        for start in range(0, len(lines), rows + 1):
            header = [int(value) for value in lines[start].split()]
            if len(header) != 8 or tuple(header[4:6]) != (128, 64):
                raise ValueError(f"unexpected SRA header in {source}: {header}")
            values = np.fromstring(" ".join(lines[start + 1:start + 1 + rows]), sep=" ")
            if values.size != SOURCE_SHAPE[0] * SOURCE_SHAPE[1]:
                raise ValueError(f"incomplete SRA record in {source} at line {start + 1}")
            header[4:6] = [192, 96]
            output.write("".join(f"{value:10d}" for value in header) + "\n")
            field = values.reshape(SOURCE_SHAPE)[np.ix_(y, x)].ravel()
            for offset in range(0, field.size, 8):
                output.write("".join(f"{value:14.6E}" for value in field[offset:offset + 8]) + "\n")


def main() -> None:
    for source in sorted(SOURCE.glob("N064_surf_*.sra")):
        target = DESTINATION / source.name.replace("N064", "N096", 1)
        remap(source, target)
        print(target)


if __name__ == "__main__":
    main()
