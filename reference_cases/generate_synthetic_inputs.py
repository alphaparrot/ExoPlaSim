#!/usr/bin/env python
"""Recreate the small, constant SRA forcing fixtures used by the case manifests."""

from pathlib import Path


HERE = Path(__file__).resolve().parent / "_assets"
GRIDS = {"T21": (64, 32), "T63": (192, 96)}


def write_sra(subdir: str, grid: str, code: int, value: float, records: int = 14) -> None:
    nlon, nlat = GRIDS[grid]
    filename = f"N{nlat:03d}_surf_{code:04d}.sra"
    destination = HERE / subdir / filename
    destination.parent.mkdir(parents=True, exist_ok=True)
    header = "".join(f"{part:10d}" for part in (code, 0, 20090000, -1, nlon, nlat, 0, 0))
    row = "".join(f"{value:14.6E}" for _ in range(8))
    with destination.open("w") as output:
        for _ in range(records):
            output.write(header + "\n")
            for _ in range(nlon * nlat // 8):
                output.write(row + "\n")


if __name__ == "__main__":
    write_sra("aqua_T21", "T21", 169, 280.0)
    write_sra("aqua_T63", "T63", 169, 280.0)
    for code, value, count in ((237, 1e-6, 140), (130, 280.0, 14),
                               (22, 1e-6, 14), (903, 10.0, 14), (709, 10.0, 14)):
        write_sra("T21", "T21", code, value, count)
