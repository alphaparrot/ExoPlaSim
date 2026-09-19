#!/usr/bin/env python
"""Summarise native MOST output and make compact case plots."""

from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path

import numpy as np

from reference_cases.native_output import logical_name, native_files, read_native


DEFAULT_CODES = ["139", "167", "142", "143", "144", "178", "179"]
MM_PER_DAY = 86_400_000.0


def pyplot():
    cache = Path(tempfile.gettempdir()) / "exoplasim-reference-matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt
FEATURE_CODES = {
    "D1": ["130", "138", "155"], "D2": ["139", "129"],
    "D3": ["130", "133"], "D4": ["139", "130"],
    "D5": ["139", "178"], "D6": ["130", "138"],
    "D7": ["133", "230"], "D8": ["139", "143"],
    "P1": ["139", "178", "179"], "P2": ["152", "130", "179"],
    "P3": ["139", "318"], "P4": ["178", "50", "52"],
    "P5": ["178", "179"], "P6": ["139", "152"],
    "A1": ["178", "179", "176", "177"],
    "A2": ["179", "178"], "A3": ["265", "179"],
    "A4": ["164", "178", "179"], "A5": ["175", "174", "184"],
    "A6": ["143", "164", "130"],
    "A7": ["142", "143", "144", "182"],
    "A8": ["146", "147", "180"], "A9": ["130", "138"],
    "S1": ["129", "172", "139"], "S2": ["172", "129", "139"],
    "S3": ["139", "170"], "S4": ["140", "160", "182"],
    "S5": ["141", "175", "139"], "S6": ["139", "176", "177"],
    "S7": ["139", "176"], "S8": ["210", "211", "139"],
    "S9": ["139", "176", "177"],
    "S10": ["298", "304", "305"],
    "S11": ["232", "267", "268"],
    "S12": ["266", "319", "152"],
    "X2": ["410", "411"], "X3": ["410", "178"],
    "X4": ["326", "327", "329"],
    "X5": ["326", "329"], "X6": ["130", "139"],
    "I1": ["139", "130"], "I2": ["139", "178", "179"],
    "I3": ["139", "178"], "I4": ["178", "179", "142", "143"],
    "I5": ["139", "178", "179"],
    "Z1": ["139", "178", "179"],
    "Z2": ["139", "298", "304", "305", "232", "267", "268"],
    "Z3": ["139", "178", "179"],
    "F1": ["139", "167", "178", "179", "142", "143", "144"],
    "F2": ["139", "167", "178", "179", "142", "143", "144"],
}


def label(code: str) -> str:
    from exoplasim.pyburn import ilibrary

    entry = ilibrary.get(code)
    return f"{entry[1]} ({code}, {entry[2]})" if entry else f"raw code {code}"


def climate_plots(case_id: str, workdir: Path, arrays: dict[str, np.ndarray],
                  time: list[float], fields: dict, required: list[str], plt) -> list[Path]:
    """Plot common climate fields and one diagnostic chosen for this case."""
    if "139" not in arrays or "167" not in arrays:
        raise RuntimeError("Surface or 2 m air temperature output is missing")
    has_precipitation = all(code in arrays for code in ("142", "143", "144"))
    surface_c = arrays["139"] - 273.15
    air_c = arrays["167"] - 273.15
    if has_precipitation:
        total = (arrays["142"] + arrays["143"]) * MM_PER_DAY
        snow = arrays["144"] * MM_PER_DAY
        rain = np.maximum(total - snow, 0.0)
    namelist = (workdir / "plasim_namelist").read_text()
    step_match = re.search(r"(?im)^\s*MPSTEP\s*=\s*([0-9.eE+-]+)", namelist)
    minutes_per_step = float(step_match.group(1)) if step_match else 45.0
    days = (np.asarray(time) + 1) * minutes_per_step / 1440.0
    plot_dir = workdir / "plots"
    plot_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    maps = ((axes[0, 0], surface_c[-1], "Surface temperature (°C)", "coolwarm"),
            (axes[0, 1], air_c[-1], "2 m air temperature (°C)", "coolwarm"))
    if has_precipitation:
        maps += ((axes[0, 2], rain[-1], "Liquid rainfall (mm/day)", "Blues"),)
    for ax, data, title, cmap in maps:
        im = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap)
        fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set(xlabel="longitude index", ylabel="latitude index", title=title)
    if not has_precipitation:
        axes[0, 2].text(0.5, 0.5, "Precipitation unavailable", ha="center", va="center")
        axes[0, 2].set_axis_off()
    axes[1, 0].plot(days, surface_c.mean(axis=(1, 2)), label="surface")
    axes[1, 0].plot(days, air_c.mean(axis=(1, 2)), label="2 m air")
    axes[1, 0].set(xlabel="elapsed days", ylabel="°C", title="Grid mean temperature")
    axes[1, 0].legend()
    if has_precipitation:
        axes[1, 1].plot(days, rain.mean(axis=(1, 2)), label="liquid rain")
        axes[1, 1].plot(days, snow.mean(axis=(1, 2)), label="snowfall")
        axes[1, 1].plot(days, total.mean(axis=(1, 2)), label="total precipitation")
        axes[1, 1].set(xlabel="elapsed days", ylabel="mm/day", title="Grid mean precipitation")
        axes[1, 1].legend()
    else:
        axes[1, 1].text(0.5, 0.5, "Precipitation unavailable", ha="center", va="center")
        axes[1, 1].set_axis_off()
    feature = next((code for code in required if code not in ("139", "167")), required[0])
    feature_series = np.asarray(fields[feature]["series_value"])
    feature_ylabel = fields[feature]["series_statistic"]
    feature_title = label(feature)
    if feature in ("139", "167"):
        feature_series = feature_series - 273.15
        feature_ylabel = "°C"
        feature_title = feature_title.replace("K)", "°C)")
    axes[1, 2].plot(days, feature_series)
    axes[1, 2].set(xlabel="elapsed days", ylabel=feature_ylabel,
                   title=f"Case diagnostic: {feature_title}")
    for ax in axes[1]:
        ax.grid(alpha=0.3)
    fig.suptitle(f"{case_id} — final output interval and time series ({days[-1]:.0f} elapsed days)")
    summary = plot_dir / "summary.png"
    fig.savefig(summary, dpi=140)
    plt.close(fig)

    plots = [summary]
    if has_precipitation and days[-1] >= 359.5:
        month = np.ceil(days / 30.0).astype(int) - 1  # 360-day model calendar
        if all(np.count_nonzero(month == index) >= 2 for index in (5, 11)):
            fields_by_column = ((surface_c, "Surface temperature (°C)", "coolwarm"),
                                (air_c, "2 m air temperature (°C)", "coolwarm"),
                                (rain, "Liquid rainfall (mm/day)", "Blues"))
            monthly = [[array[month == index].mean(axis=0) for index in (5, 11)]
                       for array, _, _ in fields_by_column]
            fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
            for col, (_, title, cmap) in enumerate(fields_by_column):
                low = 0.0 if col == 2 else min(float(field.min()) for field in monthly[col])
                high = max(float(field.max()) for field in monthly[col])
                for row, month_name in enumerate(("June", "December")):
                    ax = axes[row, col]
                    im = ax.imshow(monthly[col][row], origin="lower", aspect="auto",
                                   cmap=cmap, vmin=low, vmax=high)
                    fig.colorbar(im, ax=ax, shrink=0.8)
                    ax.set(xlabel="longitude index", ylabel="latitude index",
                           title=f"{month_name}: {title}")
            fig.suptitle(f"{case_id} — first model year, monthly mean output")
            seasonal = plot_dir / "seasonal.png"
            fig.savefig(seasonal, dpi=140)
            plt.close(fig)
            plots.append(seasonal)
    return plots


def summarise(case_id: str, workdir: Path, required_codes: list[str] | None = None) -> dict:
    from exoplasim.pyburn import refactorvariable
    plt = pyplot()

    files = native_files(workdir)
    if not files:
        raise RuntimeError("No native MOST output found for plots")
    codes = list(dict.fromkeys(["139", *FEATURE_CODES[case_id], *DEFAULT_CODES]))
    series: dict[str, list[np.ndarray]] = {code: [] for code in codes}
    time: list[float] = []
    grid_shape = None
    for path in files:
        headers, variables = read_native(path)
        grid_shape = (min(headers["main"][4:6]), max(headers["main"][4:6]))
        ntimes = len(variables["time"])
        if ntimes == 0:
            raise RuntimeError(f"{path.name} has no time records; increase its runsteps")
        segment_time = np.asarray(variables["time"], dtype=float)
        if time and segment_time[0] <= time[-1]:
            segment_time += time[-1]
        time.extend(segment_time.tolist())
        for code in codes:
            if code in variables:
                field = refactorvariable(variables[code], headers[code],
                                         ntimes=ntimes, nlev=len(variables["sigmah"]))
                series[code].append(np.asarray(field, dtype=float))

    fields = {}
    arrays = {}
    for code, segments in series.items():
        if not segments:
            continue
        array = np.concatenate(segments, axis=0)
        finite = np.isfinite(array)
        if not finite.any():
            raise RuntimeError(f"{code}: no finite output")
        if not finite.all():
            raise RuntimeError(f"{code}: non-finite output values")
        spectral = array.ndim == 3 and array.shape[-2:] != grid_shape
        if spectral:
            timeline = np.sqrt(np.mean(array * array, axis=tuple(range(1, array.ndim))))
            statistic = "spectral coefficient RMS"
        else:
            timeline = array.mean(axis=tuple(range(1, array.ndim)))
            statistic = "field mean"
        fields[code] = {
            "label": label(code), "shape": list(array.shape),
            "min": float(array.min()), "max": float(array.max()),
            "mean": float(array.mean()), "std": float(array.std()),
            "last_mean": float(array[-1].mean()),
            "series_statistic": statistic, "series_value": timeline.tolist(),
        }
        arrays[code] = array

    required = required_codes or FEATURE_CODES[case_id]
    missing = [code for code in required if code not in fields]
    if missing:
        raise RuntimeError(f"feature output codes missing: {', '.join(missing)}")
    result = {"time_steps": time, "fields": fields,
              "feature_codes": required, "raw_files": [logical_name(p) for p in files]}
    (workdir / "reference_metrics.json").write_text(json.dumps(result, indent=2) + "\n")

    climate_plots(case_id, workdir, arrays, time, fields, required, plt)
    return result


def plot_case(case_id: str, case_result_dir: Path, variants: list[str]) -> None:
    plt = pyplot()
    code = FEATURE_CODES[case_id][0]
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for name in variants:
        path = case_result_dir / name / "reference_metrics.json"
        metrics = json.loads(path.read_text())
        values = np.asarray(metrics["fields"][code]["series_value"])
        if code in ("139", "167"):
            values = values - 273.15
        ax.plot(metrics["time_steps"], values,
                marker="o", markersize=2, label=name)
    ax.set(xlabel="model step",
           ylabel="°C" if code in ("139", "167") else metrics["fields"][code]["series_statistic"],
           title=f"{case_id}: {label(code).replace('(139, K)', '(139, °C)').replace('(167, K)', '(167, °C)')}")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(case_result_dir / "comparison.png", dpi=140)
    plt.close(fig)
