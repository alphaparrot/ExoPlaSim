#!/usr/bin/env python
"""Compare every native MOST field against a captured reference run."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from reference_cases.plot_reference import MM_PER_DAY, label, pyplot
from reference_cases.native_output import logical_name, native_files, read_native


RTOL = 1e-10
ATOL = 1e-12
PLOT_CODES = ("139", "167", "rainfall", "130", "138", "155")


def plot_comparison(reference: Path, current: Path, output: Path) -> None:
    from exoplasim.pyburn import refactorvariable

    ref_headers, ref_vars = read_native(native_files(reference)[-1])
    run_headers, run_vars = read_native(native_files(current)[-1])
    codes = [code for code in PLOT_CODES
             if (code in ref_vars and code in run_vars) or
             (code == "rainfall" and all(item in ref_vars and item in run_vars
                                         for item in ("142", "143", "144")))]
    if not codes:
        return
    plt = pyplot()
    fig, axes = plt.subplots(len(codes), 3, figsize=(15, 3.2 * len(codes)),
                             squeeze=False, constrained_layout=True)
    for row, code in enumerate(codes):
        def final_field(variables, headers):
            def field(item):
                array = refactorvariable(variables[item], headers[item],
                                         ntimes=len(variables["time"]),
                                         nlev=len(variables["sigmah"]))
                return np.asarray(array[-1]).reshape(-1, array.shape[-1])
            if code == "rainfall":
                return np.maximum((field("142") + field("143") - field("144"))
                                  * MM_PER_DAY, 0.0)
            result = field(code)
            return result - 273.15 if code in ("139", "167") else result

        a = final_field(ref_vars, ref_headers)
        b = final_field(run_vars, run_headers)
        if a.shape != b.shape:
            continue
        spectral = code in ("130", "138", "155")
        display_a = np.log10(np.abs(a) + 1e-12) if spectral else a
        display_b = np.log10(np.abs(b) + 1e-12) if spectral else b
        lo = float(min(display_a.min(), display_b.min()))
        hi = float(max(display_a.max(), display_b.max()))
        for col, field in enumerate((display_a, display_b)):
            im = axes[row, col].imshow(field, origin="lower", aspect="auto",
                                       cmap="viridis", vmin=lo, vmax=hi)
            fig.colorbar(im, ax=axes[row, col], shrink=0.75)
        difference = b - a
        bound = float(np.max(np.abs(difference)))
        im = axes[row, 2].imshow(difference, origin="lower", aspect="auto",
                                 cmap="coolwarm", vmin=-bound if bound else -1,
                                 vmax=bound if bound else 1)
        fig.colorbar(im, ax=axes[row, 2], shrink=0.75)
        display_label = ({"139": "Surface temperature (°C)",
                          "167": "2 m air temperature (°C)",
                          "rainfall": "Liquid rainfall (mm/day)"}.get(code, label(code)))
        for col, title in enumerate(("Reference", "Test run", "Test − reference")):
            axes[row, col].set_title(f"{display_label} — {title}")
            axes[row, col].set_xlabel("coefficient index" if spectral else "longitude index")
            axes[row, col].set_ylabel("model level" if spectral else "latitude index")
    fig.suptitle("Final output record; spectral fields show log10 magnitude in first two columns")
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output, dpi=140)
    plt.close(fig)


def compare(reference: Path, current: Path) -> dict:
    """Write complete native-field differences and a machine-readable verdict."""
    baseline_files = native_files(reference)
    test_files = native_files(current)
    report = {"rtol": RTOL, "atol": ATOL, "reference": str(reference),
              "test": str(current), "files": {}, "passed": True,
              "all_values_exact": True}
    if [logical_name(p) for p in baseline_files] != [logical_name(p) for p in test_files]:
        raise RuntimeError("Native MOST segment filenames differ from reference")
    split_differences = len(baseline_files) > 8
    differences: dict[str, np.ndarray] = {}
    if split_differences:
        (current / "field_differences").mkdir(exist_ok=True)
    for baseline_path, test_path in zip(baseline_files, test_files):
        segment_differences: dict[str, np.ndarray] = {}
        _, baseline = read_native(baseline_path)
        _, test = read_native(test_path)
        keys = set(baseline) | set(test)
        file_report = {"fields": {}, "missing_in_test": sorted(set(baseline) - set(test)),
                       "new_in_test": sorted(set(test) - set(baseline))}
        if file_report["missing_in_test"] or file_report["new_in_test"]:
            report["passed"] = False
        for code in sorted(keys & set(baseline) & set(test)):
            a = np.asarray(baseline[code], dtype=np.float64)
            b = np.asarray(test[code], dtype=np.float64)
            entry = {"reference_shape": list(a.shape), "test_shape": list(b.shape)}
            if a.shape != b.shape:
                entry["passed"] = False
                report["passed"] = False
            else:
                delta = b - a
                exact = bool(np.array_equal(a, b))
                close = np.isclose(a, b, rtol=RTOL, atol=ATOL, equal_nan=False)
                entry.update(exact=exact, passed=bool(close.all()),
                             values_outside_tolerance=int(np.count_nonzero(~close)),
                             max_absolute_error=float(np.max(np.abs(delta))) if delta.size else 0.0,
                             rms_error=float(np.sqrt(np.mean(delta * delta))) if delta.size else 0.0)
                if delta.size:
                    worst = int(np.argmax(np.abs(delta)))
                    entry["worst_index"] = [int(index) for index in np.unravel_index(worst, delta.shape)]
                    entry["reference_at_worst"] = float(a.flat[worst])
                    entry["test_at_worst"] = float(b.flat[worst])
                key = f"{logical_name(baseline_path).replace('.', '_')}_{code}"
                (segment_differences if split_differences else differences)[key] = delta
                report["passed"] &= entry["passed"]
                report["all_values_exact"] &= exact
            file_report["fields"][code] = entry
        report["files"][logical_name(baseline_path)] = file_report
        if split_differences:
            np.savez_compressed(
                current / "field_differences" / f"{logical_name(baseline_path)}.npz",
                **segment_differences,
            )
    if not split_differences:
        np.savez_compressed(current / "field_differences.npz", **differences)
    plot_comparison(reference, current, current / "plots" / "comparison.png")
    (current / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
