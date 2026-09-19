#!/usr/bin/env python
"""Write the run contract and captured baseline details for each case."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CHECKLIST = HERE.parent / "todo" / "GPU_REFERENCE_FEATURES.md"
sys.path.insert(0, str(HERE.parent))
from reference_cases.plot_reference import FEATURE_CODES, label  # noqa: E402


def description(case_id: str) -> str:
    pattern = re.compile(rf"^- \[ \] \*\*{re.escape(case_id)} — [^:]+:\*\* (.+)$", re.M)
    match = pattern.search(CHECKLIST.read_text())
    if not match:
        raise ValueError(f"No checklist description for {case_id}")
    return re.sub(r"\s+\(`[^)]*\)$", "", match.group(1)).strip()


def write_readme(case_id: str) -> None:
    case_dir = HERE / case_id
    case = json.loads((case_dir / "case.json").read_text())
    lines = [f"# {case_id} — {case['title']}", "", "## What this case is", "",
             description(case_id), "", "## What it tests", ""]
    if case.get("blocked_reason"):
        lines += [f"**Blocked:** {case['blocked_reason']}", "",
                  "No reference run exists for this feature yet.", ""]
    else:
        lines += ["| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |",
                  "| --- | --- | --- | --- |"]
        for variant in case["variants"]:
            model = variant.get("model", {})
            build = f"{model.get('resolution', 'T21')} / {model.get('layers', 10)} / {model.get('ncpus', 1)}"
            duration = (f"{variant['steps']} steps × {variant.get('segments', 1)} segments"
                        if "steps" in variant else f"{variant.get('segments', 1)} model year(s)")
            settings_parts = [f"`{key}={value}`" for key, value in
                              variant.get("configure", {}).items()]
            for namelist, fields in variant.get("namelists", {}).items():
                settings_parts.extend(f"`{namelist}.{key}={value}`"
                                      for key, value in fields.items())
            if model.get("mars"):
                settings_parts.append("Mars build")
            if model.get("precision", 8) != 8:
                settings_parts.append(f"`precision={model['precision']}`")
            if variant.get("crashtolerant") is False:
                settings_parts.append("`crashtolerant=False`")
            settings = (", ".join(settings_parts).replace("|", "\\|")) or "Earth defaults"
            lines.append(f"| {variant['name']} | {build} | {duration} | {settings} |")
        lines += ["", "Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.", ""]

        lines += ["## Outputs checked against", "",
                  "The solver's native `MOST` output is compared for **every output code and every "
                  "stored value**, including the time and vertical coordinates. Field presence and shape "
                  "must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; "
                  "`comparison.json` records exact equality, failures, maximum and RMS error, "
                  "and the worst index for each field. Compressed NumPy files keep the full "
                  "difference arrays for debugging. Native output, solver diagnostics, "
                  "and logs are retained in each run; D4, I1, and I5 references also retain "
                  "restart files. Their SHA-256 hashes are recorded for provenance, "
                  "but hashes do not determine the scientific comparison result.", "",
                  "The reference `plots/summary.png` shows surface and 2 m air temperature, "
                  "liquid rainfall, temperature and precipitation time series, and one "
                  "case-specific diagnostic. Temperature maps use Celsius and precipitation "
                  "uses mm/day; raw values retain the model's native units. Liquid rainfall "
                  "is calculated as codes 142 + 143 − 144 (large-scale and convective "
                  "precipitation minus snowfall). "
                  "Runs with monthly output covering a 360-day model year also get "
                  "`plots/seasonal.png` with June and December monthly maps. "
                  "Test runs also write `plots/comparison.png` with reference, test, "
                  "and difference panels for key fields. `reference_metrics.json` records field "
                  "shapes, ranges, moments, and per-record means or spectral RMS values.", "",
                  "The feature fields highlighted for this case are:", ""]
        for code in FEATURE_CODES[case_id]:
            lines.append(f"- {label(code)}")
        if any("check_codes" in variant for variant in case["variants"]):
            lines += ["", "Required fields vary by variant:", ""]
            for variant in case["variants"]:
                codes = variant.get("check_codes", FEATURE_CODES[case_id])
                lines.append(f"- `{variant['name']}`: " + ", ".join(codes))
        lines += ["", "The exact file inputs and any required diagnostic text are specified in "
                  "`case.json`.", ""]
        if case_id == "D1":
            lines += ["D1 covers 320 steps (10 model "
                      "days at 45 minutes per step), so it cannot show June or December. "
                      "The full-year D5 cases are intended for seasonal output.", ""]
        if case_id.startswith("Z"):
            lines += ["These long runs save one annual mean record per model year. "
                      "They show long-term drift but do not resolve seasonal phase. "
                      "A smoke run uses one 320-step segment and is not a scientific baseline. "
                      "For long comparisons, full field differences are split into one "
                      "compressed NumPy file per model year under `field_differences/`.", ""]
        if case_id.startswith("F"):
            lines += ["These precision stress cases disable automatic crash recovery. "
                      "A failed solver year must remain visible; finite output alone does not "
                      "establish numerical equivalence to an 8-byte run. F1 saves annual means "
                      "and F2 saves four 90-day means per year to expose seasonal extremes. "
                      "For long comparisons, full field differences are split into one "
                      "compressed NumPy file per model year under `field_differences/`.", ""]

        lines += ["## Reference run", ""]
        reference_dir = case_dir / "reference"
        if reference_dir.is_dir():
            complete = all((reference_dir / f"{variant['name']}.json").is_file()
                           for variant in case["variants"])
            lines += [("Captured" if complete else "Attempted") +
                      " under [`reference/`](reference/).", ""]
            provenance_path = reference_dir / "provenance.json"
            if provenance_path.is_file():
                provenance = json.loads(provenance_path.read_text())
                lines += [f"Source commit: `{provenance['commit']}`. `source.patch` records "
                          "tracked local edits beyond that commit.", ""]
            lines += ["| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |",
                      "| --- | --- | ---: | ---: | ---: | ---: |"]
            for variant in case["variants"]:
                name = variant["name"]
                record_path = reference_dir / f"{name}.json"
                log_path = reference_dir / f"{name}.log"
                if record_path.is_file():
                    record = json.loads(record_path.read_text())
                    t = record["timings_seconds"]
                    lines.append(f"| {name} | captured | {t['build']:.2f} | {t['solver']:.2f} | "
                                 f"{t['analysis']:.2f} | {t['total']:.2f} |")
                elif log_path.is_file():
                    crash = reference_dir / f"ref_{case_id}_{name}_crashed"
                    evidence = f"[log](reference/{name}.log)"
                    if crash.is_dir():
                        evidence += f", [crash files](reference/{crash.name}/)"
                    lines.append(f"| {name} | **failed** ({evidence}) | — | — | — | — |")
                else:
                    lines.append(f"| {name} | pending | — | — | — | — |")
            lines += ["", "Solver time excludes compilation and plotting. Variants run one at a time.", ""]
        else:
            lines += ["Pending. Capture with `python run_reference_cases.py " + case_id + " --reference`.", ""]
        runs = sorted(case_dir.glob("run_*"))
        if runs:
            lines += ["## Test runs", "",
                      "| Run | Variant | Solver (s) | Total (s) | Result |",
                      "| --- | --- | ---: | ---: | --- |"]
            for run in runs:
                for variant in case["variants"]:
                    record_path = run / f"{variant['name']}.json"
                    if not record_path.is_file():
                        continue
                    record = json.loads(record_path.read_text())
                    timings = record["timings_seconds"]
                    verdict = record.get("comparison", {})
                    result = ("PASS (exact)" if verdict.get("all_values_exact") else
                              "PASS (within tolerance)" if verdict.get("passed") else "FAIL")
                    lines.append(f"| [`{run.name}/`]({run.name}/) | {variant['name']} | "
                                 f"{timings['solver']:.2f} | {timings['total']:.2f} | {result} |")
            lines += [""]
    (case_dir / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    for manifest in HERE.glob("[A-Z][0-9]*/case.json"):
        write_readme(manifest.parent.name)
