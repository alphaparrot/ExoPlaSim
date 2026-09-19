#!/usr/bin/env python
"""Run the CPU reference-case manifests in reference_cases/<ID>/case.json."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = ROOT / "reference_cases"
GROUP_ORDER = {group: index for index, group in enumerate("DPASXIZF")}
ID_RE = re.compile(r"^[A-Z][0-9]+$")
RUN_DIR_RE = re.compile(r"^run_\d{8}T\d{6}(?:\d{6})?Z$")
DEFAULT_SEED = "31415926"


@contextlib.contextmanager
def capture_output(log):
    """Capture Python messages and compiler/MPI child-process output in one log."""
    sys.stdout.flush()
    sys.stderr.flush()
    stdout_fd, stderr_fd = os.dup(1), os.dup(2)
    try:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            yield
    finally:
        log.flush()
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
        os.close(stdout_fd)
        os.close(stderr_fd)


def ordered_ids() -> list[str]:
    ids = [p.name for p in CASES.iterdir() if p.is_dir() and ID_RE.fullmatch(p.name)]
    return sorted(ids, key=lambda name: (GROUP_ORDER.get(name[0], 99), int(name[1:])))


def non_reference_runs(case_dir: Path) -> list[Path]:
    """Return only local timestamped run directories, never reference data."""
    return sorted(
        path
        for path in case_dir.iterdir()
        if RUN_DIR_RE.fullmatch(path.name) and path.is_dir() and not path.is_symlink()
    )


def read_case(case_id: str) -> dict:
    path = CASES / case_id / "case.json"
    case = json.loads(path.read_text())
    if case.get("id") != case_id:
        raise ValueError(f"{path}: case ID does not match its directory")
    for variant in case.get("variants", []):
        if not re.fullmatch(r"[a-z0-9_]+", variant.get("name", "")):
            raise ValueError(f"{path}: invalid variant name")
        if variant.get("model", {}).get("precision", 8) not in (4, 8):
            raise ValueError(f"{path}: precision must be 4 or 8 bytes")
        for relative in variant.get("inputs", []):
            source = (path.parent / "inputs" / relative).resolve()
            if not source.is_file() or not source.is_relative_to(CASES.resolve()):
                raise ValueError(f"{path}: missing or external input {relative}")
        names = [Path(item).name for item in variant.get("inputs", [])]
        if len(names) != len(set(names)):
            raise ValueError(f"{path}: duplicate input filenames in one variant")
        resolve_input_values(variant.get("configure", {}), path.parent / "inputs")
    if not case.get("variants") and not case.get("blocked_reason"):
        raise ValueError(f"{path}: no variants or explicit blocker")
    return case


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance(run_root: Path, case_id: str) -> dict:
    def git(*args: str) -> str:
        result = subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True
        )
        return result.stdout.strip() if result.returncode == 0 else "unavailable"

    patch = subprocess.run(
        [
            "git",
            "diff",
            "--binary",
            "HEAD",
            "--",
            ".",
            ":!exoplasim/plasim/run/most_plasim_run",
        ],
        cwd=ROOT,
        capture_output=True,
    )
    (run_root / "source.patch").write_bytes(
        patch.stdout if patch.returncode == 0 else b""
    )
    snapshot = run_root / "source_snapshot"
    snapshot.mkdir()
    scripts = [
        ROOT / "run_reference_cases.py",
        CASES / "plot_reference.py",
        CASES / "compare_reference.py",
        CASES / "native_output.py",
        CASES / "compact_reference_data.py",
        CASES / case_id / "case.json",
    ]
    captured = {}
    for source in scripts:
        relative = source.relative_to(ROOT)
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        captured[str(relative)] = sha256(target)
    compiler = subprocess.run(["gfortran", "--version"], capture_output=True, text=True)
    return {
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "tracked_diff_sha256": sha256(run_root / "source.patch"),
        "captured_script_sha256": captured,
        "untracked_files": git(
            "ls-files", "--others", "--exclude-standard"
        ).splitlines(),
        "python": sys.version,
        "platform": platform.platform(),
        "gfortran": compiler.stdout.splitlines()[0]
        if compiler.stdout
        else "unavailable",
        "runner_sha256": sha256(ROOT / "run_reference_cases.py"),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }


def resolve_input_values(value, input_dir: Path):
    if isinstance(value, str) and value.startswith("@input/"):
        path = (input_dir / value[len("@input/") :]).resolve()
        if not path.is_file() or not path.is_relative_to(CASES.resolve()):
            raise ValueError(f"Missing case input: {value}")
        return str(path)
    if isinstance(value, dict):
        return {
            key: resolve_input_values(item, input_dir) for key, item in value.items()
        }
    if isinstance(value, list):
        return [resolve_input_values(item, input_dir) for item in value]
    return value


def run_variant(
    case_id: str,
    variant: dict,
    case_dir: Path,
    case_output: Path,
    build_cache: dict,
    smoke_steps: int | None,
) -> dict:
    import exoplasim as exo
    from reference_cases.plot_reference import summarise

    started = time.perf_counter()
    name = variant["name"]
    workdir = case_output / name
    workdir.parent.mkdir(parents=True, exist_ok=True)
    input_dir = case_dir / "inputs"
    model_args = {"resolution": "T21", "layers": 10, "ncpus": 1, "precision": 8}
    model_args.update(variant.get("model", {}))
    build_key = (
        model_args["resolution"],
        model_args["layers"],
        model_args["ncpus"],
        model_args.get("mars", False),
        model_args["precision"],
    )
    executable_key = build_key[:3]
    # Build filenames omit both planet type and precision.
    recompile = build_cache.get(executable_key) != build_key[3:]
    configure = resolve_input_values(variant.get("configure", {}), input_dir)
    long_smoke = smoke_steps is not None and variant.get("long_run", False)
    segments = 1 if long_smoke else variant.get("segments", 1)
    steps = smoke_steps if smoke_steps is not None else variant.get("steps")
    if variant.get("postprocess") and steps is not None:
        # Fewer steps have no time record in the raw output, which pyburn
        # cannot decode even though the solver itself can finish them.
        steps = max(steps, 160)
    if steps is not None:
        configure["runsteps"] = steps
    configure.setdefault("snapshots", 0)
    namelists = {
        name: dict(fields) for name, fields in variant.get("namelists", {}).items()
    }
    namelists.setdefault("plasim_namelist", {}).setdefault("seed", DEFAULT_SEED)
    if long_smoke:
        # A partial year needs an output record for pyburn's field decoder.
        namelists["plasim_namelist"]["NSTPW"] = max(1, steps)
    records = {
        "case": case_id,
        "variant": name,
        "model": model_args,
        "configure": configure,
        "steps": steps,
        "namelists": namelists,
        "segments": segments,
        "crashtolerant": variant.get("crashtolerant", True),
        "inputs": {},
        "manifest_sha256": sha256(case_dir / "case.json"),
        "coverage_verified": False,
    }

    with (workdir.parent / f"{name}.log").open("w") as log, capture_output(log):
        build_started = time.perf_counter()
        model = exo.Model(
            workdir=str(workdir),
            modelname=f"ref_{case_id}_{name}",
            crashtolerant=variant.get("crashtolerant", True),
            recompile=recompile,
            **model_args,
        )
        build_seconds = time.perf_counter() - build_started
        build_cache[executable_key] = build_key[3:]
        configure_started = time.perf_counter()
        model.configure(**configure)
        if steps is not None:
            for key in ("N_RUN_YEARS", "N_RUN_MONTHS", "N_RUN_DAYS"):
                model._edit_namelist("plasim_namelist", key, "0")
        for namelist, fields in namelists.items():
            for field, value in fields.items():
                model._edit_namelist(namelist, field, str(value))

        # Model copies plasim/run wholesale. Replace those surface files with
        # this case's complete, recorded input selection after configuration.
        for old in workdir.glob("*.sra"):
            old.unlink()
        for relative in variant.get("inputs", []):
            source = (input_dir / relative).resolve()
            target = workdir / source.name
            if target.exists():
                raise ValueError(
                    f"Duplicate input filename in {case_id}/{name}: {target.name}"
                )
            shutil.copy2(source, target)
            records["inputs"][target.name] = sha256(source)

        configure_seconds = time.perf_counter() - configure_started
        solver_started = time.perf_counter()
        model.run(years=segments, postprocess=False, clean=False)
        solver_seconds = time.perf_counter() - solver_started
        conversion_started = time.perf_counter()
        if variant.get("postprocess", False):
            raw_path = workdir / "MOST.00000"
            exo.pyburn.postprocess(
                str(raw_path),
                str(raw_path) + ".nc",
                logfile=str(workdir / "pyburn.log"),
                variables=["ts", "ntr"],
                times=1,
                timeaverage=False,
            )
        conversion_seconds = time.perf_counter() - conversion_started

    analysis_started = time.perf_counter()
    diagnostic = workdir / "MOST_DIAG.00000"
    restart = workdir / "MOST_REST.00000"
    raw = workdir / "MOST.00000"
    expected_outputs = (diagnostic, restart, raw)
    # Multi-segment runs retain one raw/diagnostic/restart set per segment.
    for segment in range(segments):
        expected_outputs += tuple(
            workdir / f"{stem}.{segment:05d}"
            for stem in ("MOST_DIAG", "MOST_REST", "MOST")
        )
    for expected in set(expected_outputs):
        if not expected.is_file() or expected.stat().st_size == 0:
            raise RuntimeError(f"{case_id}/{name}: missing or empty {expected.name}")
    diagnostic_text = "\n".join(
        p.read_text(errors="replace") for p in workdir.glob("MOST_DIAG.*")
    )
    records["diagnostic_checks"] = {
        token: token in diagnostic_text for token in variant.get("require_diag", [])
    }
    if not all(records["diagnostic_checks"].values()):
        raise RuntimeError(f"{case_id}/{name}: expected diagnostic text not found")
    records["outputs"] = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(workdir.glob("MOST*"))
        if path.is_file()
    }
    if variant.get("postprocess") and not list(workdir.glob("MOST.00000.*")):
        raise RuntimeError(f"{case_id}/{name}: processed output missing")
    metrics = summarise(case_id, workdir, variant.get("check_codes"))
    records["feature_codes"] = metrics["feature_codes"]
    records["analysis_outputs"] = {
        str(path.relative_to(workdir)): {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in (
            workdir / "reference_metrics.json",
            *sorted((workdir / "plots").glob("*.png")),
        )
    }
    analysis_seconds = time.perf_counter() - analysis_started
    records["timings_seconds"] = {
        "build": build_seconds,
        "configure": configure_seconds,
        "solver": solver_seconds,
        "conversion": conversion_seconds,
        "analysis": analysis_seconds,
        "total": time.perf_counter() - started,
    }
    records["status"] = "solver output captured"
    (workdir.parent / f"{name}.json").write_text(json.dumps(records, indent=2) + "\n")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "selection", nargs="?", help="Comma-separated IDs, e.g. I1,X3,X4"
    )
    parser.add_argument("--all", action="store_true", help="Select all case IDs")
    parser.add_argument(
        "--list", action="store_true", help="List cases without running them"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate run inputs, or preview --delete-runs without deleting",
    )
    parser.add_argument(
        "--smoke-steps",
        type=int,
        help="Override every duration for a short plumbing check",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--reference",
        action="store_true",
        help="Capture each case's reference/ directory",
    )
    mode.add_argument(
        "--compare",
        action="store_true",
        help="Run and compare against reference/ (default)",
    )
    mode.add_argument(
        "--delete-runs",
        action="store_true",
        help="Delete timestamped run_* directories for selected cases; keep reference/",
    )
    args = parser.parse_args()
    ids = ordered_ids()
    if args.list and (args.selection or args.all or args.delete_runs):
        parser.error("--list cannot be combined with a selection or --delete-runs")
    if args.list:
        for case_id in ids:
            case = read_case(case_id)
            print(f"{case_id:3} {case['title']}")
        return 0
    if (args.selection is None) == (not args.all):
        parser.error("provide either comma-separated case IDs or --all")
    if args.smoke_steps is not None and args.smoke_steps < 1:
        parser.error("--smoke-steps must be positive")
    selected = (
        ids
        if args.all
        else [item.strip().upper() for item in args.selection.split(",")]
    )
    if (
        not selected
        or len(selected) != len(set(selected))
        or any(item not in ids for item in selected)
    ):
        parser.error("unknown or duplicate case ID; use --list to see valid IDs")
    if args.delete_runs:
        if args.smoke_steps is not None:
            parser.error("--smoke-steps cannot be combined with --delete-runs")
        from reference_cases.write_case_readmes import write_readme

        count = 0
        for case_id in selected:
            case_dir = CASES / case_id
            if case_dir.is_symlink() or case_dir.resolve().parent != CASES.resolve():
                parser.error(f"case directory is not local: {case_dir}")
            paths = non_reference_runs(case_dir)
            if not paths:
                print(f"{case_id}: no timestamped test runs")
                continue
            for path in paths:
                print(f"{'Would delete' if args.dry_run else 'Deleting'} {path}")
                if not args.dry_run:
                    shutil.rmtree(path)
                count += 1
            if not args.dry_run:
                write_readme(case_id)
        print(
            f"{'Would delete' if args.dry_run else 'Deleted'} {count} test run(s); reference/ preserved"
        )
        return 0
    cases = [read_case(case_id) for case_id in selected]
    if args.dry_run:
        for case in cases:
            state = (
                "BLOCKED"
                if case.get("blocked_reason")
                else f"{len(case['variants'])} variant(s)"
            )
            print(f"{case['id']}: {state}; inputs validated")
        return 0

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    mode_name = "reference" if args.reference else "compare"
    case_outputs = {}
    for case in cases:
        case_id = case["id"]
        if case.get("blocked_reason"):
            continue
        case_output = (
            CASES / case_id / ("reference" if args.reference else f"run_{stamp}")
        )
        if case_output.exists():
            parser.error(f"output directory already exists: {case_output}")
        if not args.reference and not (CASES / case_id / "reference").is_dir():
            parser.error(f"reference data missing for {case_id}; run --reference first")
        case_outputs[case_id] = case_output
    build_cache = {}
    failures = []
    for case in cases:
        case_id = case["id"]
        if case.get("blocked_reason"):
            print(f"{case_id}: BLOCKED — {case['blocked_reason']}")
            failures.append(case_id)
            continue
        case_failed = False
        case_output = case_outputs[case_id]
        case_output.mkdir(parents=True)
        (case_output / "provenance.json").write_text(
            json.dumps(
                {**git_provenance(case_output, case_id), "mode": mode_name}, indent=2
            )
            + "\n"
        )
        for variant in case["variants"]:
            label = f"{case_id}/{variant['name']}"
            print(f"Running {label} ...", flush=True)
            try:
                record = run_variant(
                    case_id,
                    variant,
                    CASES / case_id,
                    case_output,
                    build_cache,
                    args.smoke_steps,
                )
                if args.reference:
                    from reference_cases.compact_reference_data import compact_variant

                    compact_variant(case_output / variant["name"], record)
                else:
                    from reference_cases.compare_reference import compare

                    baseline = CASES / case_id / "reference" / variant["name"]
                    if not baseline.is_dir():
                        raise RuntimeError(f"missing reference variant: {baseline}")
                    comparison = compare(baseline, case_output / variant["name"])
                    record["comparison"] = {
                        "passed": comparison["passed"],
                        "all_values_exact": comparison["all_values_exact"],
                    }
                    (case_output / f"{variant['name']}.json").write_text(
                        json.dumps(record, indent=2) + "\n"
                    )
                    if not comparison["passed"]:
                        raise RuntimeError(
                            f"field comparison failed; see {case_output / variant['name'] / 'comparison.json'}"
                        )
                print(f"  finished: {case_output / (variant['name'] + '.json')}")
            except Exception as exc:
                print(exc)
                failures.append(label)
                case_failed = True
                print(f"  FAILED: {exc}", file=sys.stderr)
        if not case_failed:
            try:
                from reference_cases.write_case_readmes import write_readme

                write_readme(case_id)
            except Exception as exc:
                failures.append(case_id + "/summary")
                print(f"  FAILED case summary: {exc}", file=sys.stderr)
    for case_id, output in case_outputs.items():
        if output.exists():
            print(f"Outputs for {case_id}: {output}")
    if failures:
        print("Incomplete: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
