# CPU reference cases

The [capture status](REFERENCE_RUN_STATUS.md) lists successful references,
failed variants, and pending Z long runs.

This directory holds one folder for each of the 51 feature IDs in
[`todo/GPU_REFERENCE_FEATURES.md`](../todo/GPU_REFERENCE_FEATURES.md). Each folder has a
`case.json` manifest and an `inputs/` tree. The input links point to local copies
under `_assets/`, so the cases do not depend on the original data, spectrum, or
aerosol directories after repository cleanup. The runner copies only the listed
`.sra` files into the model work directory after `Model.configure` has finished.

From the repository root:

```sh
python run_reference_cases.py --list
python run_reference_cases.py --dry-run --all
python run_reference_cases.py D1 --reference
python run_reference_cases.py D1 --compare
python run_reference_cases.py D1,D2 --delete-runs --dry-run
python run_reference_cases.py D1,D2 --delete-runs
python run_reference_cases.py --all --delete-runs
python reference_cases/compact_reference_data.py D1 --dry-run
python reference_cases/compact_reference_data.py --all
# Later, after review:
# python run_reference_cases.py I1,X3,X4 --reference
# python run_reference_cases.py --all --reference
```

The default runs use the durations in the manifests. `--reference` writes to each
case's `reference/` directory and refuses to overwrite an existing baseline.
`--compare` (also the default mode) writes to a new `run_<UTC timestamp>/` directory
and compares each native output field, value by value, with the reference. The
comparison uses `rtol=1e-10`, `atol=1e-12`, and records whether values are *exactly*
equal. The full difference arrays are saved as compressed NumPy data, alongside a
field-by-field JSON report and side-by-side plots. Every run retains native output,
diagnostics, logs, runtimes, source commit, and a patch of tracked local edits.
Reference runs retain saved restart files for D4, I1, and I5; fresh test runs keep
their restarts. Both result directory patterns are ignored by Git. The source commit
alone may not recreate a baseline if `source.patch` is nonempty.

### Reference storage

The 43 completed reference variants, plus four failed-case folders, now use about
**0.34 GiB**, down from **3.57 GiB** before compaction. The native `MOST` output was
binary already, so converting it to NumPy `.npz` would require translating the
model's record format. We instead store it as lossless `MOST.00000.gz` (and likewise
compress `ice_output` and `ocean_output`). Saved restarts for D4, I1, and I5 are
compressed too. The comparison and plotting scripts read compressed or uncompressed
`MOST` files directly. Every retained output was verified against its original byte
count and SHA-256 hash before restart pruning. This preserves every output field and
the original native bytes.

Each completed work directory discards copies of compiled executables, the selected
`.sra` inputs already stored under the case's `inputs/`, and `plasim_restart` when
it exactly matches the saved final `MOST_REST` segment. Other cases discard
`MOST_REST.*` after the runner validates it; the field comparison never reads it.
Failed-case folders discard
only copied executables; their logs, inputs, and other crash evidence remain. The
per-variant `storage_manifest.json` records the latest compaction actions. Future
`--reference` runs compact each successful variant automatically. The standalone
compaction command above is safe to rerun and can preview changes with `--dry-run`.
Test runs remain uncompressed to make inspection easy; use `--delete-runs` to remove
them when no longer needed. To feed a retained restart back to the solver manually,
decompress `MOST_REST.00000.gz` first.

Every completed variant's `plots/summary.png` now shows surface and 2 m air
temperature in °C, liquid rainfall and precipitation in mm/day, and a diagnostic
specific to that case. Runs covering a 360-day model year also have
`plots/seasonal.png` with June and December monthly means. Raw model values and
field-by-field comparisons remain in native units. Existing plots can be
regenerated from stored `MOST` files without rerunning the solver; plot refresh
metadata records the plotting script and leaves the original solver runtimes intact.

`--delete-runs` removes only timestamped `run_*` directories under the selected
case folders and refreshes their READMEs. It keeps `reference/`, `case.json`,
and `inputs/`. Add `--dry-run` to preview the directories before deleting them;
`--all` selects every case.

`--smoke-steps 4` overrides most variants with a four-step plumbing check; NetCDF
conversion needs at least 160 steps to contain a raw time record. Smoke checks do
**not** establish scientific coverage and should not be used for a reference run.
For Z long runs, a smoke check uses just one segment and adjusts the write cadence
to its step count. Use at least 160 steps to exercise native output decoding.

The 95 variants in 51 cases include T21, T42 and a T63 aquaplanet, serial and two-process MPI,
several full-year runs, and optional physics switches. They can take substantial
time and disk space. The runner builds on first use of each resolution/layer/process
combination and recompiles when switching Earth and Mars builds or 4- and 8-byte
precision, because the upstream build gives these configurations the same
executable filename.

`X1` is explicitly blocked: the active solver allocates only one tracer slot and
uses it for humidity. `--all` reports this and exits nonzero even if the other runs
finish. Enabling a separate passive tracer requires a solver change and a new
initial-condition fixture. Other variants are **candidate** reference runs: some
specialist modes may need revised forcing or additional diagnostics before they
prove the feature is active. In particular, a short storm-capture run cannot count
as coverage unless a storm actually triggers, and a short glacier run cannot show
long-term ice evolution.

## Input provenance

- `_assets/T21` and `_assets/T42` contain the ten bundled Earth boundary files
  selected for their respective grids: codes 129, 169, 172, 173, 1730, 174, 210,
  212, 229 and 232. The T21 and T42 SST files are the original monthly fields.
- `_assets/T21_mars` contains the bundled Mars code-129 topography.
- `_assets/stellar` and `_assets/aerosol` contain copied stellar spectra and aerosol
  optical constants.
- `_assets/aqua_T21` and `_assets/aqua_T63` contain synthetic 280 K monthly SST
  grids (64×32 and 192×96). T63 has no bundled Earth boundary fields.
- `_assets/T63_earth` contains the ten T42 Earth boundary fields remapped to the
  192×96 T63 grid by nearest neighbour. They are generated by
  [`generate_t63_earth_inputs.py`](generate_t63_earth_inputs.py) and preserve the
  T42 source's spatial information, so T63 Earth results are a resolution and
  solver-path test rather than a high-resolution observed-data benchmark.
- Synthetic T21 codes 22, 130, 237, 709 and 903 provide fixed forcing for
  nudging, file ozone and ocean/ice flux-correction paths. These are deliberately
  simple plumbing fixtures, not observed climatologies. Each is a 64×32 text SRA
  with 14 records, except the ozone file, which has 10 levels × 14 records.

All `.sra` input bytes are hashed in the result JSON.
[`generate_synthetic_inputs.py`](generate_synthetic_inputs.py) recreates the constant
fields with eight-integer SRA headers and the values above. Replace these with
physically representative data before using these cases for scientific comparison.

## Long-run Z cases

Z1 requests 100 model years of default Earth at T21, T42, and T63. Z2 requests
10 years with dynamic vegetation, sea ice, and glaciers. Z3 requests 10 years
on a 45° obliquity, 0.1 eccentricity Keplerian orbit. The cases save one annual
mean field record per year (`NSTPW=11520` at 45 minutes for T21/T42 and
`NSTPW=17280` at 30 minutes for T63),
keeping complete annual fields while avoiding monthly or daily output over a
century. They are **not captured references yet**. Do not include them in a
routine `--all --reference` run unless you intend to run for many hours.

The Z1 T21 and T42 320-step smoke checks passed. T63 Earth crashed with SIGILL
at the default 45-minute step, but a 30-minute T63 run completed 320 steps and
decoded its output, so all Z T63 variants use 30 minutes. Z2 and Z3 T21 320-step
smoke checks passed, but
Z2's glacier configuration remains a longer-term candidate. The successful
short checks are outside `reference/` and are not reference data. A full-year
T21 check produced exactly one decoded annual field record, and a two-segment
short check produced two distinct `MOST` files. Long comparisons split their
full difference arrays into one compressed NumPy archive per model year under
`field_differences/`.

The short-run solver timings suggest that a 100-year Z1 reference would take
roughly 6 hours at T21, 27 hours at T42, and 4–5 days at T63 when run
sequentially on this machine. These are simple step-count extrapolations from
320-step smoke runs, not measured long-run runtimes. Plan disk space and run
them individually rather than launching the whole Z group at once.

## Four-byte precision stress cases

F1 repeats the Z1 T21 Earth configuration for 100 years with `precision=4`,
one annual mean per year, and automatic crash recovery disabled. F2 uses T42,
45° obliquity and 0.1 eccentricity for 30 years with `precision=4` and four
90-day mean records per year. The latter covers repeated seasonal extremes at
the higher grid resolution. Both cases record temperature, precipitation,
and top-of-atmosphere radiation and retain every native output field. A crash,
non-finite field, or implausible drift is evidence of instability; a completed
run alone does not establish 4-byte scientific accuracy. Compare F1 with the
matching 8-byte Z1 T21 trajectory and F2's first decade of annual means with
the 8-byte Z3 T42 case when those baselines are available.

These cases have **not been run**. Run them in a separate checkout after Z1
finishes, because changing precision recompiles an executable with the same
filename as the current 8-byte build. They are CPU precision tests, not GPU
tests. `--all` includes them, so select cases deliberately for long captures.
