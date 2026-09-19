# S11 — Glaciers and evolving topography

## What this case is

if enabled, cover initial ice sheet or persistent snow, glacier formation/ablation and its effect on orography and atmospheric fields. A short case may need deliberately chosen initial conditions to trigger this path.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| initial_ice_sheet | T21 / 10 / 1 | 1 model year(s) | `glaciers={'toggle': True, 'mindepth': 2.0, 'initialh': 10.0}` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- glacier_cover (232, 1)
- ground_geopotential (267, m2 s-2)
- glacier_geopotential (268, m2 s-2)

The exact file inputs and any required diagnostic text are specified in `case.json`.

## Reference run

Attempted under [`reference/`](reference/).

Source commit: `bb46316ca917e907b129b0d764b11a9119e0088e`. `source.patch` records tracked local edits beyond that commit.

| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| initial_ice_sheet | **failed** ([log](reference/initial_ice_sheet.log), [crash files](reference/ref_S11_initial_ice_sheet_crashed/)) | — | — | — | — |

Solver time excludes compilation and plotting. Variants run one at a time.
