# D1 — Spectral dynamics

## What this case is

vorticity, divergence, temperature, surface pressure and humidity evolve through the gridpoint/spectral transforms, semi-implicit step, and time filter. Compare 3-D fields and kinetic-energy or spectral diagnostics over multiple steps, not just end-of-year means.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| baseline | T21 / 10 / 1 | 320 steps × 1 segments | `plasim_namelist.seed=31415926` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- air_temperature (130, K)
- atm_relative_vorticity (138, s-1)
- divergence_of_wind (155, s-1)

The exact file inputs and any required diagnostic text are specified in `case.json`.

D1 covers 320 steps (10 model days at 45 minutes per step), so it cannot show June or December. The full-year D5 cases are intended for seasonal output.

## Reference run

Captured under [`reference/`](reference/).

Source commit: `bb46316ca917e907b129b0d764b11a9119e0088e`. `source.patch` records tracked local edits beyond that commit.

| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| baseline | captured | 23.85 | 6.10 | 0.51 | 30.47 |

Solver time excludes compilation and plotting. Variants run one at a time.

## Test runs

| Run | Variant | Solver (s) | Total (s) | Result |
| --- | --- | ---: | ---: | --- |
| [`run_20260917T122657841018Z/`](run_20260917T122657841018Z/) | baseline | 5.93 | 30.10 | PASS (exact) |
