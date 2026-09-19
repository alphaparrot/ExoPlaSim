# S12 — Carbon–silicate weathering

## What this case is

cover off/on, runoff/temperature-dependent weathering and supply-limited weathering. For `evolveco2`, cross the annual update and check CO2 and pressure changes; a sub-year case cannot test it.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| weathering_off | T21 / 10 / 1 | 1 model year(s) | Earth defaults |
| weathering_on | T21 / 10 / 1 | 1 model year(s) | `co2weathering=True` |
| supply_limited_evolving_co2 | T21 / 10 / 1 | 2 model year(s) | `co2weathering=True`, `evolveco2=True`, `erosionsupplylimit=1.0` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- local_weathering (266, W_earth)
- weatherable_precipitation (319, mm day-1)
- log_surface_pressure (152, 1)

The exact file inputs and any required diagnostic text are specified in `case.json`.

## Reference run

Attempted under [`reference/`](reference/).

Source commit: `bb46316ca917e907b129b0d764b11a9119e0088e`. `source.patch` records tracked local edits beyond that commit.

| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| weathering_off | captured | 0.19 | 207.52 | 2.36 | 210.09 |
| weathering_on | captured | 0.19 | 217.56 | 3.72 | 221.50 |
| supply_limited_evolving_co2 | **failed** ([log](reference/supply_limited_evolving_co2.log), [crash files](reference/ref_S12_supply_limited_evolving_co2_crashed/)) | — | — | — | — |

Solver time excludes compilation and plotting. Variants run one at a time.
