# S10 — Vegetation

## What this case is

distinguish fixed/initial forest cover (`NVEG=0`), diagnostic SimBA vegetation (`NVEG=1`) and fully coupled vegetation (`NVEG=2`). For coupled mode, check biomass/soil carbon, forest fraction, albedo/roughness and soil-water effects.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| fixed_forest | T21 / 10 / 1 | 1 model year(s) | `vegetation=0` |
| diagnostic_simba | T21 / 10 / 1 | 1 model year(s) | `vegetation=1` |
| coupled_simba | T21 / 10 / 1 | 1 model year(s) | `vegetation=2` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- forest_cover (298, 1)
- vegetation_plant_carbon (304, kg C m-2)
- vegetation_soil_carbon (305, kg C m-2)

Required fields vary by variant:

- `fixed_forest`: 298
- `diagnostic_simba`: 298, 299
- `coupled_simba`: 298, 304, 305

The exact file inputs and any required diagnostic text are specified in `case.json`.

## Reference run

Attempted under [`reference/`](reference/).

Source commit: `bb46316ca917e907b129b0d764b11a9119e0088e`. `source.patch` records tracked local edits beyond that commit.

| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| fixed_forest | captured | 0.16 | 211.47 | 2.51 | 214.17 |
| diagnostic_simba | **failed** ([log](reference/diagnostic_simba.log)) | — | — | — | — |
| coupled_simba | captured | 0.13 | 206.82 | 2.12 | 209.10 |

Solver time excludes compilation and plotting. Variants run one at a time.
