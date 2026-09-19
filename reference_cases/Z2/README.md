# Z2 — Coupled vegetation, sea ice and glaciers

## What this case is

integrate dynamic vegetation, sea ice, and evolving glaciers for 10 years at T21, T42, and T63. Check long-term surface temperature, forest and soil state, sea ice, and ice sheet evolution. The glacier path currently has a short-case crash (S11), so this remains a candidate until that failure is resolved.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| t21_10y | T21 / 10 / 1 | 10 model year(s) | `vegetation=2`, `seaice=True`, `glaciers={'toggle': True, 'mindepth': 2.0, 'initialh': -1.0}`, `plasim_namelist.NSTPW=11520` |
| t42_10y | T42 / 10 / 1 | 10 model year(s) | `vegetation=2`, `seaice=True`, `glaciers={'toggle': True, 'mindepth': 2.0, 'initialh': -1.0}`, `plasim_namelist.NSTPW=11520` |
| t63_10y | T63 / 10 / 1 | 10 model year(s) | `vegetation=2`, `seaice=True`, `glaciers={'toggle': True, 'mindepth': 2.0, 'initialh': -1.0}`, `timestep=30.0`, `plasim_namelist.NSTPW=17280` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- surface_temperature (139, K)
- forest_cover (298, 1)
- vegetation_plant_carbon (304, kg C m-2)
- vegetation_soil_carbon (305, kg C m-2)
- glacier_cover (232, 1)
- ground_geopotential (267, m2 s-2)
- glacier_geopotential (268, m2 s-2)

Required fields vary by variant:

- `t21_10y`: 139, 298, 304, 305, 232, 267, 268
- `t42_10y`: 139, 298, 304, 305, 232, 267, 268
- `t63_10y`: 139, 298, 304, 305, 232, 267, 268

The exact file inputs and any required diagnostic text are specified in `case.json`.

These long runs save one annual mean record per model year. They show long-term drift but do not resolve seasonal phase. A smoke run uses one 320-step segment and is not a scientific baseline. For long comparisons, full field differences are split into one compressed NumPy file per model year under `field_differences/`.

## Reference run

Pending. Capture with `python run_reference_cases.py Z2 --reference`.
