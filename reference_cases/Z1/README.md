# Z1 — Century Earth baseline

## What this case is

integrate Earth for 100 years at T21, T42, and T63, with annual mean output. T63 uses a 30-minute step because its 45-minute Earth-input smoke run traps in the native solver. Check long-term drift, energy balance, water cycle, and resolution dependence. The T63 Earth boundaries are nearest-neighbour remaps of bundled T42 data and must not be treated as independent high-resolution observations.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| t21_100y | T21 / 10 / 1 | 100 model year(s) | `plasim_namelist.NSTPW=11520` |
| t42_100y | T42 / 10 / 1 | 100 model year(s) | `plasim_namelist.NSTPW=11520` |
| t63_100y | T63 / 10 / 1 | 100 model year(s) | `timestep=30.0`, `plasim_namelist.NSTPW=17280` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- surface_temperature (139, K)
- toa_net_shortwave_flux (178, W m-2)
- toa_net_longwave_flux (179, W m-2)

Required fields vary by variant:

- `t21_100y`: 139, 167, 178, 179
- `t42_100y`: 139, 167, 178, 179
- `t63_100y`: 139, 167, 178, 179

The exact file inputs and any required diagnostic text are specified in `case.json`.

These long runs save one annual mean record per model year. They show long-term drift but do not resolve seasonal phase. A smoke run uses one 320-step segment and is not a scientific baseline. For long comparisons, full field differences are split into one compressed NumPy file per model year under `field_differences/`.

## Reference run

Pending. Capture with `python run_reference_cases.py Z1 --reference`.
