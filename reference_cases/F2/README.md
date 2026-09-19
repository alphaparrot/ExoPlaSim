# F2 — Four-byte seasonal orbit stability

## What this case is

run a 45° obliquity, 0.1 eccentricity Keplerian Earth at T42 for 30 years with 4-byte default reals and four 90-day means per year. Check seasonal extremes, long-term drift, radiation and water cycle; compare the first decade with Z3's 8-byte T42 case after reducing the quarterly records to annual means. Disable automatic crash recovery so numerical failures remain visible.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| t42_f32_30y | T42 / 10 / 1 | 30 model year(s) | `obliquity=45.0`, `eccentricity=0.1`, `keplerian=True`, `plasim_namelist.NSTPW=2880`, `precision=4`, `crashtolerant=False` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- surface_temperature (139, K)
- air_temperature_2m (167, K)
- toa_net_shortwave_flux (178, W m-2)
- toa_net_longwave_flux (179, W m-2)
- lwe_of_large_scale_precipitation (142, m s-1)
- convective_precipitation_rate (143, m s-1)
- lwe_of_snowfall_amount (144, m s-1)

Required fields vary by variant:

- `t42_f32_30y`: 139, 167, 142, 143, 144, 178, 179

The exact file inputs and any required diagnostic text are specified in `case.json`.

These precision stress cases disable automatic crash recovery. A failed solver year must remain visible; finite output alone does not establish numerical equivalence to an 8-byte run. F1 saves annual means and F2 saves four 90-day means per year to expose seasonal extremes. For long comparisons, full field differences are split into one compressed NumPy file per model year under `field_differences/`.

## Reference run

Pending. Capture with `python run_reference_cases.py F2 --reference`.
