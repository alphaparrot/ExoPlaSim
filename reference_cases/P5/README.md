# P5 — Stellar spectrum

## What this case is

solar/default forcing versus a different blackbody temperature and/or a supplied stellar spectrum. Compare the two shortwave-band weights and absorbed shortwave flux.

## What it tests

| Variant | Grid / layers / MPI ranks | Duration | Distinguishing settings |
| --- | --- | --- | --- |
| solar | T21 / 10 / 1 | 320 steps × 1 segments | Earth defaults |
| cool_blackbody | T21 / 10 / 1 | 320 steps × 1 segments | `startemp=3000.0` |
| wolf_spectrum | T21 / 10 / 1 | 320 steps × 1 segments | `starspec=@input/stellar/wolf.dat` |

Inputs are listed exactly in [`case.json`](case.json) and linked in `inputs/`.

## Outputs checked against

The solver's native `MOST` output is compared for **every output code and every stored value**, including the time and vertical coordinates. Field presence and shape must match. The comparison uses `rtol=1e-10` and `atol=1e-12`; `comparison.json` records exact equality, failures, maximum and RMS error, and the worst index for each field. Compressed NumPy files keep the full difference arrays for debugging. Native output, solver diagnostics, and logs are retained in each run; D4, I1, and I5 references also retain restart files. Their SHA-256 hashes are recorded for provenance, but hashes do not determine the scientific comparison result.

The reference `plots/summary.png` shows surface and 2 m air temperature, liquid rainfall, temperature and precipitation time series, and one case-specific diagnostic. Temperature maps use Celsius and precipitation uses mm/day; raw values retain the model's native units. Liquid rainfall is calculated as codes 142 + 143 − 144 (large-scale and convective precipitation minus snowfall). Runs with monthly output covering a 360-day model year also get `plots/seasonal.png` with June and December monthly maps. Test runs also write `plots/comparison.png` with reference, test, and difference panels for key fields. `reference_metrics.json` records field shapes, ranges, moments, and per-record means or spectral RMS values.

The feature fields highlighted for this case are:

- toa_net_shortwave_flux (178, W m-2)
- toa_net_longwave_flux (179, W m-2)

The exact file inputs and any required diagnostic text are specified in `case.json`.

## Reference run

Attempted under [`reference/`](reference/).

Source commit: `bb46316ca917e907b129b0d764b11a9119e0088e`. `source.patch` records tracked local edits beyond that commit.

| Variant | Status | Build (s) | Solver (s) | Analysis (s) | Total (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| solar | captured | 0.14 | 6.01 | 0.23 | 6.39 |
| cool_blackbody | captured | 0.19 | 5.89 | 0.22 | 6.32 |
| wolf_spectrum | **failed** ([log](reference/wolf_spectrum.log), [crash files](reference/ref_P5_wolf_spectrum_crashed/)) | — | — | — | — |

Solver time excludes compilation and plotting. Variants run one at a time.
