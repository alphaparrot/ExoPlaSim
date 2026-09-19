# F1 f32 assessment

## Comparison

F1 completed 30 T21L10 model years in four-byte precision with automatic crash
recovery disabled. The executable used bounds checking and traps for invalid,
divide-by-zero, and overflow floating-point operations. All 30 annual outputs
were present. A scan of 6,167,640 stored numeric values across every output
field found no NaN or infinity. No runtime numerical failure was found in the
log.

The comparison uses F1 and the first 30 years of the matching Z1 f64 run. They
have the same source commit, tracked source patch, boundary input hashes,
resolution, layers, timestep, output cadence, and random seed. Annual grid
means below use Gaussian latitude weights. Years 11--30 reduce the influence of
the strongest initial adjustment, although neither trajectory is fully
equilibrated by year 30.

| Quantity, years 11--30 | f32 | f64 | f32 - f64 | f64 annual standard deviation |
| --- | ---: | ---: | ---: | ---: |
| Surface temperature | 283.822 K | 283.798 K | +0.025 K | 0.283 K |
| 2 m air temperature | 283.427 K | 283.405 K | +0.021 K | 0.269 K |
| Net TOA flux (shortwave + longwave) | -0.713 W m-2 | -0.652 W m-2 | -0.061 W m-2 | 0.337 W m-2 |
| Total precipitation | 2.833 mm day-1 | 2.825 mm day-1 | +0.008 mm day-1 | 0.018 mm day-1 |
| Snowfall | 0.281 mm day-1 | 0.278 mm day-1 | +0.003 mm day-1 | 0.008 mm day-1 |
| Cloud cover | 0.5677 | 0.5671 | +0.0005 | 0.0017 |
| Sea-ice fraction | 0.09977 | 0.09985 | -0.00008 | 0.00421 |
| Sea-ice thickness | 0.2449 m | 0.2448 m | +0.0001 m | 0.0184 m |

The years 11--30 climatological maps also remain close:

| Field | Spatial RMSE | Pattern correlation |
| --- | ---: | ---: |
| Surface temperature | 0.299 K | 0.999935 |
| 2 m air temperature | 0.276 K | 0.999924 |
| Net shortwave radiation | 1.270 W m-2 | 0.999908 |
| Net longwave radiation | 1.147 W m-2 | 0.999696 |
| Total precipitation | 0.125 mm day-1 | 0.998488 |
| Sea-ice fraction | 0.0086 | 0.999729 |

The annual maps do not match point by point. Surface-temperature RMS difference
is already about 1 K in the first annual mean and remains roughly 1--2 K. This
is expected after small arithmetic differences perturb a chaotic circulation;
climate statistics and conservation behavior are the useful comparison.

Both runs are still cooling and growing sea ice at year 30. Over years 11--30,
the surface-temperature trends are -0.56 K per decade in f32 and -0.44 K per
decade in f64. The final-decade f32 mean is only 0.044 K cooler than f64, but a
30-year transient cannot rule out a small precision-dependent long-term drift.
The 100-year f64 trajectory reaches a cooler later climate than either run's
years 11--30, confirming that the comparison window is still spin-up.

![Area-weighted f32 and f64 trajectories](reference/f32_vs_f64_climate.png)

## Performance observed on this CPU

F1 used 4,750.97 seconds for 30 years, or 158.4 seconds per model year. Z1 used
20,078.32 seconds for 100 years, or 200.8 seconds per model year. This is a
1.27x throughput increase, or about 21% less solver time. The final diagnostics
reported approximately 23.5 MB for f32 and 43.4 MB for f64.

These timings are useful for this compiler and machine but do not predict the
CUDA speedup.

## Implications for T170

F1 supports using f32 as a serious production candidate. It does not establish
that a T170L30 CPU or CUDA solver is stable or scientifically equivalent.
Several numerical pressures increase with resolution:

- Direct Legendre transforms accumulate over many more latitudes and spectral
  coefficients. The longest sums grow from roughly 22 terms at T21 to 171 at
  T170, while latitude reductions grow from 16 paired rows to 128.
- Thirty vertical layers introduce longer column reductions and larger
  semi-implicit systems than this ten-layer case.
- High-wavenumber coefficients can be small relative to the leading modes, so
  f32 loses their relative accuracy first.
- A higher-resolution integration normally needs a smaller timestep and
  resolution-appropriate diffusion. A stability failure caused by those
  settings could be mistaken for a precision failure.
- CUDA changes fused operations, transcendental implementations, and reduction
  order, so CPU f32 and CUDA f32 will not follow the same trajectory.

There is encouraging existing evidence in the implementation. `legmod.f90`
states that its transforms were tested in 32- and 64-bit builds from T21 to
T341. It generates the Gaussian grid and Legendre recurrence in explicit f64,
then stores the transform tables in the configured default precision. The
transform accumulations themselves use the default precision and therefore
remain the main T170 concern.

Before accepting T170 f32 for production:

1. Compare f32 and f64 transform round trips at T170, with error reported by
   total wavenumber and for derivatives as well as scalar transforms.
2. Run short paired T170L30 integrations and inspect finite values, extrema,
   conservation residuals, spectral tails, and sensitivity to timestep.
3. Run a multi-decade paired climate case with seasonal output and compare
   energy balance, water balance, trends, variability, and spatial spectra.
4. Repeat operation-level and climate validation for CUDA f64 and CUDA f32,
   because CPU f32 is not a substitute for validating GPU arithmetic.

If pure f32 is weak at T170, retain f32 state and physics while using f64 for
Legendre or global accumulations, conservation totals, and small implicit
solves. Those operations contain much less data-parallel work than the full
gridpoint physics, so selective wider accumulation should preserve most of the
GPU benefit.

## Evidence from documentation and the literature

The direct ExoPlaSim/PlaSim evidence is limited:

- The [ExoPlaSim documentation](https://github.com/alphaparrot/ExoPlaSim/blob/master/docs/index.rst)
  supports 4- and 8-byte builds, uses 8-byte precision by default, and warns
  that 4-byte precision may reduce stability.
- The [ExoPlaSim paper](https://arxiv.org/abs/2107.07685) uses an explicit
  `precision=8` setting in its worked THAI configuration. It does not report a
  controlled f32-versus-f64 validation experiment.
- The upstream [PlaSim quick-run instructions](https://gogs.elic.ucl.ac.be/TECLIM/PlaSim/src/5b126acd46267e15bcb69af4236c05d306b36703/README.md)
  tell users to enable double precision. This is a recommended configuration,
  rather than evidence that f32 is unsuitable.
- The header of [`legmod.f90`](https://github.com/alphaparrot/ExoPlaSim/blob/master/exoplasim/plasim/src/legmod.f90)
  says its transforms were tested at T21--T341 in 32- and 64-bit builds. This is
  encouraging for T170, but the comment gives no error thresholds, duration,
  compiler details, or climate-validation results. The same source computes
  Gaussian-grid and Legendre-table setup with explicit 8-byte intermediates,
  then stores and applies the tables at the configured model precision.

I did not find a paper that directly compares complete PlaSim or ExoPlaSim
climates in f32 and f64. The most relevant evidence therefore comes from other
spectral atmospheric models:

- [ECMWF's single-precision IFS study](https://www.ecmwf.int/en/newsletter/148/meteorology/single-precision-ifs)
  found comparable forecast skill and small one-year climate differences at
  TL399. It retained f64 for the one-time Legendre and vertical-operator setup,
  plus a few sensitive physics calculations. Long integrations required a
  pressure mass fixer. Single precision reduced step cost by 37%, rising to
  40.7% after loop-size tuning.
- [Paxton et al. (2022)](https://arxiv.org/abs/2104.15076) tested the
  coarse-resolution spectral SPEEDY climate model and found f32 more than
  sufficient under their statistical rounding-error test. This is strong
  supporting evidence for simplified spectral GCMs, but it does not establish
  T170 PlaSim behavior.
- [Kimpson et al. (2023)](https://arxiv.org/abs/2207.14598) ran 100-year SPEEDY
  climate-change experiments. Relative to f64, their f32 global-mean surface
  temperature and precipitation mean-bias errors were about 0.0002 K and
  0.00008 mm per six hours. Their result covers a long, changing climate, but
  still uses a different coarse-resolution model.
- [Cotronei and Slawig (2020)](https://arxiv.org/abs/2001.01214) converted ECHAM
  radiation incrementally. Most of it worked in f32, a small section retained
  higher precision, radiation sped up by about 40%, and the complete model sped
  up by up to 18%.

Taken together, the literature supports f32 as the main CUDA storage and
arithmetic format, with explicit f64 islands where tests show a need. It also
supports the present validation plan: test conservation and long-term climate,
not trajectory identity, and validate the high-resolution spectral transforms
separately. Our paired F1 result is more direct evidence for ExoPlaSim than any
published f32/f64 comparison I found.

## Assessment

No significant negative effect is visible in this T21L10 30-year experiment.
The f32 run is stable, its mean climate is very close to f64, and its spatial
climatology has very high pattern agreement. The result is encouraging enough
to continue with f32 and mixed-precision CUDA development. T170 requires its
own transform and climate validation before f32 becomes the default production
mode.
