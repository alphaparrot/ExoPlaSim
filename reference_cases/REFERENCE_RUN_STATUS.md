# Reference capture status

Status on 2026-09-17. The D, P, A, S, X, and I groups have been attempted. The Z long-run and F four-byte groups have manifests but no captured reference yet. Each successful variant has its original solver timing, provenance, plots, diagnostics, and native field data under its ignored `reference/` directory.

| Group | Captured | Planned variants |
| --- | ---: | ---: |
| D | 16 | 17 |
| P | 9 | 12 |
| A | 18 | 18 |
| S | 20 | 23 |
| X | 4 | 7 |
| I | 7 | 7 |
| Z | 0 | 9 |
| F | 0 | 2 |

**74 of 86 non-Z variants are captured.** X1 is separately blocked because the active solver has no independent passive-tracer slot.

| Case | Captured | Planned | Missing variants |
| --- | ---: | ---: | --- |
| D1 | 1 | 1 | — |
| D2 | 3 | 3 | — |
| D3 | 2 | 2 | — |
| D4 | 1 | 2 | [`mpi2`](./D4/reference/mpi2.log) |
| D5 | 2 | 2 | — |
| D6 | 3 | 3 | — |
| D7 | 2 | 2 | — |
| D8 | 2 | 2 | — |
| P1 | 1 | 1 | — |
| P2 | 2 | 2 | — |
| P3 | 3 | 3 | — |
| P4 | 1 | 2 | [`fixed_orbit`](./P4/reference/fixed_orbit.log) |
| P5 | 2 | 3 | [`wolf_spectrum`](./P5/reference/wolf_spectrum.log) |
| P6 | 0 | 1 | [`mars`](./P6/reference/mars.log) |
| A1 | 1 | 1 | — |
| A2 | 3 | 3 | — |
| A3 | 3 | 3 | — |
| A4 | 2 | 2 | — |
| A5 | 3 | 3 | — |
| A6 | 1 | 1 | — |
| A7 | 1 | 1 | — |
| A8 | 2 | 2 | — |
| A9 | 2 | 2 | — |
| S1 | 1 | 1 | — |
| S2 | 3 | 3 | — |
| S3 | 1 | 1 | — |
| S4 | 1 | 1 | — |
| S5 | 2 | 2 | — |
| S6 | 2 | 2 | — |
| S7 | 2 | 2 | — |
| S8 | 2 | 2 | — |
| S9 | 2 | 2 | — |
| S10 | 2 | 3 | [`diagnostic_simba`](./S10/reference/diagnostic_simba.log) |
| S11 | 0 | 1 | [`initial_ice_sheet`](./S11/reference/initial_ice_sheet.log) |
| S12 | 2 | 3 | [`supply_limited_evolving_co2`](./S12/reference/supply_limited_evolving_co2.log) |
| X1 | 0 | 0 | The compiled solver sets NTRACE=1 and reserves it for humidity; a separate passive tracer cannot be initialized or checked without a solver change. |
| X2 | 0 | 1 | [`aerosol_transport`](./X2/reference/aerosol_transport.log) |
| X3 | 1 | 2 | [`transport_only`](./X3/reference/transport_only.log) |
| X4 | 1 | 1 | — |
| X5 | 0 | 1 | [`permissive_capture`](./X5/reference/permissive_capture.log) |
| X6 | 2 | 2 | — |
| I1 | 2 | 2 | — |
| I2 | 1 | 1 | — |
| I3 | 2 | 2 | — |
| I4 | 1 | 1 | — |
| I5 | 1 | 1 | — |
| Z1 | 0 | 3 | `t21_100y` (pending), `t42_100y` (pending), `t63_100y` (pending) |
| Z2 | 0 | 3 | `t21_10y` (pending), `t42_10y` (pending), `t63_10y` (pending) |
| Z3 | 0 | 3 | `t21_10y` (pending), `t42_10y` (pending), `t63_10y` (pending) |
| F1 | 0 | 1 | `t21_f32_100y` (pending) |
| F2 | 0 | 1 | `t42_f32_30y` (pending) |

A missing variant with a log was attempted but did not complete validation. Read its log and, where present, the corresponding `ref_*_crashed/` directory before changing inputs or the solver. Z variants have not been run to their declared duration; short Z smoke checks are described in [README.md](README.md). F variants have not been run at all.
