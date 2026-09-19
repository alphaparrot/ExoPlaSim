# Z1 setup checks

These checks used the current CPU solver on 2026-09-17. They are short plumbing
checks outside `reference/`; they are **not** 100-year reference data.

| Configuration | Model duration | Solver result | Solver time |
| --- | --- | --- | ---: |
| T21 Earth, 45-minute step | 320 steps | Native output decoded and plotted | 6.06 s |
| T42 Earth, 45-minute step | 320 steps | Native output decoded and plotted | 27.28 s |
| T63 Earth, 45-minute step | 320 steps | SIGILL in native solver during early timesteps | — |
| T63 Earth, 30-minute step | 320 steps | Native output decoded and plotted | 71.82 s |
| T21 Earth, annual write cadence | 360 days | One decoded annual field record and plot | 212.29 s |
| T21 Earth, two short segments | 320 steps each | `MOST.00000` and `MOST.00001` decoded | — |

The T63 SIGILL was reproduced by running the compiled executable directly. Its
diagnostic shows that it read the 192×96 Earth boundary fields and entered the
first timesteps. The Fortran backtrace did not identify a source line, so the
precise floating-point operation is unresolved. A 30-minute step completed the
short check; the T63 Z manifests use that timestep and `NSTPW=17280` for one
annual mean record. Long-term stability at T63 remains unverified.
