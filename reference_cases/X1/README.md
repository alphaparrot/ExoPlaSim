# X1 — Passive tracer transport

## What this case is

if generic tracers remain in scope, initialize a nonuniform tracer and verify advection, positivity and global mass behavior.

## What it tests

**Blocked:** The compiled solver sets NTRACE=1 and reserves it for humidity; a separate passive tracer cannot be initialized or checked without a solver change.

No reference run exists for this feature yet.
