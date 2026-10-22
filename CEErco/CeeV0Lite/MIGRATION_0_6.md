> 本页记录旧版本迁移；其中 Bz 表述已由 2026-09-20 By 补丁更正。当前 CEE 磁场沿 y，请同时阅读 `MIGRATION_BY_20260920.md`。

# Migrating from CeeV0Lite 0.5.x to 0.6.x

Version 0.6.0 intentionally changes two statistical meanings and the ROOT
schema. Do not merge 0.5.x and 0.6.x candidate trees under the same branch
definitions without an explicit conversion.

## 1. `topologicalChi2Ndf` is now the total

In 0.5.x, the public field called `topologicalChi2Ndf` contained only the
two-dimensional mother-to-PV constraint divided by 2. In 0.6.x it follows the
STAR-like accumulated fit convention:

```text
topologicalChi2    = geometricChi2 + primaryVertexConstraintChi2
topologicalNdf     = geometricNdf + 2
topologicalChi2Ndf = topologicalChi2 / topologicalNdf.
```

The former PV-only meaning is available explicitly as

```text
primaryVertexConstraintChi2
primaryVertexConstraintNdf        // 2
primaryVertexConstraintChi2Ndf
```

For an old QA plot that meant “PV constraint only”, switch from
`topologicalChi2Ndf` to `primaryVertexConstraintChi2Ndf`. For a STAR-like
candidate topology cut, keep `maximumTopologicalChi2Ndf`, but retune its value
against 0.6.x simulation because it now acts on the total. Even the PV-only
distribution need not be bitwise identical after the finite-field POCA change.

## 2. Daughter-to-PV cut uses raw chi2

The canonical 0.6.x API is

```cpp
cuts.minimumDaughterChi2FromPv = value;
```

and is applied to `Candidate::daughterChi2FromPv[i]`. The inspected STAR
`GetDeviationFromVertex()` value is raw chi2; the active STAR interface default
is `18.6`, which is the value used by the 0.6.x `star` preset.

The following legacy API is retained:

```cpp
cuts.minimumDaughterDeviationFromPv = value;
```

It still means `sqrt(chi2)`, as it did in 0.5.x. The old CLI
`--min-daughter-pv-deviation` has the same deprecated square-root meaning. New
steering should use

```text
--min-daughter-pv-chi2 VALUE
```

Do not pass `18.6` to the legacy square-root option. If an old steering value
was `d`, its equivalent raw threshold is `d*d` (subject to retuning and normal
floating-point boundary effects).

`Candidate::daughterDeviationFromPv[i]` remains as a deprecated alias of
`daughterPvSignificance[i] = sqrt(daughterChi2FromPv[i])` so source code can be
migrated gradually.

## 3. ROOT output format 3

`TNamed("CeeV0AnalysisConfig")` now records `outputFormat=3`. The candidate tree
retains format-2 mother-state/covariance fields and adds:

- `geometricChi2`, `geometricNdf`, `geometricChi2Ndf`;
- `topologicalChi2`, `topologicalNdf`, `topologicalChi2Ndf`;
- `primaryVertexConstraintChi2`, `primaryVertexConstraintNdf`,
  `primaryVertexConstraintChi2Ndf`;
- `daughterPvChi20`, `daughterPvChi21`;
- `daughterPvSignificance0`, `daughterPvSignificance1`.

The old ROOT names `daughterDeviation0/1` remain as deprecated aliases of
`daughterPvSignificance0/1`; they keep the square-root meaning and must not be
used for a raw-chi2 cut.

Use the `daughterPvChi2*` branches for STAR-like raw-chi2 selections. Use the
`daughterPvSignificance*` branches only for old square-root-significance plots.
Analysis scripts should read the metadata first and reject an unexpected format
instead of assuming that a branch with the same name has the same semantics.

## 4. Finite-field point-POCA

At non-zero uniform `Bz`, point-to-helix closest approach now uses exact
transverse periodicity to reduce the global search to at most two periods
around the longitudinal optimum, enumerates every stationary point in that
window, rejects local maxima, includes both original endpoints, and selects the
global distance minimum. This can move
the chosen path, fitted SV, and downstream metrics relative to 0.5.x when more
than one helix turn is available.

This change does **not** implement a simultaneous global two-helix DCA. It also
does not add a non-uniform field map, material effects, an extra pointing cut,
or a maximum-decay-length cut.

## Recommended migration check

1. Keep separate 0.5.x and 0.6.x output files.
2. Run the same simulated sample and compare accepted counts before tightening
   any topology cut.
3. Plot `geometricChi2Ndf`, `primaryVertexConstraintChi2Ndf`, and the new total
   `topologicalChi2Ndf` separately.
4. Scan raw daughter thresholds using `daughterPvChi20/1`; do not infer the new
   optimum only by squaring an old cut.
5. For finite-field samples, inspect chosen path/vertex residuals in the
   multi-turn and path-boundary phase space.
6. Record `outputFormat=3`, cut preset, and all overrides with the physics
   result.
