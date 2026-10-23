# 0.7.0 — configurable analysis layer

Added CeeKFParticleAnalysisMaker, Interface, format-neutral events, a CEE ROOT reader, selectable/custom ROOT columns, role-aware daughter DCA accessors and a direct-call macro. Core fit equations and the v0.6.5 default magnetic track–PV DCA reference are unchanged. See ../docs/AnalysisMaker_v0.7.0.md.

# Changelog

## 0.6.3 — 2026-09-30

- Added explicit-charge/By production-vertex fitting for charged mothers with a free arc-length S.
- Added momentum rotation and full S derivatives in both endpoint transports; preserve charge/By in the state.
- Retained the neutral fit and guarded the old overload against omitted fields for charged mothers.
- Use a square-root covariance transform for stable semidefinite endpoint transport.
- ROOT outputFormat=6 adds motherTopoCharge and motherTopoByTesla; default Lambda/V0 channels and cuts remain unchanged.
- Added charged GLS/ODE tests and an integration example for hypertriton -> helium-3 + pion using the formal C++ core.
- See MIGRATION_0_6_3.md for local-flight limits, rigidity conventions, missing upstream correlations and remaining experiment validation.


## 0.6.2 — 2026-09-28

- Add a separate neutral 8D production-vertex constrained state, S and the full covariance.
- Fit a correlated 10D mother/PV prior; re-linearize against one fixed prior.
- Propagate the full TOF solve's PV response into both marginal and cross covariance.
- Add signed flight length, cτ, correlated errors and idempotent endpoint transport.
- Preserve original seven-dimensional mass/cuts; invalid auxiliary states are flagged/NaN.
- ROOT outputFormat=5 adds motherTopo* branches; --no-production-state disables the extra fit.
- Charged composites and cascade reconstruction remain unsupported. See MIGRATION_0_6_2.md.


## 0.6.1 CEE By correction (2026-09-20)

- Corrected the CEE field axis to y; the beam remains along z. The previous
  0.6.1 implementation incorrectly used a uniform z field.
- Corrected x-z helix transport, the Lorentz-force sign, point-POCA derivatives,
  circle centres, and the y-parallel multi-turn search interval.
- Renamed the CEE API to `magneticFieldByTesla` and CLI to `--by-tesla`.
  Legacy `--bz-tesla 0` remains a warned zero-field compatibility option;
  nonzero `--bz-tesla` is rejected. State/covariance ordering is unchanged.
- Output metadata records `ByTesla`, field axis y, beam axis z and this patch
  identifier. Added independent Lorentz-force ODE/covariance and CLI tests.
- Corrected the magnetic-field pages of the existing 78-page Chinese manual.
  Version remains 0.6.1 and ROOT outputFormat remains 4.

See `MIGRATION_BY_20260920.md` for API and run-command changes.

## 0.6.1

- Added an explicit, opt-in QA path for candidates whose nominal POCA,
  vertex/TOF solution and mother central values succeed but at least one of the
  24 outer numerical-Jacobian perturbation fits fails.  Strict rejection
  remains the default for every cut preset.
- Added `kMotherCovarianceJacobianPerturbationFailed`, 24-side attempted/
  succeeded/failure counts, first failed column/side/status, and separate 12-bit
  plus/minus failure masks.
- Added `--allow-incomplete-mother-jacobian`. Retained QA rows keep central
  observables that do not require the final 7x7 mother covariance; the final
  covariance, its errors, length significance and mother-covariance-dependent
  topology statistics are `NaN`. Geometric and daughter-PV chi2 still use their
  input track/PV covariances and are not generically “covariance-independent.”
- Added ROOT branches for the complete Jacobian diagnostics and bumped
  `CeeV0Analysis` metadata to `outputFormat=4`.
- Added deterministic unit and end-to-end tests for strict rejection, explicit
  retention, all-24-side accounting and `NaN` invalidation.

See `MIGRATION_0_6_1.md` before enabling the QA retention option.

## 0.6.0

- Reworked finite-`Bz` point-to-helix closest approach: exact transverse
  periodicity reduces the global search to at most two periods around the
  longitudinal optimum; every stationary point there is considered, local
  maxima are rejected, and the original interval endpoints are also compared.
- Changed `Candidate::topologicalChi2Ndf` to the STAR-like total definition
  `(geometricChi2 + primaryVertexConstraintChi2) /
  (geometricNdf + 2)`.  The PV-only contribution is now exposed separately as
  `primaryVertexConstraintChi2{,Ndf}` for QA and migration studies.
- Added the canonical raw daughter-to-PV cut
  `minimumDaughterChi2FromPv`; the `star` preset uses the active STAR interface
  value `18.6`.  The old `minimumDaughterDeviationFromPv` and
  `--min-daughter-pv-deviation` retain their historical `sqrt(chi2)` meaning as
  deprecated compatibility controls.
- Bumped `CeeV0Analysis` ROOT metadata to `outputFormat=3` and added explicit
  raw chi2, NDF, PV-only topology, and daughter raw-chi2/significance branches.
- Added regression checks for the topology identities, raw-versus-square-root
  daughter cuts, cut boundaries, and finite-field point-POCA minimum selection.

See `MIGRATION_0_6.md` before comparing 0.6.x candidate trees or cut values
with 0.5.x output.

## 0.5.1

- Made ordinary simulation use the exact generator time origin `t=0` without
  checking, binding, or reading `CeeT0Hit`.
- Changed the lightweight-core defaults to zero start-time uncertainty, zero
  additional common-time uncertainty, and disabled shared-T0 propagation.
- Kept the independent iTOF/eTOF hit-time resolution active in simulation.
- Made real-data mode explicitly enable measured-T0 input and common-time
  covariance, subject to the user-supplied time-unit conversion and calibration.
- Centralized mode-dependent defaults in `ConfigureDataMode()` and reject
  measured-hit-only T0 options when the selected source is fixed.
- Added configuration-policy and end-to-end regression tests that prevent a
  shared T0 uncertainty from being silently reintroduced into simulation.

## 0.5.0

- Added the seven-parameter mother state `(x,y,z,px,py,pz,E)` and a stored dense
  7x7 covariance under the implemented input-error model.
- Replaced the order-dependent sequential vertex update with a symmetric batch
  fit of the two transverse track measurements.
- Eliminated the two free path parameters before applying the single effective
  two-track DCA residual (`ndf=1`).
- Made closest-approach covariance propagation recompute the path for every
  state perturbation.
- Added accepted-candidate `12D -> 7D` numerical covariance propagation through
  the currently implemented POCA/vertex estimator.
- Replaced independent daughter mass-error addition with propagation from the
  mother `(P,E)` covariance block.
- Removed event T0 from the two daughter marginal covariances and added it once
  as a shared fixed-point nuisance; added full-step/half-step stability checks.
- Included mother momentum and SV--momentum covariance in topological chi2.
- Added SPD/PSD checks, explicit covariance status values and output-format-2
  ROOT branches including packed `motherCovariance[28]`.
- Added daughter-order, reference-point, finite-field, exact shared-T0 rank-one,
  indefinite-covariance and small-error toy-MC regression tests.

Known model boundaries are documented in
`docs/STAR_LIKE_MOTHER_COVARIANCE.md`.
