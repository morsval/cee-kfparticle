# Migration to CeeV0Lite 0.6.1 / outputFormat 4

Version 0.6.1 does not change the nominal Lambda central-value estimator or the
default accepted sample. Every `legacy`, `loose`, and `star` preset still
requires all 24 outer numerical-Jacobian perturbation fits to succeed.

## Optional QA retention

The new CLI option

```text
--allow-incomplete-mother-jacobian
```

retains one narrowly defined failure class: the nominal POCA and vertex/TOF
central solution produced a usable mother state, but at least one of the twelve
`+h/-h` central-difference pairs could not be fitted. It never overrides a
nominal POCA or vertex-fit failure. By default `requireTofConvergence=true`, so
a non-converged TOF fixed point is also rejected. The separate diagnostic option
`--allow-nonconverged-tof` can relax only that TOF requirement; combining the
two diagnostic options can therefore retain a usable final TOF iterate whose
convergence threshold was not reached and whose mother Jacobian is incomplete.
Such rows are QA data, not the default physics sample.

To recover the default strict sample from a file produced with the option:

```text
motherJacobianComplete == 1
```

or equivalently:

```text
motherCovarianceValid == 1
```

Use `motherJacobianPerturbationFailed == 1` to isolate the retained QA rows.

## New branches

- `motherJacobianComplete`
- `motherJacobianPerturbationFailed`
- `motherJacobianFitsAttempted`
- `motherJacobianFitsSucceeded`
- `motherJacobianFailureCount`
- `motherJacobianPlusFailureMask`
- `motherJacobianMinusFailureMask`
- `motherJacobianFailedColumn`
- `motherJacobianFailedSide`
- `motherJacobianFirstFailureFitStatus`
- `covarianceDependentCutsEvaluated`

The twelve mask bits and column indices use

```text
0:x1 1:y1 2:z1 3:px1 4:py1 5:pz1
6:x2 7:y2 8:z2 9:px2 10:py2 11:pz2
```

`motherJacobianFailedSide` is `+1` or `-1`, and is zero when no side failed.
The attempted/succeeded counters are normally `24/24`.

## Invalid observables on a retained QA row

The following remain valid central observables: mass, SV, mother momentum and
energy, daughter DCA, geometric chi2, decay length, pointing, and daughter-PV
central metrics.

The final 7x7 mother covariance and every quantity that needs it are invalid and
written as `NaN`: `sigmaMassGeV`, SV errors, `sigmaDecayLengthCm`,
`decayLengthSignificance`, PV-constraint chi2, and topological chi2. The flag
`covarianceDependentCutsEvaluated` is zero because the online finder does not
apply length-significance or topology cuts to undefined values.

Never replace these `NaN` values with zero: zero would incorrectly mean an
infinitely precise measurement.
