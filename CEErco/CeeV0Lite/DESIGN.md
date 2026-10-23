# Mathematical notes

The fitted daughter state is

```text
x = (x, y, z, px, py, pz),  covariance C[6][6].
```

`CeeV0LiteTransport` propagates the parameters to a signed path length `s`.
CEE coordinates have beam along z and uniform field `B=(0,By,0)`.
Momentum magnitude and `py` stay constant while `(px,pz)` rotate. For
`theta=omega*s`, positive charge moving along +z bends toward -x in +By:

```text
omega = q * 0.00299792458 * By[T] / p[GeV/c]   [rad/cm]
px(s) = px*cos(theta) - pz*sin(theta)
py(s) = py
pz(s) = pz*cos(theta) + px*sin(theta)
x(s) = x0 + (px*sin(theta) - pz*(1-cos(theta))) / (p*omega)
y(s) = y0 + s*py/p
z(s) = z0 + (pz*sin(theta) + px*(1-cos(theta))) / (p*omega)
```

For a fixed path, the state covariance is propagated with a central-difference
Jacobian.  `TransportToClosestApproach` uses a stronger Jacobian: every perturbed
state re-solves the closest path, so `ds_closest/dx` is included.  This is needed
for invariance under moving a reference state along the same trajectory.
Optional position and angular process noise are exposed in `TransportConfig`;
the final outer mother-covariance propagation currently requires both noise
rates to remain zero because their independent stochastic modes are not yet
part of that outer Jacobian.

For closest approach of one helix to a fixed point at finite `By`, the distance
squared is treated as a bounded one-dimensional function of path.  Exact
periodicity in the x-z plane proves that a global minimizer is represented within
one period of the field-parallel (y) optimum (a window at most two periods wide), or
at an original interval endpoint.  The solver enumerates every stationary
point in that reduced window, rejects local maxima, evaluates both original
endpoints, and selects the smallest distance.  This is important because a helix can provide several stationary
paths and a single Newton start can converge to the wrong turn or to a local
maximum.  This point-to-helix operation must not be described as a simultaneous
global minimization over both path parameters of two helices.

At each vertex iteration, both tracks are transported to their closest approach
to the same frozen current vertex. A track constrains only the two coordinates
transverse to its local direction.  For each track, with transverse basis `Hi`,
transported position `mi` and projected covariance `Ri`, one symmetric batch
update solves

```text
N = sum_i Hi^T Ri^-1 Hi
b = sum_i Hi^T Ri^-1 Hi mi
v = N^-1 b,     Cov(v)=N^-1.
```

Both daughters are linearized before this solve, so swapping their order does
not change the result.  After convergence, the two free path parameters leave
one effective two-track residual.  With local unit directions `u1,u2`, define
`n=unit(u1 cross u2)` and constrain only

```text
x = (r1, p1, r2, p2),  Cx = diag(C1,C2),
q = n dot (r2-r1),
H = (-n,0,+n,0),
S = H Cx H^T + R,
x_fit = x - Cx H^T S^-1 q,
C_fit = Cx - Cx H^T S^-1 H Cx.
```

This one-dimensional residual is consistent with `ndf=1`; constraining all three
raw position components would incorrectly treat the two path nuisances as
measurements. `C_fit` supplies correlated fitted daughter marginals and the
central mother state

```text
y = (SVx, SVy, SVz, Px, Py, Pz, E),
P = p1+p2,
E = sqrt(p1^2+m1^2) + sqrt(p2^2+m2^2),
```

For the covariance saved in an accepted candidate, the code does not rely on
the fixed-linearization `C_fit`.  It treats the currently implemented POCA plus batch-vertex
estimator as a map `F: R^12 -> R^7`; every input position/momentum component is
perturbed and the complete solve is repeated:

```text
J = dF/d(r1,p1,r2,p2),
Cy = J diag(C1,C2) J^T.
```

Here `C_proton,pion` would mean the 6x6 cross-covariance
`Cov((r_proton,p_proton),(r_pion,p_pion))`; it is not a scalar, an interaction
term, or a named STAR class member. The input model assumes this block is zero
upstream. The local `C_fit` develops an effective non-zero block from the common
vertex residual, while the production calculation compresses all correlations
induced by the current closest-path/common-SV map into one `Cy` without exposing
or storing that block. This is a scalar reference implementation with a
high-level data flow analogous to KFParticle; it is deliberately slower than
STAR's analytic SIMD Jacobians.
Mass and mass uncertainty are calculated only from the mother `(P,E)` block,

```text
M^2 = E^2-P^2,
dM/d(P,E) = (-Px/M,-Py/M,-Pz/M,E/M),
sigmaM^2 = g Cy g^T.
```

For real data, the two TOF estimates share event T0. `sigmaStartTimeNs` is
therefore excluded from both independent daughter covariances. The complete TOF fixed-point is
rerun at `T0+/-h` and `T0+/-h/2`.  The two derivatives must agree within
`sharedTimeDerivativeStabilityTolerance`; with the half-step result `b=dy/dT0`,
the single common contribution is

```text
Cy <- Cy + b (sigmaT0^2 + sigmaDecayTimeModel^2) b^T.
```

This includes the T0-induced daughter correlation without double counting it.
Ordinary simulation instead uses the exact generator origin `t=0`, both common
sigmas equal to zero, and shared-time propagation disabled. Its independent
stop-hit timing resolution is still propagated per daughter.
`daughterDcaCm` remains the separation of the two states transported to their
respective point-to-helix closest approaches to the fitted vertex; it is not a
separate simultaneous global two-helix minimization.

The topological pointing residual uses the full mother `(SV,P)` covariance,
including the SV--momentum cross block and the direction derivative proportional
to decay length.  Let its two-dimensional PV-constraint contribution be
`chi2PV`, with `ndfPV=2`, and let the common-vertex fit provide `chi2Geo` and
`ndfGeo`.  Candidate-level topology uses the STAR-like accumulated definition

```text
topologicalChi2    = chi2Geo + chi2PV,
topologicalNdf     = ndfGeo + 2,
topologicalChi2Ndf = topologicalChi2 / topologicalNdf.
```

The PV-only quantities remain available as
`primaryVertexConstraintChi2`, `primaryVertexConstraintNdf`, and
`primaryVertexConstraintChi2Ndf`; they are QA components, not synonyms for the
total.  The canonical daughter-to-PV metric is raw chi2
`daughterChi2FromPv[i]`, so the production cut is
`minimumDaughterChi2FromPv`.  `daughterPvSignificance[i]` and the deprecated
`minimumDaughterDeviationFromPv` preserve the historical `sqrt(chi2)` convention
for compatibility only.

The compact output still does not retain one full joint
covariance containing
PV, SV, both daughter states and T0.  Upstream track--track/PV--track correlations
which are absent from the input cannot be reconstructed, and the daughter-to-PV
metrics do not carry the shared-T0 cross block.  In TOF mode, the base 12D
Jacobian begins at the converged track states; only the common-time nuisance is
currently rerun through the entire TOF fixed point.  PV covariance feedback into
the TOF momentum is also not yet propagated.  Therefore simulation pull/coverage
is still required rather than copied STAR cut values.

Closest-approach, vertex and zero-field TOF fixed-point iterations expose their
convergence state. Production defaults reject a path-limit solution, degenerate
parallel/anti-parallel pair, or a candidate that reaches an iteration limit. The
non-strict switches exist only for diagnostic studies.

This is deliberately a two-prong vertex Kalman fitter, not a full drop-in
replacement for every KFParticle decay-chain feature. The extension seams are:

- another magnetic-field transport implementation;
- material-dependent process noise;
- more PID hypotheses/channels in Finder;
- a multi-daughter or cascade fitter for Xi/Omega;
- detector-specific adapters outside the mathematical core.

Non-uniform field maps, material transport, a simultaneous two-helix global
DCA, and additional pointing/maximum-decay-length selection are not supplied by
the finite-`By` point-POCA change above.

For the current CEE zero-field Lambda analysis, each input track also carries a
source-dependent raw PV-DCA boundary. The default role assignment accepts the
(anti)proton hypothesis below or at that boundary and the pion hypothesis above
it. This reproduces the existing AnaMaker selection and is deliberately kept
separate from the covariance-normalized PV-deviation topology cut.

## v0.6.2 neutral production state

`FitProductionVertex` fits the joint prior `(SV,p,E,PV)` with two smooth pointing
constraints. With fixed initial unit direction n and transverse bases e_a,
`g_a = e_a.d - (n.d)/(n.p) e_a.p`, where d=SV-PV. Its analytic derivative uses
`w_a=e_a - (e_a.p)/(n.p)n`: H_a=(w_a,-S_plane*w_a,0,-w_a).
Every iteration applies `x_next=x0-C H^T(H C H^T)^-1 [g(x)+H(x0-x)]` to the
same x0,C. No repeated conditioning is performed. At convergence the covariance
is `(I-KH) C (I-KH)^T`, and S=p.d/p² is mapped with its full endpoint derivative.
The state is then transported to production using r_prod=r_decay-S*p; the
Jacobian includes d(r_prod)/dS=-p. All eight rows/columns are propagated.

Finder builds the TOF joint input as `[[C_cond+D V D^T,D V],[V D^T,V]]`; D comes
from full SolvePair PV perturbations, checked at two step sizes. This accounts
for algorithm-induced PV correlation under independent raw track/PV/T0 inputs,
not unknown upstream correlations. The existing conditional track covariance
model is retained. See MIGRATION_0_6_2.md for assumptions and output semantics.

The new state is auxiliary: legacy cuts and seven-dimensional outputs are kept.
This is a neutral CEE joint fit, not a numerical clone of STAR.


## v0.6.3 charged production state

For charge q and uniform CEE By, impose g=rD-rP-A(S)pD=0, where S is arc length / momentum.
A(S) is the integrated backward helix and R=dA/dS rotates decay momentum to production.
Linearizing gives H dx-v dS+g=0, v=R pD. Project the linearized system onto two directions
perpendicular to v, apply the constrained update against the same joint 10D prior, and recover dS
from the remaining longitudinal equation. Only the projected 2x2 covariance is inverted.

The full prior-to-posterior derivative is Jx=I-K H_transverse. The flight row is
Js=v^T H Jx/(v^T v). Produce rP from the fitted PV rows directly, and transform momentum with
R Jx_p+(dR/dS pD) Js. This supports zero-covariance known PVs without artificial jitter.
Endpoint covariance transport uses the Gram form (J L)(J L)^T, C=L L^T.
See MIGRATION_0_6_3.md for matrices, explicit API and local-flight validity limits.
