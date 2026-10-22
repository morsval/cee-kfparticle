# v0.7.0: start here

New analyses should use `macros/runCeeKFParticle.C`, which directly calls `CeeKFParticleAnalysisMaker`. See [the Chinese guide](../docs/AnalysisMaker_v0.7.0.md) and [candidate variables](../docs/CandidateVariables_v0.7.0.md). Build the shared library with `make shared` in a configured CEE environment.

The older Analyzer/CLI and runCeeV0Reco.C remain available for the fixed Lambda tree and legacy QA/truth histograms described below. They are a compatibility path, not aliases of the new configurable writer.

# CEE lightweight Lambda analysis

**v0.6.5（2026-10-22）**：有磁模式默认将每条输入 track 在内存中转换到与当前重建 PV 的 DCA 参考点，完整输运位置、动量与 6×6 协方差；配对前只做一次。详见根目录 `MIGRATION_0_6_5.md` 与 `v0.6.5_review`。切换 `--magnetic-track-reference stored` 可重放旧参考约定。零场 TOF 不变。更新后需要清理并重编核心与分析程序。

**v0.6.3（2026-09-30）**：正式 C++ 核心新增均匀 By 下带电母粒子的八维 PV 约束，
同时处理位置、动量旋转及 S 的完整相关误差。保留原中性分支与 V0 选择。
ROOT 输出格式为 6，新增 motherTopoCharge 和 motherTopoByTesla。
超氚集成示例位于 CeeV0Lite/examples/hypertriton_by；默认分析入口仍为原 Λ/V0 通道。
接口、边界及重新编译步骤见 `../CeeV0Lite/MIGRATION_0_6_3.md`。



> 2026-09-20 磁场方向修正：CEE 束流沿 z，均匀磁场沿 y；使用带符号 `--by-tesla` / `magneticFieldByTesla`。详见 `../CeeV0Lite/MIGRATION_BY_20260920.md`。
This directory is a complete analysis layer for `CeeV0Lite`. It is independent
of the old `MyTrack/MyLambda/AnaMaker` reconstruction and never uses MC truth as
a reconstruction selection.

## What it does

- reads a `ceesim` chain from one ROOT file, a `.list`, or a directory;
- reads `CeeTpcTrack`, `CeeTpcVertex`, `CeeiTOFHit`, `CeeEtofHit`; `CeeT0Hit`
  is checked/bound/read only in real-data measured-T0 mode;
- uses the highest-rank valid primary vertex and its stored errors;
- atomically validates/selects the `CeeTpcTrack` vertex-seed reference state when
  valid, with a visible QA fallback to legacy `GetX/Y/Z/GetP*` for old
  zero-field files. Magnetic tracking retains the full 6x6 seed covariance;
  zero-field TOF uses its 3x3 position block plus a compressed scalar direction
  error, not its anisotropic direction/cross blocks;
- performs maximum-cardinality, minimum-normalized-residual, event-level unique
  track-to-TOF assignment. One track and one physical hit are each used at most
  once;
- runs Lambda-only, TPC-TPC `CeeV0Lite::Finder` reconstruction;
- writes one `LambdaCandidates` entry per candidate plus basic QA histograms;
- in simulation mode optionally labels candidates from `MCTrack`; real mode
  never enables, binds or reads the `MCTrack` branch.

In magnetic tracking mode the analysis intentionally does not attach TOF with
the zero-field straight-line projector. Tracking supplies momentum; a future
TOF cross-check in field must use a helical detector-surface matcher.

Simulation output also contains generated-Lambda count, pT, and rapidity QA
histograms, so the first closure has a truth denominator without letting truth
enter reconstruction or selection.

The candidate tree retains the assigned TOF detector (`1=iTOF, 2=eTOF`), hit
index, residual, raw/converted time, all main vertex/cut variables, and exact
truth role checks. Output format 4 retains the format-3 `energyGeV`,
`motherCovarianceValid`, `motherCovarianceStatus`, and the 28-element packed
lower triangle of the mother covariance ordered as `(x,y,z,px,py,pz,E)`, and
all explicit metric components. It additionally writes:

- `motherJacobianComplete` and the explicit alias
  `motherJacobianPerturbationFailed`;
- `motherJacobianFitsAttempted`, `motherJacobianFitsSucceeded`, and
  `motherJacobianFailureCount`;
- `motherJacobianPlusFailureMask` and `motherJacobianMinusFailureMask`, whose
  lower 12 bits follow `(x1,y1,z1,px1,py1,pz1,x2,y2,z2,px2,py2,pz2)`;
- `motherJacobianFailedColumn`, `motherJacobianFailedSide`, and
  `motherJacobianFirstFailureFitStatus` for the first failed side;
- `covarianceDependentCutsEvaluated`.

The metric components retained from format 3 are:

- `geometricChi2`, `geometricNdf`, `geometricChi2Ndf`;
- `topologicalChi2`, `topologicalNdf`, `topologicalChi2Ndf`, where the last is
  `(geometricChi2+primaryVertexConstraintChi2)/(geometricNdf+2)`;
- `primaryVertexConstraintChi2`, `primaryVertexConstraintNdf`,
  `primaryVertexConstraintChi2Ndf` for the PV-only QA component;
- `daughterPvChi20/1` as raw daughter-to-PV chi2 and
`daughterPvSignificance0/1` as `sqrt(chi2)`.

For source-file compatibility, `daughterDeviation0/1` are also written as
deprecated aliases of `daughterPvSignificance0/1`; they are not raw chi2.

Every preset still requires a complete 24-side mother Jacobian by default, so
ordinary output has `motherJacobianComplete=1` and
`motherCovarianceValid=1`. The explicit diagnostic option
`--allow-incomplete-mother-jacobian` retains only the exact case where the
nominal candidate succeeded but at least one outer `+/-h` refit failed. Such a
row has `motherJacobianPerturbationFailed=1`; its final covariance,
`sigmaMassGeV`, SV errors, `sigmaDecayLengthCm`, `decayLengthSignificance`,
PV-constraint statistics and topological statistics are `NaN`. Central mass,
SV, momentum, daughter DCA, geometric chi2, decay length and pointing remain
valid. Recover the strict sample offline with
`motherJacobianComplete==1` (equivalently `motherCovarianceValid==1`). Other
covariance failures remain hard rejections. The format number is stored as
`outputFormat=6` in `TNamed("CeeV0AnalysisConfig")`, not as a tree branch.
`truthMatched=1` requires daughter 0 to be MC proton,
daughter 1 to be MC pi-, and both to share a PDG 3122 mother. `-1` means truth
could not be resolved and is not counted as background.

## Required installation order

Before compiling this analysis, merge the repository versions of
`CeeTpcTrack`, `CeeTpcVertex`, and `CeeT0Hit` into the real `VMCWORKDIR` source,
then rebuild their CEE libraries and ROOT dictionaries. Merely copying this
analysis directory is insufficient: calls such as `HasVertexSeed()` and
`GetVertexSeedCov()` otherwise produce missing-member or undefined-symbol
errors.

Old ROOT files remain readable after a correct ROOT schema evolution build, but
new fields will normally contain no vertex seed and empty vertex-error vectors.
The adapter then uses configured zero-field fallbacks. Check
`hTrackReferenceSource`, `hReferenceAdapterStatus`, and `hPvStoredErrors`; do not
assume the full 6x6 seed covariance entered a zero-field TOF state merely
because the program ran.

For future magnetic-field data the tracking producer must call
`CeeTpcTrack::SetVertexSeed()` with a consistent `(x,y,z,px,py,pz)` covariance,
store fitted charge, and fill vertex errors before producing new ROOT files.
Magnetic mode deliberately refuses the legacy fallback.

## Build

```bash
cd CEErco/CeeV0Analysis
export VMCWORKDIR=/path/to/your/ceeroot
export FAIRROOTPATH=/path/to/fairroot
make check
make -j4
```

The Makefile uses the same CEE/FairRoot libraries as the existing analysis.
If your local library names differ, use the overrides below or adjust only the
include/library section of this Makefile.
For a site-specific FairSoft include or extra library, prefer an override such
as `make EXTRA_INCLUDES=-I/software/apps/FairSoft/include` or
`make EXTRA_LDLIBS='-lYourLibrary'`.

## First simulation closure

The recommended first run uses the `legacy` preset and no tight topological
cuts. Simulation mode automatically uses the exact generator time origin
`t=0`, with no `CeeT0Hit` and no shared T0 uncertainty:

```bash
./CeeV0Analysis.bin \
  --input /data/cee_sim.root \
  --output lambda_legacy.root \
  --mode simulation \
  --momentum tof --by-tesla 0 \
  --cut-preset legacy --events 10000 --debug 1
```

Input may instead be `files.list` (one ROOT path per line) or a directory.
The default TOF raw-time scale is `1`; change it with
`--tof-time-scale-to-ns` if the persisted unit is not ns.

To make the central-success/Jacobian-failure population available for an
offline QA cut, add:

```bash
--allow-incomplete-mother-jacobian
```

Do not use its `NaN` covariance-dependent branches in physics cuts. The
standard run above intentionally omits this option and keeps strict behavior.

After the closure works, compare:

```bash
./CeeV0Analysis.bin --input files.list --output lambda_loose.root \
  --mode simulation --cut-preset loose

./CeeV0Analysis.bin --input files.list --output lambda_star_like.root \
  --mode simulation --cut-preset star
```

`star` is only a STAR-like approximate mapping. Its topology now uses the same
accumulation pattern as STAR—geometric fit chi2/NDF plus the two-dimensional PV
production constraint—and its canonical daughter cut is the raw chi2
`minimumDaughterChi2FromPv=18.6`, matching the active STAR interface default.
The reduced covariance/input model and decay-length definition are still not
statistically identical to STAR KFParticle, and STAR does not apply every mapped
quantity as one uniform hard Lambda cut. The preset does not silently add an
extra pointing or maximum-decay-length cut.

The zero-field proton/pion role assignment defaults to the existing CEE rule:
proton candidate DCA <= 1.2 cm and pion candidate DCA > 1.2 cm. It can be
changed with `--lambda-role-dca-cm` or disabled for a control study.
With no measured charge this labels a `p + pi-` mass hypothesis; topology alone
cannot experimentally distinguish Lambda from anti-Lambda charge conjugates.

For an individual raw daughter-PV cut use
`--min-daughter-pv-chi2 VALUE`. The old
`--min-daughter-pv-deviation VALUE` is retained for compatibility and still
means `sqrt(chi2)`; it is deprecated and is not where the STAR value `18.6`
belongs. `--max-topo-chi2ndf` applies to the total topology quantity above;
use the `primaryVertexConstraint*` output branches when studying only the PV
component. Because these are semantic changes from 0.5.x, do not mix old and
new candidate trees without following `../CeeV0Lite/MIGRATION_0_6.md`.

The fallback uncertainties used only when old files lack a valid vertex seed
or stored PV errors are also CLI options:
`--fallback-track-sigma-{x,y,z}-cm`,
`--fallback-direction-sigma-rad`, and
`--fallback-pv-sigma-{x,y,z}-cm`. Their defaults are placeholders, not a CEE
calibration. Use `--disable-legacy-track-fallback` and `--strict-pv-errors` to
turn either fallback into a hard input-integrity check. `--max-tracks` and
`--max-candidates` expose the two event-size safety limits; increase them if
the corresponding output metadata shows a limit inappropriate for your run.

## T0 policy

The hit class provides `GetTime()`, but the source is not sufficient to settle
the persisted unit/calibration: the original header comment said ps while
existing numerical usage looks compatible with ns. Therefore:

- simulation is strict: generator `t=0`, `sigma(T0)=0`, additional common-time
  sigma `=0`, and shared-T0 propagation disabled;
- simulation does not check, bind, or read a `CeeT0Hit` branch. A numerical
  `startTimeNs=0` is only the coordinate origin needed by the TOF equation;
- the independent iTOF/eTOF hit resolution `--tof-time-sigma-ns` remains active;
- `--t0-source hit`, a non-zero fixed/T0/common-time value, or a T0 scale is
  rejected in simulation rather than silently contaminating MC closure;
- hit-only scale, offset, and multi-hit policy options are rejected whenever
  the selected source is fixed, because they would otherwise be silently
  ignored;
- real-data hit mode requires an explicit `--t0-scale-to-ns`, even if it is
  `1`;
- conversion is `t0_ns = raw * scale + offset`;
- `first`, `earliest`, `mean`, and default `median` multi-hit policies exist;
- inspect `hT0Raw`, `hT0`, and `hTofMinusT0` before interpreting momentum.

The real-mode defaults `sigma(T0)=0.10 ns` and
`additionalDecayTimeSigmaNs=0.10 ns` are placeholders retained for the first
integration, not CEE calibration results.  Because these are independent
common terms, their default combined width is
`sqrt(0.10^2+0.10^2)=0.141 ns`.  Set the additional term to zero when there is
no separately justified common decay-time-model uncertainty, and replace the
T0 resolution from calibration before physics use. Programmatic callers should
call `ConfigureDataMode(config, kRealData)` before applying their overrides;
the CLI does this automatically.

Example only (the `1e-3` must be confirmed from CEE calibration):

```bash
./CeeV0Analysis.bin --input real.list --output lambda_real.root \
  --mode real --t0-source hit --t0-scale-to-ns 1e-3 \
  --t0-offset-ns 0 --t0-policy median --cut-preset legacy
```

Use `./CeeV0Analysis.bin --help` for all individual cut overrides.

## ROOT macros

The run wrapper calls the compiled executable:

```bash
root -l -b -q 'macros/runCeeV0Analysis.C("/data/sim.root","lambda.root")'
```

The wrapper rejects T0 hit/scale arguments in simulation. It deliberately does
not pass a T0 scale by default, so a real-data `CeeT0Hit` run still stops until
the scale is explicitly confirmed,
for example the sixth macro argument `1.0` when the persisted value is already
in ns.
The ninth macro argument enables the incomplete-Jacobian QA stream and defaults
to `false`, preserving strict output.

The QA macro overlays all/truth/background mass, daughter DCA, decay length,
pointing cosine, geometric chi2, total topological chi2, the PV-only constraint,
and raw daughter-PV chi2. Truth overlays are omitted for real data or unresolved
truth:

```bash
root -l -b -q 'macros/plotCeeV0.C("lambda.root","lambda_QA.pdf")'
```

The matching/T0 helper is ROOT-independent. `make check` runs it and also
parses the complete analyzer against compile-shape API stubs. This catches C++
regressions on a machine without ROOT, but only `make` in the real CEE
environment verifies actual dictionaries, schema evolution, linking, and ROOT
I/O.
