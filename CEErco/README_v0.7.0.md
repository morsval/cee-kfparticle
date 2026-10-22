# CEE KFParticle 0.7.0

- Start: [AnalysisMaker 使用与迁移](docs/AnalysisMaker_v0.7.0.md)
- Macro: [runCeeKFParticle.C](CeeV0Analysis/macros/runCeeKFParticle.C)
- Variables: [中文对照表](docs/CandidateVariables_v0.7.0.md) / [CSV](docs/CandidateVariables_v0.7.0.csv)
- Full manual source: [docs/manual_source](docs/manual_source)

Build in a CEE environment: `make -C CeeV0Analysis shared`.
Standalone regression: `make -C CeeV0Lite test` and `make -C CeeV0Analysis check`.
ROOT-only integration: `make -C CeeV0Analysis root-writer-test root-source-test`.
The latter uses synthetic CEE-shaped fixture dictionaries and does not link real FairRoot/CEE EDM.

The v0.6.5 magnetic track–PV DCA reference convention and all fit equations are preserved.
