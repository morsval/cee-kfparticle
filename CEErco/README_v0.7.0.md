# v0.7.0 通道与用户手册更新

在仓库根目录或 CEErco 执行 `make -j4`，一次编译核心静态库、分析 CLI 和 Maker 共享库。先加载 ROOT/CeeRoot/FairRoot 环境。`make check` 运行独立回归；旧 Analyser 用 `make legacy`。

- 原两体：`Make/Run`；已有 TOF 与 tracking 算法保持。
- 自定义 tracking 通道：`Channels().channels/pid` + `MakeChannels/RunChannels`。
- 两体、共同 SV 三体、复合 daughter 级联；通道依赖排序和共享末级 track 排除。
- 动态 daughter 角色、输入/拟合状态、DCA、candidate/track 链接及 ROOT 输出。
- 示例：`CeeV0Analysis/macros/runCeeChannels.C`。默认 PID 不伪造标定；最后 true 是显式 charge-only 流程演示。
- BenoDst：`examples/BenoDstAdapterExample.h` 要求经过验证的 5D→6D 转换 callback；不是完整即用 reader。具体字段核对见手册。
- 手册生成入口：`docs/manual_source/build_user_book_v070.py`，13 章、动态目录、连续排版。

**范围：**通用 Channels 只支持 tracking 动量和均匀 By；零场 TOF 三体/级联尚未实现。中间态使用原始 7D，不使用碰撞 PV 约束副本；没有自动逐级 8D 飞行链。真实 CEE/FairRoot 数据闭合仍需实验环境。

以下为原两体 AnalysisMaker 的兼容说明：

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
