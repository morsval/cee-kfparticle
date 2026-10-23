# v0.7.0：BenoDst 输入、通道与用户手册

在仓库根目录或 CEErco 执行 `make -j4`，一次编译核心静态库、分析 CLI 和 Maker 共享库。先加载 ROOT/CeeRoot/FairRoot 环境。`make check` 运行独立回归；旧 Analyser 用 `make legacy`。

- 原两体：`Make/Run`；已有 TOF 与 tracking 算法保持。
- 自定义 tracking 通道：`Channels().channels/pid` + `MakeChannels/RunChannels`。
- 两体、共同 SV 三体、复合 daughter 级联；通道依赖排序和共享末级 track 排除。
- 动态 daughter 角色、输入/拟合状态、DCA、candidate/track 链接及 ROOT 输出。
- 示例：`CeeV0Analysis/macros/runCeeChannels.C`。默认 PID 不伪造标定；最后 true 是显式 charge-only 流程演示。
- BenoDst：`CeeKFParticleBenoSource` 直接读取 ROOT/list；零场复用保存的 track–TOF 匹配及已扣 T0 的时间。5D 的前四维几何协方差保留相关项；tracking 默认支持已核对的单位电荷参数化。无需用户重写 reader。
- macro 独立选择 `"raw"/"benodst"` 与 `"tof"/"tracking"`；旧的数值 By 参数调用兼容。
- 检测到 CeeRoot Beno EDM 头文件时，顶层 `make` 自动加入 Beno 支持；可显式设置 `WITH_BENODST=1/0`。
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

## 当前零场 BenoDst 的运行入口

```bash
cd CEErco/CeeV0Analysis
root -l -b -q 'macros/runCeeKFParticle.C("input.list","mothers.root","benodst","tof",0.,-1,true)'
```

`true` 指实验数据。修改 macro 中的 T0/TOF 分辨率示例为实际标定值；Beno 的时间 scale=1、offset=0，不能再次减 T0。Beno 缺少 NDF/fitConverged，macro 不应用这两种基础 cut。`requireUniqueMatch`、共享 hit 排除和双 TOF 优先级在 `BenoInputOptions` 中设置。

新输入经过真实上游 Beno EDM + ROOT 合成事件回读测试，尚未使用真实实验文件闭合。核素的电荷/刚度转换用 `SetStateConverter` 明确提供。三体和级联仍只支持 tracking。
