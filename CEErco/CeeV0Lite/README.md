# CeeV0Lite

**v0.6.5（2026-10-22）**：有磁模式默认将每条输入 track 在内存中转换到与当前重建 PV 的 DCA 参考点，完整输运位置、动量与 6×6 协方差；配对前只做一次。详见根目录 `MIGRATION_0_6_5.md` 与 `v0.6.5_review`。切换 `--magnetic-track-reference stored` 可重放旧参考约定。零场 TOF 不变。更新后需要清理并重编核心与分析程序。

**v0.6.4（2026-10-14）**：保守优化固定路径解析导数、协方差矩阵乘法和 POCA 临时数组分配。共同 SV、POCA 数值导数、TOF 迭代和物理 cuts 保留。详情见 `MIGRATION_0_6_4.md`；验证与补丁见 `v0.6.4_review`。下方 v0.6.3 说明保留为历史功能背景。

**v0.6.3（2026-09-30）**：正式 C++ 核心新增均匀 By 下带电母粒子的八维 PV 约束，
同时处理位置、动量旋转及 S 的完整相关误差。保留原中性分支与 V0 选择。
ROOT 输出格式为 6，新增 motherTopoCharge 和 motherTopoByTesla。
超氚集成示例位于 CeeV0Lite/examples/hypertriton_by；默认分析入口仍为原 Λ/V0 通道。
接口、边界及重新编译步骤见 `MIGRATION_0_6_3.md`。



> 2026-09-20 磁场方向修正：CEE 束流沿 z，均匀磁场沿 y；使用带符号 `--by-tesla` / `magneticFieldByTesla`。详见 `MIGRATION_BY_20260920.md`。
`CeeV0Lite` 是面向 CEE 的轻量级二级顶点重建代码。目前支持

- `Lambda -> p + pi-`
- `anti-Lambda -> anti-p + pi+`（默认关闭）
- `K0S -> pi+ + pi-`（代码保留、测试覆盖，但默认关闭）

它不依赖 STAR 的 `KFParticleBase.h/.cxx`、`KFParticleDef.h`、SIMD、
`KFParticleFinder` 或其他 KFParticle 类，也不依赖 ROOT。源码使用 C++11，
Finder 是串行的 `O(N^2)` 二体组合。

## 设计边界

主要模块如下（另有 CeeV0LiteMath 提供矩阵运算）：

| 模块 | 作用 |
|---|---|
| `CeeV0LiteTypes` | 输入 track、主顶点、六维 track 状态、七维 mother state 和 V0 candidate |
| `CeeV0LiteTransport` | 在均匀 `By` 中输运六维状态和 6x6 协方差；`By=0` 自动成为直线 |
| `CeeV0LiteVertexFitter` | 两条 daughter 的公共三维顶点拟合；联合 Kalman 更新后构造 `(SV,P,E)` mother state |
| `CeeV0LiteProductionVertex` | 中性／带电母粒子 PV 约束、S、八维协方差及端点输运 |
| `CeeV0LiteFinder` | TOF 动量迭代、共享 T0 传播、`Lambda/K0S` 组合、质量/拓扑量和 cuts |
| `CeeV0LiteCeeAdapter` | 将 `CeeTpcTrack/CeeTpcVertex` 的 getter 转成轻量输入 |

这里借鉴的是 KFParticle 的核心思路——状态、协方差、输运、逐次 Kalman
更新和拓扑约束——而不是照搬 STAR 的类层次。

## 本次零磁场数据

将 `FinderConfig::momentumMode` 设为 `kTofMomentum`，并令
`magneticFieldByTesla = 0`。每条 track 需要：

- TPC 上一点 `(x,y,z)` 和直线方向；
- 匹配后的 iTOF 或 eTOF hit 位置；
- TOF hit 时间；
- 时间原点：模拟取生成器的精确 `t=0`，真实数据才来自 T0；
- 位置、方向、TOF 时间和路径长度的不确定度。

对质量假设 `m`，daughter 动量由

```text
Delta t = t_TOF - t_decay
L       = (r_TOF - r_decay) dot unit_direction
beta    = L / (c * Delta t)
p       = m * beta / sqrt(1 - beta^2)
```

得到。`r_decay/t_decay` 起初未知，因此 Finder 会做

```text
直线交会 seed -> TOF daughter 动量 -> Kalman 顶点
-> parent 从主顶点到衰变点的飞行时间 -> 重算 TOF 动量
```

直到顶点和时间收敛。生产默认会拒绝达到最大迭代次数仍未收敛的候选；近平行/反平行
track 和最近点输运未收敛也会被拒绝。`CeeiTOFHit` 和 `CeeEtofHit`
只要最终都提供 `(x,y,z,time)`，在本库中不需要分成两个类；探测器几何投影和事件级
唯一匹配由 `CeeV0Analysis` 处理。真实数据中 `CeeT0Hit` 的时间用于
`EventInput::startTimeNs`，不是 daughter TOF hit。

模拟数据不需要 `CeeT0Hit`。这里仍写 `startTimeNs=0`，只是声明 MC TOF 时间所采用的
坐标原点，并非加入了一条 T0 测量。模拟必须同时使用
`sigmaStartTimeNs=0`、`additionalDecayTimeSigmaNs=0` 和
`propagateSharedTimeUncertainty=false`。每个 TOF hit 自身的 `sigmaTofTimeNs` 仍应按
digitization 分辨率填写，因为它和 T0 不是同一个误差来源。

每条 track 的边缘协方差只包含该 track 自己的 TOF 时间误差。在真实数据模式，事件 T0 误差和
`additionalDecayTimeSigmaNs` 被视为同一 candidate 两条 daughter 共享的时间 nuisance：
Finder 用 `T0 +/- sharedTimeDerivativeStepNs` 完整重跑上述 fixed-point，得到七维 mother
state 对共同时间的导数，再只加入一次 rank-one covariance。这样不会把同一个 T0
当成两次独立测量误差。

零磁场直线本身不能给出电荷正负。若当前数据也没有别的电荷信息，设置
`chargeKnown=false`；Finder 会测试允许的质量/电荷假设，因此组合背景会比有场时高。
应尽量用 TOF、dE/dx 或已有 PID 信息收紧 `hypothesisMask`。

当前默认按现有 CEE `AnaMaker` 的经验规则给 Lambda daughter 分角色：相对主顶点
三维 DCA 不大于 `TrackInput::lambdaRoleDcaBoundaryCm` 的 track 可作 proton，DCA
大于该边界的 track 可作 pion。TPC 默认边界为 1.2 cm；若接入 MWDC track，应给该
track 设置 5 cm。这个 cut 是运动学角色分配，不是基本粒子鉴别，会改变接受度，必须用
模拟评估效率和错配率。可用 `FinderConfig::useLambdaDcaRoleAssignment=false` 关闭。

## 未来有磁场数据

若 TPC tracking 已经从曲率给出动量与电荷，使用：

```cpp
config.momentumMode = CeeV0Lite::kTrackingMomentum;
config.magneticFieldByTesla = measuredByTesla;
```

输运器使用均匀 `By` 的螺旋解析轨迹，并用数值 Jacobian 做
`C_out = F C_in F^T + Q`；最近点输运的每个扰动都会重新求 path，因此包含
`ds_closest/dstate`。有限 `By` 下的“轨迹到给定点”POCA 利用横向螺旋的严格周期性，
先把全局搜索严格缩到纵向最优点附近至多两个周期，再枚举该窗口内全部驻点、排除局部
极大值，并与原有界区间端点一起选距离的全局最小值；这避免了
Newton 迭代落到错误螺旋圈或极大值。它仍不是“两条螺旋同时变化 path 的全局 mutual
DCA”。接口无需改变。当前版本只实现均匀 `By`；将来若使用
三维非均匀场图，只需新增 transport backend，Finder 和顶点 fitter 不必重写。
有场情况下默认要求 tracking 给出可靠电荷；默认不会混用“TOF 给模长、曲率给轨迹”的
半状态。

## 最小调用方式

```cpp
#include "CeeV0LiteFinder.h"

CeeV0Lite::EventInput event;
event.primaryVertex.positionCm = CeeV0Lite::Vec3(pvx, pvy, pvz);
event.startTimeNs = 0.;       // simulation generator time origin
event.sigmaStartTimeNs = 0.; // no CeeT0Hit in simulation

std::vector<CeeV0Lite::TrackInput> tracks;
// 每个重建 track 只加入一次；先在 AnaMaker 中完成 iTOF/eTOF 匹配和 PID。
// 填 positionCm、direction、tofPositionCm、tofTimeNs、误差和 hypothesisMask。

CeeV0Lite::FinderConfig config;
config.momentumMode = CeeV0Lite::kTofMomentum;
config.magneticFieldByTesla = 0.;
config.reconstructK0Short = false;
config.reconstructAntiLambda = false;
config.useLambdaDcaRoleAssignment = true;
config.additionalDecayTimeSigmaNs = 0.;
config.propagateSharedTimeUncertainty = false;

CeeV0Lite::Finder finder(config);
std::vector<CeeV0Lite::Candidate> candidates;
finder.Reconstruct(tracks, event, candidates);
```

`CeeV0LiteCeeAdapter.h` 会优先读取 `CeeTpcTrack::HasVertexSeed()` 对应的同一参考面
状态。有场 tracking 模式使用完整 6x6 covariance，并严格要求 seed、covariance 和
拟合电荷均有效。零场 TOF 模式保留 seed 的 position block，并把 momentum block
压缩成 `sigmaDirectionRad` 后用 TOF 重建动量模；当前尚未把 seed 中各向异性的方向误差
及 position--direction cross block 完整变换到新的 TOF state。主顶点 adapter
选择最高 rank 的有效顶点并使用 `GetXErr/YErr/ZErr()`，不再直接访问未检查的 `[0]`。
详细接口见 `docs/CEE_ADAPTER.md`。
这些字段必须由 tracking/vertex producer 真正写入 ROOT 文件；当前数据的兼容回退和
未来有场数据的 producer 要求见 `docs/CEE_INPUT_PRODUCER.md`。

完整事件读取、模拟/真实模式、T0、事件级唯一 TOF 匹配及 ROOT 输出由同级目录
`../CeeV0Analysis` 提供。不要用 MC track id 做重建匹配；MC 信息只能用于模拟样本的
效率和纯度标签。

## 配置和 cuts

默认质量窗为 `0.35--0.65 GeV/c^2`（`K0S`）和
`1.05--1.25 GeV/c^2`（`Lambda`）。`ChannelCuts` 还包含 daughter DCA、
几何 chi2/ndf、拓扑 chi2/ndf、衰变长度、长度显著性、pointing cosine 和
daughter 相对主顶点的 raw chi2。其定义为

```text
topologicalChi2Ndf =
  (geometricChi2 + primaryVertexConstraintChi2) / (geometricNdf + 2).
```

`primaryVertexConstraintChi2{,Ndf}` 另行保存，便于观察 mother 与 PV 相容性的单独贡献；
它不是 `topologicalChi2Ndf` 的别名。canonical daughter cut 是 raw chi2
`minimumDaughterChi2FromPv`；`star` preset 使用 STAR interface 的活动值 `18.6`。
旧 `minimumDaughterDeviationFromPv` 仍按 `sqrt(chi2)` 解释，仅用于兼容旧 steering。
因此不要把这两个 cut 的数值直接互换。

这些默认值只用于让流程工作，不是 CEE 最终物理 cuts。正式使用前必须用真实重建/模拟
标定：

- `positionCovarianceCm2` 和 `sigmaDirectionRad`；
- iTOF/eTOF/T0 时间分辨和共同时间偏移；
- `sigmaPathLengthCm`；
- tracking 动量的完整 6x6 协方差；
- 各拓扑量和不变质量窗。

如果输入协方差不真实，代码仍可能运行，但 chi2、显著性和质量误差没有可信的统计含义。

从 0.5.0 起，两条 track 在同一个冻结线性化点上做对称 transverse vertex batch fit；
消去两个 path nuisance 后只使用一个有效的二轨 DCA residual。对已通过廉价中心值 cuts
的 candidate，程序把当前实现的 POCA/vertex estimator 当作 `12D -> 7D` 映射，逐个
扰动输入 state 并重新求 path/vertex，得到 mother covariance，顺序固定为
`(x,y,z,px,py,pz,E)`。`sigmaMass` 直接由其中的 `(P,E)` 4x4 block 传播。这与 STAR
“先构造 mother、再由 mother covariance 求质量误差”的数据流一致，不再用两块独立
daughter momentum covariance 简单相加。CEE 特有的共享 T0 也以上述共同 nuisance
方式进入同一 mother covariance。

为避免把数值扰动花在明显本底上，上述 mother covariance 只在 candidate 已通过质量、
daughter DCA、几何量、pointing 和角色分配等廉价中心值检查后计算；最终
`L/sigmaL` 与上述 total topological chi2 仍使用加入共享 T0 后的 covariance 再做 cut。

默认情况下，12 个输入列的正负扰动必须全部拟合成功，即
`jacobianFitsAttempted=24`、`jacobianFitsSucceeded=24`。如果中心 POCA、顶点、TOF
fixed-point 和 mother 中心值均成功，但某一个 `+h/-h` 扰动拟合失败，状态记为
`kMotherCovarianceJacobianPerturbationFailed`。所有 preset 默认仍拒绝它。只有显式设置
`ChannelCuts::requireCompleteMotherJacobian=false`（分析命令行
`--allow-incomplete-mother-jacobian`）才保留这类 QA candidate。保留时中心质量、SV、
动量、DCA、geometric chi2、decay length 和 pointing 可用；最终 covariance、质量/顶点/
长度误差、`L/sigmaL`、PV constraint 与 topological chi2 全部为 `NaN`，并且
`covarianceDependentCutsEvaluated=false`。中心解失败、最终矩阵非 PSD、共享时间导数失败和
不支持的 process noise 仍然无条件拒绝。

但这仍不是 STAR 的逐行复刻：库没有长期保存 `PV + SV + daughters + T0` 的完整联合
状态；上游输入中不存在的 track--track、PV--track correlation 也无法恢复；TOF 模式
目前只有共同 T0 会完整重跑 fixed-point，PV error 和其他独立误差尚未通过整个
fixed-point 反馈到 mother momentum；均匀场、材料模型与 STAR 也不同。非零 transport
process-noise 配置会被 mother covariance 明确拒绝，直到这些独立 noise modes 被加入
外层传播。因此在“输入两条 track 的上游误差相互独立”这一前提下，`sigmaMass` 已纳入
当前 vertex estimator 诱导的有效 daughter cross term，并解决了 double-T0 问题，但其
统计覆盖仍必须用 CEE simulation pull/closure 验证后再定正式 cut，不能直接复制 STAR
数值。

0.6.0 改变了 `topologicalChi2Ndf` 的语义和 ROOT schema；从 0.5.x 升级时请先读
`MIGRATION_0_6.md`。
0.6.1 新增了可选的 Jacobian 失败 QA 行与 output-format 4；启用前请读
`MIGRATION_0_6_1.md`。

## 编译和接入 `Makefile(1)`

独立验证：

```bash
make
make test
```

建议把整个目录命名为 `CeeV0Lite/`，放在 `analysis.cxx` 和实际 Makefile 同级。
把 `integration/Makefile.fragment` 中的内容合并进 `Makefile(1)`。它只添加 include
路径和五个 object；你原来的 ROOT/FairRoot/CEE 链接参数保持不变，也不需要再链接任何
KFParticle 库。

## 已覆盖测试

- `By=0` 的正反向直线状态输运；
- `By=1 T` 的螺旋位置/动量和反向输运；
- 两条 track 的公共顶点 Kalman 拟合；
- 最近点路径边界、迭代不收敛和近平行退化保护；
- 有限 `By` point-POCA 的多驻点、错误螺旋圈/局部极大值反例和有界全局最小选择；
- tracking 动量模式下的 `Lambda` 和 `K0S`；
- 零场 TOF 迭代模式下的 `Lambda` 和 `K0S`。
- `By=0.5 T` 的有场 Lambda Finder 端到端测试；
- mock `CeeTpcTrack/CeeTpcVertex` adapter 测试。
- 七维 mother covariance 的 PSD、非相交 daughter 交换不变性、沿轨迹 reference-point
  平移不变性、质量误差传播和共同 nuisance rank-one 测试；
- 共享 T0 的增量逐元素等于一次 `dY/dT0 * sigma(T0)^2 * dY/dT0^T`，并且
  `sigma(T0)` 加倍时增量变为四倍的端到端测试。
- total/PV-only topology 的 numerator--NDF 恒等式，以及 raw daughter chi2 与
  `sqrt(chi2)` 兼容 cut 的边界测试。
- 中心 candidate 成功但固定的 Jacobian 正/负扰动失败时，24 次诊断完整记录、默认严格拒绝、
  显式 QA 保留以及全部 covariance-dependent 输出为 `NaN` 的回归测试。

测试是解析合成事例，用于防止程序和数学实现回归；真实实验分辨、效率和本底仍需 CEE
数据/模拟验证。

有场合成测试证明的是当前均匀 `By`、有限路径情形的代码闭合，不代表低动量或长路径
生产数据已经完成验证：初始顶点 seed 仍由两条局部切线给出，之后才进行螺旋输运迭代。
point-to-helix POCA 已在有界区间内选全局最小，但 `daughterDcaCm` 仍不是两条螺旋的
同步全局 DCA；非均匀 field map、材料效应和正式场图 closure 也不在本轮实现范围。
这些相空间必须在正式有场数据前单独做 simulation closure。

与 STAR KFParticle 及现有 CEE cuts 的逐项对应见 `CUTS_STAR_CEE.md`。
本版 mother covariance、质量误差和共享 T0 的公式、ROOT packing 与仍存在的边界见
`docs/STAR_LIKE_MOTHER_COVARIANCE.md`。
