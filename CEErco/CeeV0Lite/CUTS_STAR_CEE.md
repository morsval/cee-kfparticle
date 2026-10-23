# CeeV0Lite、STAR KFParticle 与现有 CEE cuts 对照

v0.6.2 新增的 `motherTopo*` 是独立联合 PV 约束的输出，目前不参与本表的 cut。
本表原 `topologicalChi2Ndf` 仍指七维母粒子的残差统计量。
新旧字段请分开比较，不能直接交换阈值；详见 `MIGRATION_0_6_2.md`。

本表已经按仓库中的 `KFParticleFinder.h/.cxx`、
`StKFParticleInterface.cxx` 和 `StKFParticleAnalysisMaker.h/.cxx` 重新核对。
“对应”表示物理意义接近，不表示两个实现的数值定义或候选保存阶段完全相同。

| CeeV0Lite 配置 | STAR Finder 对应 | STAR 构造/有效值 | STAR AnalysisMaker 名称与实际状态 |
|---|---|---:|---|
| `maximumDaughterDcaCm` | `fDistanceCut` | `1 cm` | 仅近似对应：轻量版是两轨迹分别到拟合顶点最近点的间距，STAR 是 Finder 拟合前的 daughter--daughter 距离；两者都不是 `DcaProton/DcaPion` 到 PV 的 DCA |
| `maximumGeometricChi2Ndf` | `fCuts2D[1]` | `10` | `ld_chi2ndf` / `mLdChi2ndf`；`PassLambdaSelection` 中 `<=5`，但调用被注释 |
| `minimumDaughterChi2FromPv`（raw chi2） | STAR interface 的 `fChiPrimaryCut` / `GetDeviationFromVertex(PV)` | 活动 interface 默认值 `18.6` | `chi2primary_proton/pion` 也是 raw chi2；`PassLambdaSelection` 中可选 `>=5`，但该调用被注释。interface 仍用 `18.6` 将输入 track 分到 primary/secondary 列表 |
| `minimumDaughterDeviationFromPv`（deprecated，`sqrt(chi2)`） | 无直接 STAR 数值对应，仅为 CeeV0Lite 0.5.x 兼容 | `star` preset 中为 `0` | 旧 `--min-daughter-pv-deviation` 仍控制它；不能把 raw `18.6` 填到这里 |
| `minimumDecayLengthSignificance` | `fCuts2D[2]`，后续还有 `fSecCuts[2]` | `5`；普通构造的 `fSecCuts[2]=-10` | `ld_ldl` / `mLdLdl`；可选 `>=2`，但调用被注释 |
| `minimumDecayLengthCm` | 近似 `fLCut` | `5 cm` | `ld_l` / `mLdL`，无活动 Maker cut；STAR 使用带方向判断的 line distance，轻量版使用 `|SV-PV|` 并另算 pointing cosine |
| `maximumTopologicalChi2Ndf` | `SetProductionVertex(PV)` 后的 mother total `Chi2/NDF`，数值上对应 `fSecCuts[1]`/candidate selection | `5` | `KFParticleTopoReconstructor::SelectParticleCandidates()` 中存在活动的 `<5` candidate-competition 选择；Maker 的 `ld_chi2topo` / `mLdChiTopo` 也由 production-vertex 约束后的 total `Chi2/NDF` 得到，`PassLambdaSelection` 中重复的 `<=5` 调用被注释 |
| `minimumMassGeV/maximumMassGeV` | 仅近似 `fSecCuts[0]` | `3 sigma` mass pull | 原始 Lambda 仍可先进入 `Particles`；Maker 的活动硬 cut 是 `1.08 <= ld_mass <= 1.18 GeV/c^2` |
| `minimumPointingCosine` | 无一一对应公开 setter | -- | 与 PV 拓扑约束有关，但数值不能等同于 `chi2Topo`；头文件中的 `V0DirectionCut` 当前没有活动调用 |
| `daughterDcaToPvCm` | 无标准 V0 一一对应 cut | -- | `DcaProton/DcaPion` / `mProtonDCA/mPionDCA` |
| `lambdaRoleDcaBoundaryCm` | 无 STAR 对应 | -- | CEE 零场特有的 proton/pion 角色分配 |

## STAR 接口中可见的设置

当前头文件公开：

```cpp
SetMaxDistanceBetweenParticlesCut(cut);
SetLCut(cut);
SetChiPrimaryCut2D(cut);
SetChi2Cut2D(cut);
SetLdLCut2D(cut);
SetSecondaryCuts(sigmaMass, chi2Topo, ldl);
```

`Set2DCuts(chi2Prim=3, chi2=3, ldl=5)` 实际只写入 `fCuts2D[0]`，后两行
被注释。`SetSecondaryCuts` 的头文件默认参数是 `(3,5,10)`；这与普通构造函数的
`fSecCuts=(3,5,-10)` 不同。因此，只有显式无参数调用 `SetSecondaryCuts()` 才会把
第三项从 `-10` 改成 `+10`。当前 `StKFParticleAnalysisMaker` 没有这样调用。
这里的 Finder `fCuts2D[0]` 与 `StKFParticleInterface::fChiPrimaryCut` 是两个不同变量；
本版 `star` preset 的 raw daughter cut `18.6` 对应的是后者在输入 track 分类中的活动
默认值，不是把 `SetChiPrimaryCut2D(18.6)` 写进 Finder。

Maker 中虽然定义了完整的 `PassLambdaSelection()`，但调用整行被注释。Maker 层当前
真正活动的 Lambda 后选择只有 `1.08 <= ld_mass <= 1.18 GeV/c^2`；Finder 上游的
距离、几何质量、方向和衰变长度选择仍然会在候选到达 Maker 前执行。此外，
`StKFParticleInterface` 的活动 `fChiPrimaryCut=18.6` 在输入分类阶段使用 raw
`GetDeviationFromVertex()` chi2，`KFParticleTopoReconstructor::SelectParticleCandidates()`
也会在 candidate competition 中使用 production-vertex 约束后的 total `Chi2/NDF<5`。

0.6.0 中轻量版明确保存相同的“几何贡献再加 PV 约束贡献”数据流：

```text
topologicalChi2    = geometricChi2 + primaryVertexConstraintChi2
topologicalNdf     = geometricNdf + 2
topologicalChi2Ndf = topologicalChi2 / topologicalNdf
```

PV-only 的 `primaryVertexConstraintChi2{,Ndf}` 另存为 QA，不能当作 total topology。
daughter 的 canonical 配置是 raw chi2 `minimumDaughterChi2FromPv`；`star` preset
将其设为 `18.6`。旧 `minimumDaughterDeviationFromPv` 和旧 CLI 仍是
`sqrt(chi2)` 兼容量，默认不参与 `star` preset。

## CeeV0Lite 的其他可调选择

这些不是 STAR V0 拓扑 cut 的直接对应项：

- `momentumMode`：TOF、tracking 或逐 track 自动选择；
- `magneticFieldByTesla`：均匀 `By`；
- `maximumTofIterations`、`tofVertexConvergenceCm`、
  `tofTimeConvergenceNs`、`requireTofConvergence`：TOF—顶点外层迭代；
- `maximumBeta`：拒绝非物理或过于接近 1 的 TOF 解；
- `additionalDecayTimeSigmaNs`：衰变时间模型的附加误差；
- `maximumCandidates`：事件级保护上限；
- `reconstructK0Short/reconstructLambda/reconstructAntiLambda`：通道开关；
- `requireMeasuredChargeInMagneticField` 和
  `allowTofMomentumInMagneticField`：有场数据的输入完整性保护。

当前 `star` preset 只映射本库有明确定义且已在本轮核对的 cuts。它没有新增 STAR
可能使用的额外 pointing 或最大衰变长度选择。有限场 transport 的 point-POCA 利用
横向周期性把全局搜索缩到纵向最优点附近至多两个周期，枚举该窗口内全部驻点、排除
局部极大值并与原区间端点比较后选全局最小，但这仍不是两螺旋同步的
全局 DCA，也不包含非均匀 field map 或材料模型。

## 与现有 AnaMaker 的 DCA 角色分配

现有 CEE 代码的三维直线 DCA 是

```text
DCA = |(r_track - r_PV) cross unit_direction|.
```

选择规则是：

| track 来源 | DCA 边界 | proton 角色 | pion 角色 |
|---|---:|---:|---:|
| TPC | 1.2 cm | `DCA <= 1.2` | `DCA > 1.2` |
| MWDC | 5.0 cm | `DCA <= 5.0` | `DCA > 5.0` |

CeeV0Lite 用每条输入 track 的 `lambdaRoleDcaBoundaryCm` 表示这个来源相关边界，
并将计算结果保存到 `Candidate::daughterDcaToPvCm[2]`。这里 daughter 0 是
(anti)proton 假设，daughter 1 是 pion 假设。

该规则会偏向 proton 更接近 Lambda 飞行方向、pion 张角更大的事件。它适合复现当前
CEE 分析选择，但不是普适 PID；应分别研究正确分配效率、交换分配率和相空间偏差。
