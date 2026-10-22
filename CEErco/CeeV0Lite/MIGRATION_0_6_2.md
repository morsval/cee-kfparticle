# v0.6.2：中性母粒子的 S、PV 约束与八维状态

基线：2026-09-20 修正磁场方向后的 v0.6.1。CEE 束流沿 z、磁场沿 y 的约定保留。

## 这次增加了什么

原来的 `Candidate::mother` 仍是衰变点处的七维母粒子。新增加
`Candidate::motherTopo`，保存独立的 PV 约束副本：

`(x,y,z,px,py,pz,E,S)`，配套完整的 8×8 协方差。

`S = p·(SV-PV_fit)/|p|²`，其单位为 cm/(GeV/c)。这是有符号的飞行参数，
不是时间。有效状态默认位于拟合后的产生端；`S` 连接同一对象的产生端和衰变端。
`FlightObservables` 提供 `L=S|p|`、`ctau=S*m` 以及对应误差。
`ctau` 的单位是 cm；除以 29.9792458 cm/ns 才得到固有时间 ns。
质量使用约束副本当前的重建质量，没有施加 PDG 质量或寿命约束。

加入 S 的价值是把端点关系与相关误差放进同一个可继续输运的状态。
仅为了画一个 cτ 中心值分布，无须这套八维拟合。

## 实际计算过程

1. 先执行原来的候选构造、七维协方差计算和全部选择。
2. 对通过原选择的候选，建立 `(SV,p,E,PV)` 的十维联合先验。
3. 对“两端连线与母粒子动量平行”施加两个独立约束。
4. 迭代重新线性化，但每轮都使用同一份先验，避免重复计入 PV。
5. 从联合后验计算 S 及八维协方差，再输运到产生端保存。
6. 计算 `topoFlight`；无效结果保留状态码，数值为 NaN。

这里的 PV 是拟合中的带误差变量。得到约束副本时，母粒子的 SV、动量、能量及
联合拟合中的 PV 估计都可能调整；输入 `EventInput::primaryVertex` 不被覆盖。
返回的衰变端属于新副本，也不应混同于原七维母粒子的 SV。

数学上采用 CEE 的中性粒子联合约束拟合，借鉴 STAR 的状态语义；
不是逐行移植 STAR `SetProductionVertex()`，不承诺逐候选数值一致。

## PV 与 TOF 动量的相关性

零场 TOF 求解已经使用 PV 来估计母粒子的飞行时间。新实现对 PV 的三个坐标
做正负扰动，重新执行完整 `SolvePair()`；同时比较 h 与 h/2 的导数，失败或
不稳定时不输出有效八维状态。默认最多增加 12 次 TOF 求解。

记 `D = d(mother)/d(PV)`、`V = Cov(PV)`、原七维误差为 `C_cond`，则新拟合使用：

```
Cov(mother)_marginal = C_cond + D V D^T
Cov(mother, PV)      = D V
Cov(PV)             = V
```

因此不能只给母粒子增加 PV 方差，却漏掉交叉块。共同平移应当在两端差值中抵消。
这个扩展仅用于新八维副本；原七维输出及其选择含义保持不变。

跟踪动量模式中，当前输入模型假定原始轨迹测量与 PV 独立。TOF 模式补入的是
算法本身引入的 PV 响应；上游若用同一批子轨迹拟合 PV，还会有额外相关性。
Finder 没有这些相关性的输入，不能宣称已经消除它们。独立接口
`ProductionVertexInput::motherVertexCovariance` 支持调用者提供已知交叉块；
此时 `mother.covariance` 必须是已经包含 PV 影响的边缘协方差。

原七维误差仍沿用现有轨迹误差近似，共享 T0 按原开关处理。
这次并没有重做所有原始 TOF hit、路程、方向误差的完整联合统计模型。

## 新接口与输出

- `CeeV0LiteProductionVertex.h/.cxx`：联合 PV 拟合、八维状态和端点输运。
- `FitProductionVertex()`：接受联合输入，不修改原对象；只支持中性母粒子。
- `TransportToDecayVertex()` / `TransportToProductionVertex()`：中性端点输运，
  同时更新协方差。S 保持不变；重复请求当前端点不会重复移动。
- `GetFlightObservables()`：计算质量、L、cτ 及一阶误差，包含 S 与 p/E 的交叉项。
- `FinderConfig::buildProductionState` 默认 true。命令行可用
  `--no-production-state` 关闭新计算，以比较结果和运行时间。

ROOT 输出格式从 4 升为 **5**。原有分支保留，新增：

| 分支 | 含义 |
| --- | --- |
| `motherTopoValid`, `motherTopoStatus` | 八维拟合有效性与失败原因 |
| `motherTopoCorrelationModel` | 0：外部联合输入；1：tracking 独立 PV；2：TOF PV 响应 |
| `motherTopoAtProductionVertex`, `motherTopoIterations` | 当前端点与约束迭代次数 |
| `motherTopoParameters[8]` | x,y,z,px,py,pz,E,S；默认位于产生端 |
| `motherTopoCovariance[36]` | 对应八维协方差，下三角逐行存储 |
| `motherTopoS`, `motherTopoSigmaS` | S 与误差；S 同时也是 Parameters[7] |
| `motherTopoPvChi2`, `motherTopoPvNdf` | 新联合约束的 χ²；有效时 NDF=2 |
| `motherTopoChi2`, `motherTopoNdf`, `motherTopoChi2Ndf` | 几何 χ² 加新 PV χ² 及总自由度 |
| `motherTopoFlightValid` | 下列派生量是否有效 |
| `motherTopoMassGeV`, `motherTopoSigmaMassGeV` | 约束副本的质量与误差 |
| `motherTopoDecayLengthCm`, `motherTopoSigmaDecayLengthCm` | 有符号 L=S|p| 与误差 |
| `motherTopoDecayLengthSignificance` | 上述 L/σL；σL=0 时为 NaN |
| `motherTopoCtauCm`, `motherTopoSigmaCtauCm` | cτ 与误差，单位 cm |

原 `massGeV` 仍来自未约束副本。原 `decayLengthCm` 是欧氏距离，不能直接
当作新字段的别名。原 `decayTimeNs` 是实验室系衰变时刻。
原 `--max-topo-chi2ndf` 等选择仍使用旧字段；新 `motherTopo*` 不自动参与选择。
若后续要使用新量做 cut，需要单独研究效率与偏差。

## 失败状态

| 值 | `ProductionFitStatus` | 处理建议 |
| --- | --- | --- |
| 0 | NotRun | 尚未调用 |
| 1 | Valid | 可使用八维状态 |
| 2 | Disabled | 配置关闭 |
| 3 | InvalidMother | 原母粒子或其协方差无效 |
| 4 | PvResponseFailed | TOF 未收敛，或 PV 导数失败／不稳定 |
| 5 | InvalidJointCovariance | 联合协方差非有限、非对称或非半正定 |
| 6 | SingularConstraint | 约束方差不可逆，或方向落出本地参数化范围 |
| 7 | NotConverged | 迭代没有满足终止条件 |
| 8 | UnphysicalState | 拟合后质量／误差／输运无效 |
| 9 | UnsupportedCharge | 输入为带电母粒子，本版明确拒绝 |
| 10 | InvalidConfiguration | 配置或几何拟合统计量非法 |

新增计算失败不会删掉已经通过原选择的候选，也不会拿 S=0 伪装成有效寿命。
应先检查 `motherTopoValid`；做寿命图时再检查 `motherTopoFlightValid`。

## 编译与检查

请同时替换并重新编译两个模块；旧对象文件不能复用：

```bash
make -C CeeV0Lite clean
make -C CeeV0Lite all test
make -C CeeV0Analysis clean
make -C CeeV0Analysis check
# 配置 CEE、ROOT、FairRoot 环境后编译真实分析程序
make -C CeeV0Analysis all
```

可选独立数学对照（需要 NumPy/SciPy）：

```bash
make -C CeeV0Lite reference-test
```

默认运行命令不需要新增参数，By 沿用原来的带符号 T 单位。
检查覆盖解析加权结果、相关 PV 的抵消、S/L/cτ 误差、端点往返、异常状态、
TOF 与有场 tracking 集成，以及独立八参数广义最小二乘对照。

## 当前边界

本版面向 Λ、反 Λ、K⁰ₛ 这样的中性母粒子。有场时子粒子仍按均匀 By 螺线输运；
中性母粒子自身不受磁场弯转。没有增加 Ξ⁻ 的带电母粒子输运、质量约束、
Ξ/Ω 候选组合或整条级联链的反向重拟合。
八维状态是这些功能的基础之一，不等于已经完成它们。

当前环境已做独立数值测试和分析层无 ROOT 检查；实际 CEE/FairRoot 链接、
ROOT 文件读写、真实数据与模拟 pull/效率验证仍需在实验环境中进行。
