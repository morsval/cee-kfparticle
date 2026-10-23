# v0.6.3：带电母粒子在 CEE By 中的产生顶点约束

基线为已发布的 v0.6.2。CEE 束流沿 z，磁场为 `(0,By,0)`，By 使用带符号的 T。
本版把带电母粒子的螺线 PV 约束和八维端点输运合入正式 C++ 核心。
原中性分支、七维候选及其选择保留。

## 核查结论

超氚 `³ΛH → ³He + π⁻` 的电荷依次为 `+1 → +2 + (-1)`。
其母粒子在磁场中也会弯转，不能将子轨迹按螺线处理后，又把母粒子当成直线。

v0.6.2 的 `FitProductionVertex()` 会拒绝非零母粒子电荷。
2026-09-29 的独立超氚测试用 `charged_topology.py` 另行完成了螺线约束，
没有修改正式包。因此此前问题是正式 C++ 缺少带电能力；
不能描述成正式 v0.6.2 已悄悄给超氚施加了直线约束。

## S 仍然表示什么

`S = L_arc / |p|`，其中 `L_arc` 是产生端到衰变端的有符号三维弧长。
在均匀磁场、无材料能损的模型中，动量模长不变，因而仍有：

```
L_arc = S * |p|
ctau  = S * sqrt(E*E - |p|*|p|)
```

带电有场时，`L_arc` 一般不等于两端之间的直线距离，
S 也不能再直接用 `p·(SV-PV)/p²` 作为最终答案。
代码只将这个投影用作迭代初值，最终同时求解螺线几何与 S。
S 的单位为 cm/(GeV/c)，L 和 cτ 的数值单位为 cm。
未施加母粒子质量约束，也没有输入已知寿命参与拟合。

## 螺线关系与误差

定义 `k=0.00299792458*q*By`、`a=k*S`，使用衰变端动量 pD：

```
rP = rD - A(S) pD
pP = R(S) pD,  R(S) = dA/dS

A(S) = [ sin(a)/k       0     (1-cos(a))/k ]
       [     0         S            0     ]
       [-(1-cos(a))/k  0        sin(a)/k   ]
```

实现中用 sinc 的展开避免 k 或 a 很小时相消，零场连续回到直线。
正电荷向 +z 运动、By>0 时，正向飞行向 -x 弯转。

联合输入仍为 `X=(rD,pD,E,PV)`，含完整 10×10 协方差及母粒子-PV 交叉块。
约束为 `g=rD-PV-A(S)pD=0`。S 是待求的自由参数，三条方程提供两个独立约束。
每轮在当前点线性化，沿 `v=(dA/dS)pD` 消去 ΔS，在垂直 v 的两个方向求更新。
所有迭代都使用同一份原先验，避免重复计入 PV。

误差也从同一联合先验传播：先求受约束状态的导数 Jx，再由
`H dx-v dS=0` 得到 `Js=vᵀ H Jx/|v|²`。
产生端的位置直接取联合拟合中的 PV 估计，动量旋转及 S 行一起映射成八维协方差。
输入事件 PV 对象保持不变。精确已知的 PV 可以作为零协方差输入，不强加虚构误差。

端点输运同时包含 `d(position)/dS` 和 `d(momentum)/dS`，包括 S 与所有变量的交叉项。
协方差采用 `C=L Lᵀ → (J L)(J L)ᵀ` 计算，减少强相关项相消造成的负方差。
这仍是局部一阶误差传播，中心值的螺线公式本身是精确的均匀场解。

## 如何调用

新重载要求显式传入电荷和磁场：

```cpp
ProductionVertexInput input;
input.mother = fittedMother;  // SV 处七维状态及边缘协方差
input.vertex = primaryVertex;
// 如有已知相关性，还需填写 input.motherVertexCovariance[7][3]。

ProductionState topo;
bool ok = FitProductionVertex(input, ProductionVertexConfig(),
                              geometricChi2, geometricNdf,
                              topo, +1, +0.5);  // 超氚，CEE By=+0.5 T
if (ok) {
  // 默认在产生端，state 中同时保存 charge 和 magneticFieldByTesla。
  FlightObservables flight;
  GetFlightObservables(topo, flight);
  ProductionState decay = topo;
  TransportToDecayVertex(decay);
}
```

旧的无磁场参数重载继续支持中性粒子；给它传非零电荷会返回状态 11，
避免调用者遗漏 By 后把带电轨迹默认为直线。真正的零场请显式传 `charge, 0.0`。
`charge` 是以 e 为单位的整数。轨迹输入必须是物理动量 p；氦三的刚度 p/2
不能直接当作动量，相关协方差也必须采用相同的物理动量约定。

## 输出与状态码

`ProductionState` 新增 `charge` 和 `magneticFieldByTesla`。
保存或重建对象时，必须连同八维参数、协方差和端点标志一起保留这两项。
ROOT 输出格式升级到 6，新增 `motherTopoCharge`、`motherTopoByTesla`；
原分支及下三角 36 项的顺序保留。

新增状态：

| 值 | 状态 | 含义 |
| --- | --- | --- |
| 11 | MissingMagneticField | 旧重载传入非零电荷但未提供 By |
| 12 | FlightRangeExceeded | 带电拟合超出允许弧长或弯转角范围 |

原状态 9 `UnsupportedCharge` 保留编号供旧结果识别，本版新入口不再用它拒绝非零电荷。
其他 0–8、10 的含义沿用 v0.6.2。无效状态使用 NaN，不能把失败解释成零寿命。

## 当前边界

- 核心支持均匀 By 下的带电母粒子及中性母粒子；材料效应、非均匀场和场误差未建模。
- 带电拟合求局部衰变飞行解。默认 `maximumAbsFlightPathCm=1000`，
  `maximumAbsBendingAngleRad=π`；到达边界返回失败，不猜测绕行圈数。
  该范围不等于已解决所有局部极值或子轨迹多交点问题。
- V0 Finder 与默认 CEE ROOT 分析仍按原物种配置运行。没有新增完整超氚 PID、
  氦三输入适配、Ξ/Ω 组合或整条级联的反向重拟合。
- 包内超氚示例直接调用通用两轨顶点核心和新带电 PV 接口，不经过旧 V0 的直线
  pointing/DCA 选择。不能把原中性母粒子的这些 cut 原样套到带电母粒子上。
- tracking 与 PV 的独立假设、TOF 的 PV 响应以及上游相关性的限制沿用 v0.6.2。
  有场超氚示例使用 tracking 动量；没有扩展零场 TOF 反推动量的物种逻辑。

## 编译与验证

请清理旧对象文件并重新编译两个模块，不能混用旧 `ProductionState` 的内存布局：

```bash
make -C CeeV0Lite clean
make -C CeeV0Lite all test
make -C CeeV0Lite reference-test  # 需要 NumPy/SciPy
make -C CeeV0Analysis clean
make -C CeeV0Analysis check
# 配好 CEE、ROOT、FairRoot 后：
make -C CeeV0Analysis all
```

`tests/test_charged_production.cxx` 覆盖正负电荷/磁场/S、完整端点 Jacobian、
相关 PV、零方差 PV、零场极限及失败状态。
`tests/verify_charged_reference.py` 用独立八参数广义最小二乘和洛伦兹力 ODE 核对。
`examples/hypertriton_by` 则把原独立超氚测试接到本版正式 C++ 核心。
实际 CEE/ROOT 链接、ROOT 文件读写及探测器级效率、pull 闭合仍需实验环境验证。
