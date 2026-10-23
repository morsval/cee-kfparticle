# 母粒子 / daughter cut 变量表（v0.7.0）

下面 89 个名字可直接用于 `writer.SetVariables({...})`，或 `CandidateValue(name,candidate,event)`。除固定的 `eventId`（int64）和 `candidateIndex`（int32），可选标量都写为 double。`DcaProton(c)`、`DcaPion(c)`、`Chi2PrimaryProton(c)`、`Chi2PrimaryPion(c)` 另有直接函数接口。

`mother` 是原 7D 状态，`motherTopo` 是附加 PV 约束 8D 副本；不能混用两阶段的质量、长度和 chi2。PV 约束失败时，相应数值由查询器写为 NaN，状态标志仍保留。S 是飞行参数，不是统计显著度。下表列出的无量纲比值也不是信号质量谱的 S/sqrt(S+B)。

|输出/查询名|单位|阶段|含义|
|---|---|---|---|
|`candidateType`|1|标识|母粒子 PDG 类型：3122、-3122、310|
|`trackId0`|1|标识|daughter 0 的原始输入 track ID|
|`trackId1`|1|标识|daughter 1 的原始输入 track ID|
|`daughterPdg0`|1|标识|daughter 0 的质量假设 PDG|
|`daughterPdg1`|1|标识|daughter 1 的质量假设 PDG|
|`massGeV`|GeV/c^2|约束前 / 输入|原 7D 母粒子不变质量|
|`sigmaMassGeV`|GeV/c^2|约束前 / 输入|原 7D 质量误差|
|`daughterDcaCm`|cm|约束前 / 输入|共同顶点拟合中的两条 daughter 相互 DCA|
|`geometricChi2`|1|约束前 / 输入|共同 SV 拟合原始 chi2|
|`geometricNdf`|1|约束前 / 输入|共同 SV 拟合自由度|
|`geometricChi2Ndf`|1|约束前 / 输入|共同 SV 拟合 chi2 / NDF|
|`topologicalChi2`|1|约束前 / 输入|原候选几何 chi2 + 原 PV 相容性 chi2|
|`topologicalNdf`|1|约束前 / 输入|原候选几何 NDF + 2|
|`topologicalChi2Ndf`|1|约束前 / 输入|原候选合计 chi2 / NDF；不是 8D 副本|
|`primaryVertexConstraintChi2`|1|约束前 / 输入|原母粒子与 PV 相容性的纯 chi2|
|`primaryVertexConstraintNdf`|1|约束前 / 输入|原 PV 相容性检验自由度|
|`primaryVertexConstraintChi2Ndf`|1|约束前 / 输入|原 PV 相容性 chi2 / NDF|
|`decayLengthCm`|cm|约束前 / 输入|约束前无符号直线距离 &#124;SV-PV&#124;|
|`sigmaDecayLengthCm`|cm|约束前 / 输入|上述直线长度的传播误差|
|`decayLengthSignificance`|1|约束前 / 输入|约束前无符号直线 L / sigmaL|
|`pointingCosine`|1|约束前 / 输入|SV-PV 与原母粒子动量的夹角余弦|
|`DcaProton`|cm|约束前 / 输入|质子或反质子输入 track 到 PV 的三维 DCA|
|`DcaPion`|cm|约束前 / 输入|唯一 pion 输入 track 到 PV 的三维 DCA|
|`Chi2PrimaryProton`|1|约束前 / 输入|质子或反质子到 PV 的原始 chi2|
|`Chi2PrimaryPion`|1|约束前 / 输入|唯一 pion 到 PV 的原始 chi2|
|`motherDcaToPvCm`|cm|约束前 / 输入|原中性母粒子直线到 PV 的三维 DCA（派生量）|
|`decayTimeNs`|ns|约束前 / 输入|零场 TOF 自洽解中的衰变时刻|
|`daughterDcaToPv0Cm`|cm|约束前 / 输入|指定 daughter 输入状态到 PV 的三维 DCA|
|`daughterChi2FromPv0`|1|约束前 / 输入|指定 daughter 到 PV 的原始 chi2|
|`daughterPvSignificance0`|1|约束前 / 输入|指定 daughter 的 sqrt(chi2 to PV)|
|`daughter0PxGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`daughter0PyGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`daughter0PzGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`daughterDcaToPv1Cm`|cm|约束前 / 输入|指定 daughter 输入状态到 PV 的三维 DCA|
|`daughterChi2FromPv1`|1|约束前 / 输入|指定 daughter 到 PV 的原始 chi2|
|`daughterPvSignificance1`|1|约束前 / 输入|指定 daughter 的 sqrt(chi2 to PV)|
|`daughter1PxGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`daughter1PyGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`daughter1PzGeV`|GeV/c|约束前 / 输入|共同顶点 Kalman 更新后指定 daughter 的动量分量|
|`svXCm`|cm|约束前 / 输入|未加母粒子 PV 约束的 SV 坐标|
|`pvXCm`|cm|事件|输入 PV 坐标|
|`pxGeV`|GeV/c|约束前 / 输入|原母粒子动量分量|
|`svYCm`|cm|约束前 / 输入|未加母粒子 PV 约束的 SV 坐标|
|`pvYCm`|cm|事件|输入 PV 坐标|
|`pyGeV`|GeV/c|约束前 / 输入|原母粒子动量分量|
|`svZCm`|cm|约束前 / 输入|未加母粒子 PV 约束的 SV 坐标|
|`pvZCm`|cm|事件|输入 PV 坐标|
|`pzGeV`|GeV/c|约束前 / 输入|原母粒子动量分量|
|`energyGeV`|GeV|约束前 / 输入|原母粒子总能量|
|`pGeV`|GeV/c|约束前 / 输入|原母粒子动量模长|
|`ptGeV`|GeV/c|约束前 / 输入|相对 z 束流轴的横动量|
|`rapidity`|1|约束前 / 输入|沿 z 束流轴的快度|
|`motherCovarianceValid`|1|状态|原 7D 协方差是否有效|
|`motherCovarianceStatus`|1|状态|原协方差状态枚举|
|`motherJacobianComplete`|1|状态|母粒子数值 Jacobian 各列是否成功|
|`covarianceDependentCutsEvaluated`|1|状态|依赖原协方差的 cut 是否实际执行|
|`tofIterations`|1|状态|TOF 外层迭代次数|
|`tofConverged`|1|状态|TOF 外层是否收敛|
|`motherTopoValid`|1|状态|8D 母粒子 PV 约束是否成功|
|`motherTopoStatus`|1|状态|PV 约束状态枚举；失败原因见 ProductionFitStatus|
|`motherTopoCorrelationModel`|1|状态|PV 与母粒子关联模型枚举|
|`motherTopoAtProductionVertex`|1|状态|当前 8D 状态是否位于产生端点|
|`motherTopoIterations`|1|状态|PV 约束内部迭代次数|
|`motherTopoCharge`|1|状态|母粒子电荷，以 e 为单位|
|`topoFlightValid`|1|状态|8D 飞行长度和误差派生量是否有效|
|`motherTopoByTesla`|T|PV 约束后|约束状态对应的带符号 By|
|`motherTopoPvChi2`|1|PV 约束后|8D 拟合的纯 PV 约束 chi2|
|`motherTopoPvNdf`|1|PV 约束后|8D PV 约束自由度|
|`motherTopoPvChi2Ndf`|1|PV 约束后|8D 纯 PV 约束 chi2 / NDF|
|`motherTopoChi2`|1|PV 约束后|8D 几何与 PV 合计 chi2|
|`motherTopoNdf`|1|PV 约束后|8D 合计自由度|
|`motherTopoChi2Ndf`|1|PV 约束后|8D 合计 chi2 / NDF|
|`flightS`|cm/(GeV/c)|PV 约束后|有符号飞行参数 S = Larc / &#124;p&#124;|
|`sigmaS`|cm/(GeV/c)|PV 约束后|sqrt(C_SS)，S 的拟合误差|
|`SOverSigmaS`|1|PV 约束后|有符号 S / sigmaS；不是质量谱显著度|
|`topoEnergyGeV`|GeV|PV 约束后|PV 约束后母粒子能量|
|`topoPxGeV`|GeV/c|PV 约束后|约束后当前端点的母粒子动量分量|
|`topoXCm`|cm|PV 约束后|约束后当前端点坐标；通常在 PV，不能当 SV|
|`topoPyGeV`|GeV/c|PV 约束后|约束后当前端点的母粒子动量分量|
|`topoYCm`|cm|PV 约束后|约束后当前端点坐标；通常在 PV，不能当 SV|
|`topoPzGeV`|GeV/c|PV 约束后|约束后当前端点的母粒子动量分量|
|`topoZCm`|cm|PV 约束后|约束后当前端点坐标；通常在 PV，不能当 SV|
|`topoMassGeV`|GeV/c^2|PV 约束后|PV 约束后母粒子不变质量|
|`topoSigmaMassGeV`|GeV/c^2|PV 约束后|PV 约束后质量误差|
|`topoLCm`|cm|PV 约束后|PV 约束后有符号弧长 S &#124;p&#124;|
|`topoSigmaLCm`|cm|PV 约束后|弧长完整 8D 误差，含动量与 S 交叉项|
|`LOverSigmaL`|1|PV 约束后|PV 约束后有符号弧长 L / sigmaL|
|`ctauCm`|cm|PV 约束后|本征衰变长度 S m；不是 ns 时间|
|`sigmaCtauCm`|cm|PV 约束后|本征衰变长度的完整传播误差|

## 完整对象、矩阵和诊断

|C++ 成员或接口|可用内容|说明|
|---|---|---|
|`candidate.mother`|位置、动量、能量、7×7 协方差|顺序 x,y,z,px,py,pz,E；看 covarianceValid|
|`candidate.motherTopo`|上述量 + S、8×8 协方差、电荷、By、端点标志|先检查 valid；当前状态通常位于 PV|
|`candidate.topoFlight`|质量/误差、弧长/误差、ctau/误差|先检查 motherTopo.valid 与 topoFlight.valid|
|`candidate.decayVertexCovarianceCm2`|原 SV 的 3×3 协方差|母粒子 PV 约束前|
|`candidate.daughterAtVertex[0/1]`|各 daughter 的位置、动量、质量、PDG、电荷、6×6 协方差|共同顶点更新后的状态；不要与输入 track–PV DCA 的定义混淆|
|`candidate.mother.jacobian*`|扰动拟合次数、成功数、失败列/方向、失败状态及 bitmask|数值传播 QA；不是物理信号 cut 的默认建议|
|`candidate.daughterDeviationFromPv[0/1]`|旧版兼容别名，等于 sqrt(chi2)|新代码用 daughterChi2FromPv 或 daughterPvSignificance|
|`maker.GetSelectedTracks()`|送入 Finder 的输入 track、ID、PID mask 等|Finder 内部参考转换诊断另外查询|
|`maker.GetInterface().GetReferenceDiagnostics()`|每条 track 的参考转换状态、DCA、路径及收敛信息|用于判断默认 track–PV DCA 参考转换是否成功|

DcaProton/Pion 使用绝对 PDG 识别角色，支持反 Lambda；K0S 有两个 pion，DcaPion 因角色不唯一返回 NaN，使用显式 daughter 0/1 变量。完整 C++ 表达式保存在同目录 CSV 的最后两列。新增 motherDcaToPvCm 仅适用于当前 Finder 的中性母粒子，不为带电母粒子假造直线 DCA。
