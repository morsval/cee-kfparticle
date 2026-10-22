# v0.7.0 通用 tracking 通道

编译：在仓库根或 CEErco 执行一次 `make -j4`（已配置 ROOT/CeeRoot/FairRoot）。独立回归 `make check`；需要 ROOT 的输出测试 `make -C CEErco/CeeV0Analysis channel-writer-test`。

用户主入口：`CeeKFParticleAnalysisMaker::Channels()` 配置 `channels` 和 `pid`，调用 `MakeChannels` 或 `RunChannels`。原 `Make/Run` 保留固定两体和零场 TOF；两套母粒子选择接口互不暗中套用。

`DaughterSpec::Track(role,{pdg,charge,mass})` 引用叶 track；`Composite(role,sourceChannel)` 引用前一步原始 7D 母粒子。每通道 `select` 在候选供下游使用前执行。自动拓扑排序，拒绝循环依赖和重复末级 track。相同假设的排列去重。候选或组合上限触发时报错，清空当前事件，不静默返回截断样本。

输入要求：实测整数电荷、物理动量、与位置同参考点的完整 Cartesian 协方差。BenoDst 的电荷符号/5D 协方差不可直接冒充这些量。Beno 示例要求用户提供转换 callback；已核对上游 83e58c2f523dcfc1644068afe46c347280e9db80，并通过真实头文件语法检查。

拟合：n-body 同时拟合一个共同 SV，几何 NDF=2n-3。按原始输入反复重线性化，完整数值导数传播母粒子与 daughter 边缘协方差。复合 daughter 保留独立能量误差。均匀 CEE By，零过程噪声；不包含材料模型、质量约束或任意场图。与已有二体核并存，不宣称同一实现或提速。

级联：Lambda 使用未加碰撞 PV 约束的 7D 母粒子进入 Omega，不能把 Lambda 约束到错误产生点。parent 更新后的 daughter 状态单独存储。未自动为每一级生成相对于实际 parent SV 的 8D 飞行量。

## 常用 cut 与输出

|量|查询|定义|
|---|---|---|
|daughter-PV DCA|`c.Daughter("he3").dcaToPvCm` 或 `c.daughters.at(i).dcaToPvCm`|拟合前 daughter 到碰撞 PV 的三维无符号距离，cm|
|共同顶点质量|`c.chi2/c.ndf`|n-body 几何 χ²/NDF|
|daughter 间距|`c.maximumDaughterSeparationCm`|最终 POCA 位置最大两两间距，非分别最小化的两线 DCA|
|质量/误差|`c.massGeV`, `c.sigmaMassGeV`|原始母粒子 GeV/c²|
|PV 相容性|`c.motherTopo.pvConstraintChi2/pvConstraintNdf`|检查 valid 后使用|
|飞行显著性|`c.flight.decayLengthCm/c.flight.sigmaDecayLengthCm`|有符号弧长；先检查 valid 与误差>0|
|S|`c.motherTopo.flightParameterS`|cm/(GeV/c)，不是显著性|
|sigmaS|`sqrt(c.motherTopo.covariance.m[7][7])`|8D 协方差第八维误差|
|末级 track|`c.leafTrackIds`|去重身份与完整追溯|

通用候选尚未单独提供 daughter χ²Primary。旧两体 Candidate 的该字段仍保留。不能把 DCA 当作显著性。

`CeeKFParticleChannelWriter` 写动态 daughter 角色/PDG/电荷/ID、7D 输入/拟合状态及可选 7×7 边缘协方差。关联键是 `(eventId,candidateId)`，完整保留已接受中间态。`AddColumn` 添加用户标量；完全自定义树可实现 `CeeKFParticleChannelSink`。跨 daughter 联合协方差没有另存。

示例 `macros/runCeeChannels.C` 展示 Lambda、Omega、He3+pi、d+p+pi。默认不假造核素 PID；chargeOnlyDemo=true 仅是显式的流程测试。

## 验证与边界

- 既有核心/分析回归，2/3-body、级联不同顶点、角色、排序、共享叶 track 和错误处理通过。
- 240 个 Lambda 高斯 toy：全部接受，mass pull mean=0.125775、width=0.999246、mean chi2=0.949843。不是三体/级联完整 pull 闭合。
- ROOT 动态 daughter 图、协方差、自定义列和事件 ID 回读；Cling macro 合成输入重建一个 Lambda。
- 未完成真实 CEE/FairRoot 数据闭合。零场 TOF 三体和级联仍不支持，明确报错。

用户手册已重排为 13 章，入口 `docs/manual_source/build_user_book_v070.py`；详细原理与 cut 表见 PDF。
