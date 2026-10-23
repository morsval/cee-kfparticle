# BenoDst 输入验证（2026-10-23）

本次在 ROOT 6.40 下使用 CEESM/CeeRoot 提交
`83e58c2f523dcfc1644068afe46c347280e9db80` 的真实 Beno EDM 类生成字典、写入合成 ROOT 事件并回读。测试仅为 FairLogger 提供输出桩；Beno 字段、压缩存储和 ROOT 序列化均为上游代码。

- 顶层 `make check`：既有核心和分析回归通过。
- `make -C CEErco/CeeV0Analysis beno-source-test BENO_EDM_DIR=/path/CeeRoot/CeeBenoDst/edm`：零场 Λ、tracking Λ、固定局部面协方差、反向 globMom、公共 T0 误差、iTOF/eTOF 优先级、共享 hit、选中列表、非连续 trackId、协方差索引、文件列表、缺失分支及错误输入通过。
- 零场合成 Λ：质量约 1.11568 GeV/c²，配置下 sigmaM 约 0.00337957 GeV/c²。这是单个受控事件的结果，不是实验质量分辨率测量。
- T0 误差从 0.1 ns 改为 0：质量中心值不变，质量误差减小；MatchInfo 的已扣 T0 时间未重复扣除。
- 完整 Cling macro `runCeeKFParticle.C(...,"benodst","tof",0.,-1,true,true)`：1 个事件、2 条选中 track、1 个母粒子；ROOT 回读检查质量、daughter ID 100/107、metadata 中的 Beno cut 配置通过。
- `runCeeChannels.C(...,0.5,"lambda",false,-1,"benodst")`：读取同一套真实 EDM 的 tracking 输入并完成通用输出。
- 原始 RootSource / Writer 回归通过；tracking 不再要求读取 T0Hit 分支。

上游 `GetEntryMomentum` 可将 globMom 反向保存，而 FitParameter/VertexSeed 保留原拟合方向。零场按同一直线检查并由 TOF 定向；tracking 使用 FitParameter 方向及上游固定局部面上的误差表示。此约定针对上述生产源码，其他版本使用显式 StateConverter。

完整 CeeRoot/FairRoot 安装、真实 BenoDst 标定与大样本 pull 闭合尚未在本环境验证。TOF 路程误差模型沿用原实现，新增几何交叉项不等于对所有原始测量做联合拟合。核素需要用户明确物理电荷与刚度转换。
