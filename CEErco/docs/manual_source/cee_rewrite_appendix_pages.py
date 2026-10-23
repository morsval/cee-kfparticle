page('源码导航：生产链的十个实现文件',
    P('下表把文件名与正文的计算阶段对应起来。CEE 路径均相对 CEErco；这里列的是默认新版分析的实现文件。表中 [A] 等编号也用于各页页脚。'),
    table([['编号／目录／文件','职责与主要入口'],
      ['[C] CeeV0Analysis/src<br/>main.cxx','解析命令行；应用模式、预设与单项覆盖；启动 Analyzer。'],
      ['[C] CeeV0Analysis/src<br/>CeeV0AnalysisConfig.cxx','AnalysisConfig、ConfigureDataMode、ApplyCutPreset；建立与验证配置。'],
      ['[A] CeeV0Analysis/src<br/>CeeV0Analyzer.cxx','ROOT 输入输出、PV 与track准备、匹配、事件循环、FillRow。'],
      ['[G] CeeV0Analysis/src<br/>CeeV0AnalysisAlgorithms.cxx','ProjectionDistanceCm、AssignUniqueTofHits、T0 汇总；纯算法。'],
      ['[F] CeeV0Lite/src<br/>CeeV0LiteFinder.cxx','Reconstruct、BuildCandidate、SolvePair；组合、TOF 自洽、拓扑与 cut。'],
      ['[V] CeeV0Lite/src<br/>CeeV0LiteVertexFitter.cxx','Fit、RecalculateMotherCovariance；共同顶点、母粒子与 24 次扰动。'],
      ['[T] CeeV0Lite/src<br/>CeeV0LiteTransport.cxx','PropagateState、TransportToClosestApproach；直线／均匀 By 输运与最近点。'],
      ['[Y] CeeV0Lite/src<br/>CeeV0LiteTypes.cxx','数据默认值、EstimateTofMomentum、MotherMassAndError。'],
      ['[PV8] CeeV0Lite/src<br/>CeeV0LiteProductionVertex.cxx','联合 PV 约束、S、端点输运与飞行量误差。'],['[M] CeeV0Lite/src<br/>CeeV0LiteMath.cxx','向量矩阵运算、投影、求逆、有限性及半正定检查。']],[221,278]),
    src='[D] 适配器是 CeeV0Lite/include/CeeV0LiteCeeAdapter.h，主要逻辑定义在头文件内。')

page('源码导航：旧分析链与数据对象',
    P('这些文件均位于 CEErco 顶层。旧分析程序与新版 CeeV0Analysis 是两条不同入口；新版依然读取 CEE 数据对象，但不调用旧 MyTrack／MyLambda 进行候选重建。'),
    table([['文件','当前阅读意义'],['analysis.cxx','旧程序入口，创建 AnaMaker。'],['AnaMaker.cxx','旧事件循环、逐track寻找 TOF、调用 MyTrack／MyLambda、保存旧输出。'],['MyTrack.cxx','旧track处理；从 PV 到 TOF 的距离和时间开始估算动量。'],['MyLambda.cxx','旧 Λ 几何与动量重建；与新版对称 TOF–SV 自洽算法分别阅读。'],['CeeTpcTrack.cxx','TPC track数据对象；五维跟踪参数与六维顶点种子接口。'],['CeeTpcVertex.cxx','重建顶点数据对象；坐标、误差、rank 等由上游填入。'],['CeeT0Hit.cxx','T0 时间对象；普通模拟模式不读取它。']],[170,329]),
    Spacer(1,12),
    P('CeeTpcTrack.h 和 CeeTpcVertex.h 同样重要：ROOT 类版本、成员布局、访问器以及 Clear 行为可能定义在头文件中。判断输入是否可用，要同时检查实际 ROOT 文件、数据类字典和上游写入逻辑。'),
    P('legacy 预设不会切换到这条旧分析链。它只修改新版 Finder 的一组筛选阈值，其余输入适配、TOF 分配、顶点求解与误差检查仍然执行。'),
    src='[D] 适配器；CEE_INPUT_PRODUCER.md；CEErco 顶层数据类与旧分析文件。')

page('源码导航：测试与集成示例',
    P('这些实现文件只有在相应测试或示例入口被执行时才运行。了解每项测试覆盖的内容，可以判断某个“通过”结论到底验证了什么。'),
    table([['目录／文件','覆盖范围'],['CeeV0Lite/tests<br/>test_math_transport.cxx','矩阵运算、直线／有场输运、最近点选择、参考位置等数值性质。'],['CeeV0Lite/tests<br/>test_mother_covariance.cxx','七维协方差、扰动诊断、子粒子交换、公共 T0 传播及小误差玩具样本。'],['CeeV0Lite/tests<br/>test_end_to_end.cxx','核心零场 TOF／跟踪动量通道、筛选和收敛的合成输入检查。'],['CeeV0Lite/tests<br/>test_field_by.cxx','独立洛伦兹力积分，核对 By 方向、极性、协方差与多圈最近点。'],['CeeV0Lite/tests<br/>test_cee_adapter.cxx','种子状态、协方差、顶点误差与备用输入约定。'],['CeeV0Analysis/tests<br/>test_algorithms.cxx','探测器投影、全局唯一分配与 T0 汇总。'],['CeeV0Analysis/tests<br/>test_config.cxx','预设值、模式选择和时间配置政策。'],['CeeV0Lite/examples<br/>AnaMakerIntegrationSnippet.cxx','展示已有 Maker 如何准备核心输入；示例不是默认主程序。']],[237,262]),
    Spacer(1,10),
    P('另有 CeeV0Analysis/tests/test_field_cli.py 检查新旧场参数。合成输入检查不能替代 ROOT 字典、真实track误差、错配率与接受度验证。'),
    src='本轮已运行核心测试、分析层检查与 By 参数测试；尚未进行实际 ROOT I/O 验证。')

page('源码导航：STAR 接口与核心状态',
    P('本页和后两页路径均相对 StRoot。它们用于对照；CEE 的独立轻量核心不需要链接这些 STAR 类。'),
    table([['目录／文件','角色'],['StKFParticleAnalysisMaker<br/>StKFParticleAnalysisMaker.cxx','[S3] 分析层事件入口，调用接口并处理候选及输出。'],['StKFParticleAnalysisMaker<br/>StKFParticleInterface.cxx','[S3] 读取track、准备 PID 假设，连接 TopoReconstructor。'],['StKFParticleAnalysisMaker<br/>StKFParticlePerformanceInterface.cxx','重建结果与模拟真值的性能评估桥接；仅列于本表一次。'],['KFParticle<br/>KFParticleTopoReconstructor.cxx','[S3] 组织track、PV 与 Finder；含条件编译外部 tracker 分支。'],['KFParticle<br/>KFParticleFinder.cxx','[S2] 按通道组合，调用粒子构造与选择，管理候选副本。'],['KFParticle<br/>KFParticleBase.cxx','[S1] 标量状态初始化、测量、共同顶点、输运和产生顶点约束。'],['KFParticle<br/>KFParticle.cxx','标量 KFParticle 的场模型、距离等接口封装。'],['KFParticle<br/>KFParticleBaseSIMD.cxx','[S1] 批量状态与数学实现；同名函数不保证逐项实现完全一致。']],[254,245]),
    src='STAR 基线 Git HEAD：b6a9a93b574b5e6ca78f28f2381add3f25410263。')

page('源码导航：STAR 容器与辅助对象',
    table([['KFParticle 目录中的文件','角色'],['KFParticleSIMD.cxx','批量粒子接口及具体场模型封装。'],['KFPTrack.cxx','单track六参数、21 个压缩协方差元素、电荷与编号。'],['KFPTrackVector.cxx','适合批量处理的track容器与track提取。'],['KFPVertex.cxx','主顶点参数、协方差与统计信息容器。'],['KFVertex.cxx','与 KFParticle 顶点表示相连的封装与操作。'],['KFParticlePVReconstructor.cxx','按相应流程由track重建碰撞主顶点。'],['KFPEmcCluster.cxx','电磁量能器对象支持，不是 Λ → pπ 主线的必要测量。'],['KFParticleDatabase.cxx','供 Finder 等读取的粒子性质与质量常量。']],[223,276]),
    Spacer(1,14),
    P('单粒子类、容器和 Finder 各有职责。KFPTrack 保存输入track，并不承担 CEE 的 TOF 闭环；KFParticleDatabase 提供标准粒子性质，也不表示每个候选的质量或寿命都已被强制固定为标准值。'),
    P('查阅一段构造过程时，先找调用者传入了哪种粒子假设，再看是否显式调用质量约束或产生顶点约束。依据对象名称推断约束历史，容易把初始化、构造和额外约束混为一谈。'),
    src='[S1][S2][S3]；同一文件的角色已在相邻表中去重。')

page('源码导航：STAR 性能评估与文件计数',
    table([['KFParticlePerformance 内文件','角色'],['KFMCParticle.cxx','模拟粒子及其重建关联的记录。'],['KFMCVertex.cxx','模拟顶点记录。'],['KFPHistogram/KFPHistogramSet.cxx','性能直方图集合与相关操作。'],['KFPartEfficiencies.cxx','效率统计。'],['KFParticlePerformanceBase.cxx','性能评估的公共基础。'],['KFTopoPerformance.cxx','[S4] 拓扑重建与真值对照、参数和寿命等 QA。']],[262,237]),
    Spacer(1,13),
    P('前述检索表与新增八维章节共覆盖 50 个不重复的 .cxx：新版 CEE 生产 10 个、旧 CEE 与数据类 7 个、测试和示例 11 个、STAR 接口与核心 16 个、性能对象 6 个。性能接口已经列在 STAR 接口表中，不再重复计数。'),
    P('这里的 50 个是本书纳入的实现文件检索范围，不宣称枚举了整个仓库所有源文件。头文件中的适配器、ROOT 宏、构建文件和测试桩同样可能影响运行，应按正文中的具体问题继续查阅。'),
    note('性能评估的输出与分析选择分开','某个 QA 类会读取 cτ 或误差，说明这些量可以用于检验性能；不能反推普通分析树必定保存了它们，也不能反推它们已参与候选 cut。'),
    src='[S4] KFTopoPerformance；本表中的性能接口不与上一表重复计数。')

page('矩阵索引速查：21、28 与 36 个元素',
    P('实数对称矩阵的两半互为转置，所以输出可以只保存下三角。必须先确认状态顺序，再把行列映射成一维数组下标。内部 Mat6／Mat7／Mat8 为完整矩阵；压缩是相应输入或输出接口的存储约定。'),
    eq(r'k=\frac{i(i+1)}{2}+j,\qquad 0\leq j\leq i'),
    table([['矩阵','状态顺序','对角下标'],['6×6，21 项','x,y,z,pₓ,pᵧ,p_z','0,2,5,9,14,20'],['CEE 7×7，28 项','SVx,SVy,SVz,Pₓ,Pᵧ,P_z,E','0,2,5,9,14,20,27'],['CEE/STAR 8×8，36 项','x,y,z,pₓ,pᵧ,p_z,E,S','0,2,5,9,14,20,27,35']],[111,250,138]),
    Spacer(1,10),
    P('在 CEE 的 motherCovariance[28] 中，21、22、23 对应 E 与三个 SV 坐标的协方差；24、25、26 对应 E 与三个动量分量的协方差，27 是 E 的方差。CEE motherTopoCovariance 与 STAR 的 28-34 则是 S 与前七维的交叉项，35 是 S 的方差。'),
    code('int k = (i >= j) ? i*(i+1)/2+j : j*(j+1)/2+i;'),
    P('输出有效性标志为假时，应先处理无效结果，不能将数组中的 NaN 当作零。原七维数组只有 28 项，不能按八维的 36 项布局读取。'),
    src='[A] FillRow；[Y] MotherState；[S1] KFParticleBase.h 的 fC[36]。')

page('符号和常用术语速查',
    table([['符号／术语','本书中的含义'],['PV、SV','碰撞主顶点、当前母粒子的衰变顶点。'],['TPC、TOF、T0','时间投影室；飞行时间探测器；事件起始时间测量或时间原点。'],['X、y、C','输入状态、输出状态、协方差；维度随具体映射说明。'],['pᵢ、P、E、M','子粒子动量、母粒子总动量、总能量、重建不变质量。'],['u、n、e₁/e₂','track单位方向、两轨共同法向、垂直某个方向的两个单位基。'],['s／signedPathCm','CEE 单次输运的有符号弧长，单位 cm。'],['S、dS','S 为保存的飞行参数；STAR dS 为单次输运参数。均用长度／动量。'],['L、L_parallel','原 CEE 用距离／沿动量投影；新副本的 L=Sp 为有符号长度。'],['t_decay、cτ','实验室系衰变时刻、固有衰变长度；二者不是同一量。'],['POCA、DCA','track上的最近位置、相应最近距离；需明确比较对象。'],['Jacobian、Kalman','导数矩阵、利用残差与误差更新状态的方法。'],['NDF、QA、pull','有效自由度、质量与性能检查、按重建误差标准化的残差。']],[131,368]),
    src='英文函数名和真实输出字段保留原写法；正文术语以本表为准。')

page('版本依据与本次修订范围',
    P('本版对应 CEE KFParticle v0.7.0，以先前交付 v0.6.5 为基线。新增 Maker / Interface、可替换 reader、用户 PID / cut 与自定义 ROOT 树。旧 CLI 格式仍为 6；新树按用户选列保存。有磁 DCA 参考和所有核心拟合公式保持不变。'),
    P('v0.6.5 已将 toy 的 PV DCA 参考约定接入正式 Finder；本版继续保留。新增的 Maker 在组对前执行用户 PID 和 track 选择，在重建后执行母粒子选择，并提供按角色查询 DCA 和自定义输出接口。'),
    sub('源码内容指纹'),
    code(SNAPSHOT_HASH),
    P('指纹范围：CeeV0Lite 与 CeeV0Analysis 的 include/src 中 .h/.cxx，以及 CEErco 顶层 CeeTpcTrack、CeeTpcVertex、CeeT0Hit 的 .h/.cxx。按相对路径排序，将路径、零字节、文件字节依次拼接后计算 SHA-256。','small'),
    sub('已完成的验证与仍需验证的部分'),
    P('九项核心测试通过（新增默认参考回归）；分析层算法、配置、语法及真实 CLI 参数解析检查通过。另用独立的八参数广义最小二乘拟合，核对新状态、完整协方差和 χ²。测试还覆盖电荷和磁场正负号、洛伦兹力 ODE、公共平移抵消、零方差 PV 和端点往返；超氚示例直接调用正式 C++ 核心。'),
    P('本次新增 Maker 同输入数值一致性、PID / 子集与生命周期测试，并完成真实 ROOT 写出读回及 ROOT 加合成 CEE 形状字典的多文件 reader 测试。尚未完成真实 CEE/FairRoot 字典链接和实验数据闭合；合成字典不能替代实验 EDM。局部线性误差仍需用 pull 检验。'),
    src='2026-10-22；docs/AnalysisMaker_v0.7.0.md；源码指纹对应本次交付代码。')
