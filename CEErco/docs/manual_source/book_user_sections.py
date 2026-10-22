# User-facing chapters for the book. Executed by build_user_book_v070.py.
chapter('软件结构与阅读路线','第一次使用请依次阅读第 1–6 章；查 cut 和数据结构看第 7 章；算法原理在第 8–11 章。')
section('先分清输入、动量模式和衰变结构',
 P('CEE KFParticle 从选中的 daughter track 重建共同衰变顶点 SV（Secondary Vertex），合成母粒子，并可另存主顶点 PV（Primary Vertex）约束副本。POCA（Point of Closest Approach）是 track 到给定参考点的最近接近位置，并不是 SV 的别名。'),
 table([['独立配置','选择内容'],['输入格式','原始 CEE track + hit；经过适配的 BenoDst / 其他 DST'],['动量来源','tracking 动量；零场通过 TOF 迭代计算动量'],['衰变结构','两体；同一 SV 的三体；具有多个 SV 的级联'],['输出策略','固定两体 Writer；通用 ChannelWriter；用户自定义 Sink']],[120,CW-120]),
 P('现有两体 Finder 同时支持有磁 tracking 和零场 TOF。新 Channels 支持 tracking 动量下的自定义两体、三体和级联；不会把零场 TOF 两体的时间公式直接用于级联。BenoDst 接入需要完成真实状态与协方差映射，本书明确区分已实现接口与尚需用户验证的输入转换。'))
section('各类职责与数据流',
 table([['类 / 结构','作用'],['AnalysisEvent / AnalysisTrack','格式无关的值对象：PV、输入 track、质量与可选 PID'],['CeeKFParticleEventSource','逐事件读取接口；换输入格式时实现 Next'],['CeeKFParticleRootSource','读取现有 CEE ROOT / list，适配 TPC，必要时匹配 TOF'],['CeeKFParticleAnalysisMaker','用户分析入口：事件/track cut、PID、预选列表与调度'],['CeeKFParticleInterface','连接既有两体 CeeV0Lite::Finder'],['CeeKFParticleChannels','自定义通道依赖、组合、共同顶点与级联'],['VertexFitter / Transport','既有两体数值拟合、By 输运、POCA'],['ProductionVertex','产生独立 8D PV 约束状态，不覆盖原母粒子'],['RootWriter / ChannelWriter / Sink','树结构、变量选择、daughter 链接和输出生命周期']],[183,CW-183]),
 P('Make/Run 是原两体路径；MakeChannels/RunChannels 是通用路径。两者共用事件和 track 质量选择。通用 PID 使用 Channels().pid，每通道的母粒子 cut 使用 DecayChannel.select，不会自动调用旧 p/π mask 或 AcceptMother。'))
section('用户日常主要改哪里',
 table([['需求','修改位置'],['输入文件、磁场、事件数','macro 参数 / AnalysisConfig'],['PID、track cut、预选列表','用户 Maker 类或 callback'],['想重建哪些粒子及衰变','DecayChannel / DaughterSpec'],['中间态与末级母粒子 cut','各通道 select'],['输出哪些分支','Writer 配置或用户 Sink'],['切换 BenoDst','输入 adapter，不改拟合公式']],[160,CW-160]),
 P('推荐把自己的 UserAnalysis 和通道配置单独保存。macro 负责构造 Source、Maker、Writer 和启动事件循环；拟合核心通常不用改。小分析可以直接在 macro 中配置 callback，无需引入 STAR StMaker 框架。'))

chapter('编译与第一次运行','加载已有 CeeRoot 环境后，在仓库顶层执行一次 make；安装 CeeRoot 本身不是本项目 Makefile 的职责。')
section('一个 make 编译所需 CEE 组件',
 code('cd cee-kfparticle\nmake -j4'),
 table([['生成物','用途'],['CEErco/CeeV0Lite/libCeeV0Lite.a','独立核心静态库'],['CEErco/CeeV0Analysis/CeeV0Analysis.bin','既有 CLI 分析程序'],['CEErco/CeeV0Analysis/libCeeKFParticleAnalysis.so','供 ROOT macro 调用的共享库']],[290,CW-290]),
 P('在 CEErco 目录也可以执行 make。顶层入口传播并行编译参数，无需手动逐目录调用。root-config、VMCWORKDIR、FAIRROOTPATH 应来自你的实验环境。StRoot 仅作参考，不属于此构建目标。'),
 code('make lite     # standalone core\nmake check    # core + analysis tests\nmake legacy   # optional old Analyser.bin\nmake clean'),
 P('旧 Analyser.bin 不是新 macro 的运行依赖；需要时用 legacy 目标。旧 Makefile 中失效的 AnaMaker/ 路径已修正。完整实验库链接仍需在实际 CeeRoot/FairRoot 环境验证。'))
section('先运行已有两体 macro',
 code('cd CEErco/CeeV0Analysis\nroot -l -b -q \\\n \'macros/runCeeKFParticle.C("input.list","mothers.root",0.5)\''),
 P('ByTesla=0.5 使用 tracking 动量，0 使用零场 TOF。input 可以是单个 .root 或 .list/.txt；列表每行一个文件路径。输出默认不覆盖。PID 与 cut 数字是示例，必须按你的标定和分析工作点修改。先核对事件、track 和候选数，再判断质量谱。'))
section('通用通道 macro',
 code('root -l -b -q \\\n \'macros/runCeeChannels.C("input.list","omega.root",0.5,"omega",true)\''),
 P('末尾 true 显式开启仅按电荷组合的流程演示，不是已做标定 PID 的物理结果。正式运行用 false 并修改 ChannelPID。支持示例 lambda、omega、hypertriton2、hypertriton3。He³ 需要 q=+2 与物理动量，程序不猜测刚度转换。'))

chapter('输入层：原始 track/hit 与 BenoDst','输入格式与重建模式独立；同一事件适配接口可以供不同输入源使用。')
section('两条输入路径',
 P('原始 CEE 路径由 RootSource 读取 TPC、PV 和 TOF，必要时自行匹配；BenoDst 路径复用保存的 track、匹配与 PID，通过 adapter 形成 AnalysisEvent。用户已在 FairTask 中读取事件时，可以直接调用 Make/MakeChannels，不必再套一个文件 reader。'),
 P('AnalysisEvent 是值对象，不保留 ROOT 当前事件缓存的裸指针。预选列表只决定哪些 track 参与后续处理；不能混用 Beno 数组索引、原 trackId 和母粒子 candidateId。'))
reuse(25,26,27,28,30)
section('BenoDst 源码核对与字段映射',
 P('核对版本：CEESM/CeeRoot master，提交 83e58c2f523dcfc1644068afe46c347280e9db80。CeeBenoDstReader 依赖 FairRootManager，通过 init、event、tpcTrack(i) 访问已加载事件，不是直接接收文件路径的独立 reader。'),
 table([['字段','KFP 接入约定'],['primaryVertexIndex()','Beno Vertex 数组索引；先检查有效性'],['primaryVertex()/primaryVertexError()','PV 与 sigma，cm；构造协方差需平方'],['origin()/globMom()','输入参考位置与全局动量；不用 primMom 代替'],['trackId()','原始身份，不等于 Beno 数组下标'],['covMatrixIndex()','协方差数组索引；检查 tpcTrackIndex()==i'],['*BestMatchIndex()','对端 hit/track 数组索引，不是 MatchInfo 索引'],['nSigmaPion/Kaon/Proton()','缺失 sentinel 转 NaN；getter 存在不代表生产已填充'],['status()/trackChi2()','顶点关联标志 / 原始 χ²，不是收敛标志 / χ²/NDF']],[192,CW-192]),
 src='CeeBenoDst/edm/CeeBenoEvent.h、CeeBenoTpcTrack.h；reader/CeeBenoDstReader.h；maker/CeeBenoDstMaker.cxx。')
section('BenoDst 的 5D 协方差不能直接当 6D',
 P('CeeBenoTpcTrackCovMatrix 保存 a=(DCAxy,DCAz,φ,θ,q/p)、5 个 sigma、10 个相关系数。C5 对角为 sigma²，非对角为 rho·sigma_i·sigma_j。下三角相关项顺序：(1,0),(2,0),(2,1),(3,0),(3,1),(3,2),(4,0),(4,1),(4,2),(4,3)。'),
 eq(r'C_6=J_{6\times5}C_5J_{6\times5}^{T}',18),
 P('J 必须依据上游真正的参数定义、参考点与坐标约定计算。不能直接复制、补零或只保留位置对角。isValid() 允许 q/p sigma=0 作为无动量占位，所以它为 true 不等于可以用于有磁 tracking。'),
 P('examples/BenoDstAdapterExample.h 负责选中 track 的索引核对和 PV/PID 映射，并要求用户提供 BenoStateConverter 输出同参考点的物理状态与完整协方差。该示例已用真实 Beno 头文件做语法检查，未声称完整 DST 运行闭合。'),
 note('电荷与参考点','charge() 从 hit 数符号得到 ±1，无法单独表达 He³ 的 |q|=2。Maker 直接把 GetX/Y/Z、GetPx/Py/Pz 写入 origin/globMom，没有额外求 DCA。KFP 的默认 DCA 重参数化仍在合格输入之后进行。生产 Maker 的 nSigma 当前保留 sentinel，应核对实际数据是否后来填充。'),
 src='CeeBenoTpcTrackCovMatrix.h/.cxx；CeeBenoDstMaker::packBenoTpcTracks、fillBenoTpcTrackCovFromUpstream。')
section('BenoDst TOF 的匹配、坐标和时间',
 P('TpcItofMatchInfo / TpcEtofMatchInfo 保存 Beno track/hit 索引、匹配 χ²、TOF、beta。MatchInfo 的 hit 位置已对齐，而原 hit 分支可以仍为 ideal/raw lab，不能混用。BestMatchIndex 指向 hit，应按 trackIndex 与 hitIndex 找匹配记录。'),
 P('零场衰变重建需要到达时刻及误差，每轮扣除母粒子飞行时间并更新 daughter 路程。不能把按 primary 假设得到的 beta 直接当最终 daughter 动量；不能重复减去已经校正过的 T0。triggerTs 是 Unix 时间戳，不是本算法的碰撞时间原点。'),
 P('零场 Beno 适配还需明确时间基准、共享 T0、对齐位置、hit 冲突和缺失值。当前 Beno 示例只映射 tracking，不会隐藏这些尚需标定和验证的选择。'),
 src='CeeBenoTpcItofMatchInfo.h；CeeBenoTpcEtofMatchInfo.h；CeeBenoEvent.h；matcher/utils 时间与位置对齐。')

chapter('分析类：PID、预选 track 与 cut','先减少不需要的输入，再做组合；原始母粒子与 PV 副本的 cut 必须分开。')
reuse(8,9)
section('通用 PID 与核素假设',
 code('maker.Channels().pid =\n [](const AnalysisTrack& t, const ParticleHypothesis& h,\n    const AnalysisEvent& e) {\n   return MyPidDecision(e.eventId,t.input.id,h.pdg);\n };'),
 P('MyPidDecision 是用户实现的 PID 函数。ParticleHypothesis 包含 PDG、整数电荷和质量。允许一个 track 通过多个假设，粒子池按假设缓存，后续只组合通过 PID 的 track。AnalysisTrack 增加 nSigmaKaon；核素 PID 可放在用户类或外部映射中。'),
 P('BenoDst 没有提供本接口所需的 NDF / 拟合收敛标志时，chi2Ndf 保持 NaN，fitConverged 保持 false。应关闭相应基础 cut 或补齐真实信息，不得把 status() 硬解释成收敛。'))
section('每个通道独立选择母粒子',
 code('lambda.select =\n [](const ChannelCandidate& c,const AnalysisEvent&) {\n   if (!(c.massGeV>1.08 && c.massGeV<1.16)) return false;\n   if (!(c.chi2/c.ndf<20.)) return false;\n   return c.Daughter("proton").dcaToPvCm>0.1;\n };'),
 P('这些阈值仅展示写法。Lambda 的 select 在进入 Omega 组合之前执行，Omega 有自己的 select。Channels 不读取旧 lambdaCuts/k0Cuts，也不调用旧 AcceptMother；它复用的是磁场、参考点、输运、顶点和 PV 配置。缺失浮点值必须明确拒绝，防止 NaN 漏过。'))
reuse(10)

chapter('用户定义衰变：两体、三体与级联','daughter 的名字是角色，不是固定粒种。每个通道描述实际物理拓扑。')
section('He³+π 与 d+p+π 是两个独立通道',
 code('DecayChannel h2;\nh2.name="hypertriton2"; h2.motherPdg=1010010030;\nh2.daughters={\n DaughterSpec::Track("he3",{1000020030,2,2.8083916}),\n DaughterSpec::Track("pion",{-211,-1,0.13957039})};\nDecayChannel h3;\nh3.name="hypertriton3"; h3.motherPdg=1010010030;\nh3.daughters={\n DaughterSpec::Track("deuteron",{1000010020,1,1.8756129}),\n DaughterSpec::Track("proton",{2212,1,0.9382720813}),\n DaughterSpec::Track("pion",{-211,-1,0.13957039})};'),
 P('母粒子质量由 daughter 状态合成，motherPdg 仅作标记，不自动加质量约束。用户可统一替换质量常数。三体所有 daughter 同时参加一个 SV 拟合，不会产生人为的二体中间态。'))
section('先 Lambda，再 Lambda+K− 合成 Omega−',
 code('DecayChannel lambda;\nlambda.name="lambda"; lambda.motherPdg=3122;\nlambda.daughters={\n DaughterSpec::Track("proton",{2212,1,0.9382720813}),\n DaughterSpec::Track("pion",{-211,-1,0.13957039})};\nDecayChannel omega;\nomega.name="omega"; omega.motherPdg=3334;\nomega.daughters={\n DaughterSpec::Composite("lambda","lambda"),\n DaughterSpec::Track("kaon",{-321,-1,0.493677})};\nomega.fitCollisionPv=true;\nmaker.Channels().channels={lambda,omega};'),
 P('通道自动按依赖排序；未知依赖、重复名称、循环依赖会报错。每个候选记录所有末级 track ID，任何共享叶 track 的组合都被拒绝。相同假设的互换组合只算一次。'),
 P('用于 Omega 的 Lambda 是原始 7D 衰变态，保留能量误差和全部交叉协方差；不先约束到碰撞 PV。Lambda 的产生点是 Omega SV。fitCollisionPv 另建副本，即使生成过 Lambda 的 PV 副本，后续组合仍使用原 mother。'))
section('通用拟合的算法与范围',
 P('n 个 daughter 在当前 SV 附近分别求 POCA，以两个横向约束参与共同顶点解，再更新状态并合成母粒子。自由度为 2n−3，三体为 3。用完整求解的数值导数传播母粒子和 daughter 的边缘协方差；每个扰动重新求 POCA/SV。叶 track 是固定质量 6D 输入，中间态是含独立能量误差的 7D 输入。'),
 P('该实现不是 STAR 完整 SIMD Finder，也没有加速倍数承诺。只支持 tracking 动量（允许明确 By=0）；零场 TOF 的三体和级联时间耦合尚未实现，接口会明确拒绝。任意非均匀场与材料效应不在当前模型中。'))

chapter('macro 与 ROOT：完整事件循环','macro 装配 Source/Maker/Writer，重建核心不决定用户最终树结构。')
section('运行闭环',
 code('CeeKFParticleRootSource source(config);\nCeeKFParticleAnalysisMaker maker(config.finder);\n// configure maker.Channels().pid / channels\nCeeKFParticleChannelWriter writer("mothers.root");\nwriter.SaveCovariances(true);\nstd::string error;\nif (!source.Open(error)) throw std::runtime_error(error);\nif (!maker.RunChannels(source,writer,error))\n throw std::runtime_error(error);'),
 P('完整示例见 macros/runCeeChannels.C。在已有 Beno FairTask 循环中，用 adapter 构造 AnalysisEvent，然后 MakeChannels 与 writer.Write。Begin/Finish 放在初始化与结束，不要每事件重建输出文件。'))
section('母粒子树中保留 daughter 关系',
 table([['分支','含义'],['eventId + candidateId','事件内候选唯一键；跨事件不能只用 candidateId'],['channel / pdg / charge','通道、母粒子身份与 daughter 电荷和'],['daughterRole / daughterPdg / daughterCharge','动态 daughter 数组，按通道定义顺序'],['daughterTrackId / daughterCandidateId','叶 track 或中间候选 ID；不适用的一项为 -1'],['leafTrackIds','所有末级 track 身份，用于查重追溯'],['daughterDcaToPvCm','更新前 daughter 到碰撞 PV 的三维 DCA，cm'],['motherState / daughterInputState / daughterFittedState','7D 状态；daughter 每个占 7 项'],['*Covariance','可选 row-major 7×7；daughter 每个占 49 项']],[211,CW-211]),
 P('所有接受的中间候选保留在树中，Omega 可通过链接追溯 Lambda 再到 p/π。Omega 中 fitted Lambda 是 parent 更新后的副本，原 Lambda 行保持不变。跨 daughter 的完整联合协方差没有单独保存，不能把边缘矩阵当成独立矩阵再相加。'))
section('自主选择标量或完全自写输出',
 code('writer.ClearColumns();\nwriter.AddColumn("mass",\n [](const ChannelCandidate& c,const AnalysisEvent&) {\n   return c.massGeV;\n });\nwriter.AddColumn("dcaKaon",\n [](const ChannelCandidate& c,const AnalysisEvent&) {\n   return c.channel=="omega"\n    ? c.Daughter("kaon").dcaToPvCm\n    : std::numeric_limits<double>::quiet_NaN();\n });'),
 P('结构与关联键固定保留；标量可 AddColumn，协方差可开关。需要多个树或不同格式，实现 CeeKFParticleChannelSink::Begin/Write/Finish。PV 副本的标量也可 AddColumn，先检查 motherTopo.valid/flight.valid；默认 ChannelWriter 的数组为原始 7D，完整 8D 数组需自定义 Sink。'))
reuse(12)

chapter('数据结构与 cut 变量查表','通用 ChannelCandidate 与旧固定两体 Candidate 是不同 API，字段不能混用。')
section('通用候选查询表',
 table([['字段 / 方法','量、单位与阶段'],['massGeV / sigmaMassGeV','原母粒子质量及误差，GeV/c²'],['chi2 / ndf','共同 SV 几何 χ²；ndf=2n−3'],['maximumDaughterSeparationCm','最终 POCA 位置最大两两间距，cm'],['mother.positionCm / momentumGeV','原 SV 与动量，cm / GeV/c'],['mother.covariance','完整 7×7，含能量交叉项'],['daughters.at(i).dcaToPvCm','第 i 个 daughter 对碰撞 PV 的三维 DCA，cm'],['Daughter("he3").dcaToPvCm','按用户角色名取同一量'],['daughter.inputState / fittedState','进入 parent 前 / 后的状态及边缘协方差'],['motherTopo / flight','独立 PV 副本 / 飞行量；先检查 valid'],['leafTrackIds','所有末级输入 track，不等于直接 daughter 数']],[241,CW-241]),
 P('通用 API 尚未另存每个 daughter 的 χ²Primary，不能用 DCA 冒充。旧两体 Candidate 已有 daughterChi2FromPv。最大 POCA 间距也不是在所有 daughter 对上分别重新求两线 DCA 的最小值。'))
section('通用候选的 PV 与飞行量 cut',
 table([['字段 / 表达式','意义'],['motherTopo.valid / flight.valid','先检查拟合与派生量是否有效'],['motherTopo.pvConstraintChi2 / pvConstraintNdf','PV 新增约束 χ²/NDF；不要与几何 χ² 混用'],['flight.massGeV / sigmaMassGeV','PV 副本质量及误差，GeV/c²'],['flight.decayLengthCm','有符号弧长 L=S|p|，cm'],['flight.sigmaDecayLengthCm','完整误差传播得到 sigmaL，cm'],['L / sigmaL','飞行长度显著性；必须 sigmaL&gt;0'],['motherTopo.flightParameterS','飞行参数 S，cm/(GeV/c)，不是显著性'],['sqrt(motherTopo.covariance.m[7][7])','sigmaS；S/sigmaS 与 L/sigmaL 不必数值相同'],['flight.ctauCm / sigmaCtauCm','固有衰变长度及误差，cm']],[258,CW-258]),
 P('PV 副本可选生成；其默认位置可能在产生端，不能把 motherTopo.positionCm 当作原 SV。L 为带方向的飞行弧长；只按阈值选显著性前，先明确是否要求正向飞行。通用默认 macro 不添加 L/sigmaL cut。'))
section('已有两体的 89 个标量目录',
 P('以下表属于原两体 RootWriter。固定 [0]/[1] 仅约束旧 Candidate，不限制新 Channels。DcaProton/Pion 是按 PDG 查找的便利函数；通用角色可取 daughter1、daughter2、he3、deuteron 等任意唯一名字。'))
reuse(*range(14,25))
