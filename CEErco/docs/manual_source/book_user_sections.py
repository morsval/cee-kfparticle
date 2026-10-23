# User-facing chapters for the book. Executed by build_user_book_v070.py.
chapter('软件结构与阅读路线','第一次使用请依次阅读第 1–6 章；查 cut 和数据结构看第 7 章；算法原理在第 8–11 章。')
section('先分清输入、动量模式和衰变结构',
 P('CEE KFParticle 从选中的 daughter track 重建共同衰变顶点 SV（Secondary Vertex），合成母粒子，并可另存主顶点 PV（Primary Vertex）约束副本。POCA（Point of Closest Approach）是 track 到给定参考点的最近接近位置，并不是 SV 的别名。'),
 table([['独立配置','选择内容'],['输入格式','原始 CEE track + hit；BenoDst ROOT / list；自定义 DST Source'],['动量来源','tracking 动量；零场通过 TOF 迭代计算动量'],['衰变结构','两体；同一 SV 的三体；具有多个 SV 的级联'],['输出策略','固定两体 Writer；通用 ChannelWriter；用户自定义 Sink']],[120,CW-120]),
 P('现有两体 Finder 同时支持有磁 tracking 和零场 TOF。新 Channels 支持 tracking 动量下的自定义两体、三体和级联；不会把零场 TOF 两体的时间公式直接用于级联。BenoSource 已支持当前 BenoDst 的零场 TOF 输入及单位电荷 tracking 输入；时间和参数约定在第 3 章说明。核素或其他上游参数约定通过 StateConverter 明确扩展。'))
section('各类职责与数据流',
 table([['类 / 结构','作用'],['AnalysisEvent / AnalysisTrack','格式无关的值对象：PV、输入 track、质量与可选 PID'],['CeeKFParticleEventSource','逐事件读取接口；换输入格式时实现 Next'],['CeeKFParticleRootSource','读取原始 CEE ROOT / list，适配 TPC，必要时匹配 TOF'],['CeeKFParticleBenoSource','读取 BenoDst ROOT / list，复用保存的匹配及时间'],['ConfigureInputModes / OpenCeeInput','独立选择输入格式、动量来源，并打开对应 Source'],['CeeKFParticleAnalysisMaker','用户分析入口：事件/track cut、PID、预选列表与调度'],['CeeKFParticleInterface','连接既有两体 CeeV0Lite::Finder'],['CeeKFParticleChannels','自定义通道依赖、组合、共同顶点与级联'],['VertexFitter / Transport','既有两体数值拟合、By 输运、POCA'],['ProductionVertex','产生独立 8D PV 约束状态，不覆盖原母粒子'],['RootWriter / ChannelWriter / Sink','树结构、变量选择、daughter 链接和输出生命周期']],[183,CW-183]),
 P('Make/Run 是原两体路径；MakeChannels/RunChannels 是通用路径。两者共用事件和 track 质量选择。通用 PID 使用 Channels().pid，每通道的母粒子 cut 使用 DecayChannel.select，不会自动调用旧 p/π mask 或 AcceptMother。'))
section('用户日常主要改哪里',
 table([['需求','修改位置'],['输入文件、磁场、事件数','macro 参数 / AnalysisConfig'],['PID、track cut、预选列表','用户 Maker 类或 callback'],['想重建哪些粒子及衰变','DecayChannel / DaughterSpec'],['中间态与末级母粒子 cut','各通道 select'],['输出哪些分支','Writer 配置或用户 Sink'],['切换 BenoDst','macro 的 inputFormat；不改拟合公式']],[160,CW-160]),
 P('推荐把自己的 UserAnalysis 和通道配置单独保存。macro 负责构造 Source、Maker、Writer 和启动事件循环；拟合核心通常不用改。小分析可以直接在 macro 中配置 callback，无需引入 STAR StMaker 框架。'))

chapter('编译与第一次运行','加载已有 CeeRoot 环境后，在仓库顶层执行一次 make；安装 CeeRoot 本身不是本项目 Makefile 的职责。')
section('一个 make 编译所需 CEE 组件',
 code('cd cee-kfparticle\nmake -j4'),
 table([['生成物','用途'],['CEErco/CeeV0Lite/libCeeV0Lite.a','独立核心静态库'],['CEErco/CeeV0Analysis/CeeV0Analysis.bin','既有 CLI 分析程序'],['CEErco/CeeV0Analysis/libCeeKFParticleAnalysis.so','供 ROOT macro 调用的共享库']],[290,CW-290]),
 P('在 CEErco 目录也可以执行 make。顶层入口传播并行编译参数，无需手动逐目录调用。root-config、VMCWORKDIR、FAIRROOTPATH 应来自你的实验环境。StRoot 仅作参考，不属于此构建目标。'),
 code('make lite     # standalone core\nmake check    # core + analysis tests\nmake legacy   # optional old Analyser.bin\nmake clean'),
 P('默认 make 会检查 VMCWORKDIR/CeeBenoDst/edm/CeeBenoEvent.h；存在时自动编译 BenoSource 并链接 libCeeBenoDst。可用 make WITH_BENODST=1 显式启用，0 禁用。切换该选项后先 make clean。仅使用原始输入时，不要求安装 BenoDst。'),
 P('旧 Analyser.bin 不是新 macro 的运行依赖；需要时用 legacy 目标。旧 Makefile 中失效的 AnaMaker/ 路径已修正。完整实验库链接仍需在实际 CeeRoot/FairRoot 环境验证。'))
section('用同一个 macro 选择输入和动量',
 code('cd CEErco/CeeV0Analysis\nroot -l -b -q \\\n \'macros/runCeeKFParticle.C("input.list","mothers.root",\n "benodst","tof",0.,-1,true)\''),
 P('调用中的换行仅为排版，命令实际写在一行。参数依次为 input、output、inputFormat、momentumMode、ByTesla、maximumEvents、realData、overwrite。上例为 BenoDst 零场实验数据；时间分辨率示例需要按实际标定修改。'),
 table([['场景','inputFormat / momentumMode / ByTesla'],['原始 track+hit 零场','"raw" / "tof" / 0.'],['BenoDst 零场','"benodst" / "tof" / 0.'],['原始有磁 track','"raw" / "tracking" / 0.5'],['BenoDst 有磁 track','"benodst" / "tracking" / 0.5']],[175,CW-175]),
 P('By 是沿 y 的带符号磁场，单位 T。tof 与非零场组合会报错；tracking 可以明确使用 By=0，但输入必须确有可信动量。auto 仅按场强选择兼容模式。旧的第三参数为数值 By 的调用仍然可用。input 支持 ROOT、list/txt 与目录；输出默认不覆盖。'))
section('通用通道 macro',
 code('root -l -b -q \\\n \'macros/runCeeChannels.C("input.list","omega.root",0.5,"omega",true)\''),
 P('末尾 true 显式开启仅按电荷组合的流程演示，不是已做标定 PID 的物理结果。正式运行用 false 并修改 ChannelPID。支持示例 lambda、omega、hypertriton2、hypertriton3。He³ 需要 q=+2 与物理动量，程序不猜测刚度转换。'))

chapter('输入层：原始 track/hit 与 BenoDst','先选择输入格式，再选择动量来源；两种输入最终都形成同一个 AnalysisEvent。')
section('两条路径在哪里汇合',
 table([['步骤','raw','benodst'],['读取','CeeTpcTrack、PV、TOF hit','Event、TpcTrack、CovMatrix、MatchInfo'],['TOF 关联','按投影距离进行唯一匹配','按已保存的 best partner 查 MatchInfo'],['位置与时间','由原始 hit 和配置标定得到','直接用 MatchInfo 的已对齐位置、已校准时间'],['统一输入','AnalysisEvent / TrackInput','AnalysisEvent / TrackInput'],['重建和输出','Maker → Finder → Writer','相同 Maker → Finder → Writer']],[90,200,CW-290]),
 P('BenoDst 已提前完成 track 与 TOF 的匹配。KFP 不再调用几何匹配器，也不重复套用 TOF 时间偏移。新的 CeeKFParticleBenoSource 用 ROOT TChain 直接读取单文件或列表；无需额外建立 FairRootManager。已有 FairTask 可继续使用自己的 adapter。'))
section('先选 track，再适配',
 code('CeeKFParticleBenoSource source(config,benoOptions);\nsource.SetTrackIndexSelector(\n [](long long entry,const TClonesArray& tracks) {\n   return MySelectedBenoIndices(entry,tracks);\n });'),
 P('回调由用户实现，返回 Beno TpcTrack 数组下标；空列表表示不选任何 track，重复或越界索引报错。Source 每事件索引匹配记录一次，只适配选中的 track；没有 track×hit 的几何搜索。上游 trackId 则保存在 TrackInput.id，用于输出及 daughter 关联。'))
for i in (25,26,27,28,30):
 title,items,src,_=OLD[i-1];section('raw：'+title,*items,src=src)
section('BenoDst 分支和字段映射',
 P('源码核对版本：CEESM/CeeRoot master 83e58c2f523dcfc1644068afe46c347280e9db80。默认 treeName="ceesim"；实际树名不同就在 config.treeName 修改。读取前会检查列表中每个文件的必需分支。上游跳过的空 Event 条目会跳过，eventId 仍为链中的条目号。'),
 table([['数据','用途与约定'],['Event','primaryVertexIndex 有效；PV 位置及 sigma，sigma 平方为协方差'],['TpcTrack / TpcTrackCovMatrix','track 身份、几何/动量、5D 误差；核对双向协方差索引'],['TpcItofMatchInfo + ItofHit','iTOF 匹配记录及 partner 索引边界'],['TpcEtofMatchInfo + EtofHit','eTOF 匹配记录及 partner 索引边界'],['*BestMatchIndex()','指向 hit 数组，不是 MatchInfo 数组'],['trackId()','上游身份，可与 Beno 数组下标不同'],['nSigma* / status / trackChi2','PID 缺失转 NaN；status 是顶点关联，χ² 没有 NDF']],[180,CW-180]),
 P('零场至少有一组 TOF MatchInfo/hit 分支；tracking 不要求 TOF。PV 缺失或误差无效跳过事件；协方差无效或无可用匹配跳过该 track；坏索引、重复输入 ID 和超出数量限制会明确报错。Source 的统计给出 invalidGeometry、missingMatch、sharedMatch，避免只看到空质量谱。'))
section('零场几何：保留前四维及相关误差',
 P('原 CeeTpcTrack 本身也有 FitParameter[5]/Fitcov[5][5]。BenoDst 保存其 a=(d,z0,φ,θ,q/p)。前四维描述直线几何；无磁时 q/p 可为占位，不能拿它恢复动量。上游 fillPar5AtVertex 以过原点的 beam z 轴为参考，d=(x py-y px)/pt。'),
 eq(r'\mathbf r_0=(d\sin\phi,-d\cos\phi,z_0),\quad \mathbf u=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)',13),
 P('零场用前四维构造 g=(r0,u) 及 Cg=J6×4 C4 Jᵀ。Cg 可只有秩 4；这是直线几何的正常约束，不补一个虚假的 q/p 误差。位置、方向及二者的交叉项都保留在 geometryCovariance 中。原始适配器未提供此矩阵时，原来的位置协方差加标量角误差模型保持不变。'),
 P('每次 TOF 给出动量大小 p 和其误差后，用 A=diag(I,pI) 将 Cg 变成位置—动量协方差，再加入沿 u 的 TOF 动量误差。若方向朝向 TOF hit 时翻转，位置—方向交叉项也随之翻转。这仍沿用现有 TOF 路程/时间误差模型，不声称已做所有路程—几何相关项的联合拟合。'),
 P('会将 5D 参数重建出的几何与 origin/globMom 的横向垂足和方向对比。默认允许 0.10 cm 与 0.003 rad 的量化差异；允许方向相反，因为上游 GetEntryMomentum 会反向保存 globMom。零场最终由 TOF hit 定向；其他不一致拒绝。此检查可发现部分上游退回原始 GenFit 参数却未改标签的情况；无法替代对输入生产版本的核对。'),
 src='CEETPC_Tracker/CeeTrackingTask.cxx: fillPar5AtVertex；CeeKFParticleBenoSource.cxx: PerigeeState；CeeV0LiteTypes.cxx: BuildTofStateCovariance。')
section('tracking 输入如何使用第五维',
 P('tracking 使用完整 C5，并恢复上游 GenFit 固定局部面上的状态误差后计算 C6=J6×5 C5 Jᵀ；不能只旋转到新垂足而不更新动量。拟合方向来自 φ/θ，与原 VertexSeed 一致，不直接采用反向保存的 globMom。默认转换限定单位电荷，p=1/|q/p|，要求 q/p 误差有效、符号与 charge 一致，且动量与 globMom 大小相容。随后由 Finder 按默认配置输运到 track–PV DCA 参考点；Beno 打包本身并没有新增这一步。'),
 P('He³ 等 |q|>1 核素需要用户用 SetStateConverter 明确提供物理电荷、物理动量与相应协方差。不能从 Beno charge() 的 ±1 独自推断电荷大小，也不自动把刚度当动量。其他上游坐标/参数定义同样应使用此回调。'))
section('已有 TOF 匹配如何使用和筛选',
 P('按 track 数组索引和保存的 best hit partner 查找匹配记录，不把 BestMatchIndex 当成 MatchInfo 下标。默认接受 flag=1 或 2 的 best 记录；requireUniqueMatch=true 时只接受 flag=1。默认 rejectSharedHits=true，拒绝标记为共享或被多个 track 声明使用的 hit。不会给被拒绝的 track 重新寻找另一个几何匹配。'),
 code('BenoInputOptions options;\noptions.rejectSharedHits=true;\noptions.requireUniqueMatch=false;\noptions.preferItof=true;\noptions.maximumMatchChi2=1.e30;'),
 P('若 iTOF 与 eTOF 都有可用 best 记录，preferItof 决定优先来源；默认 iTOF，必要时退回 eTOF。这里只选择一次测量，不把同一 track 的两次 TOF 当作独立信息合并。位置取 MatchInfo 已对齐的 hitPos，hit 分支只核对 partner 索引。'))
section('时间原点已经扣 T0，误差不能漏掉',
 P('两个上游 matcher 均保存 TOF=(hitTime-offset)-T0。因此 Beno 模式设置事件时间原点 startTimeNs=0，把保存的 itof()/etof() 作为 tofTimeNs；不读原始 T0Hit，不再减 T0。triggerTs 是 Unix 时间戳，不能用作此时间原点。'),
 eq(r'\Delta t_i=t_{\mathrm{Beno},i}-t_{\mathrm{parent}},\qquad t_{\mathrm{parent}}=L_{\mathrm{PV}\to\mathrm{SV}}/(\beta_M c)',15),
 P('这里 L 是 PV→SV 在母粒子动量方向上的有符号投影。仍须扣除母粒子飞行时间，随后按 SV→TOF 路程计算子粒子速度并迭代。Beno 中按 primary 假设给出的 beta 不作为最终 daughter 动量输入。'),
 table([['配置','含义'],['config.t0.source / fixedTimeNs','kFixedT0 / 0：相对已选 T0 的坐标原点'],['config.t0.sigmaNs','上游所用 T0 的共同测量误差；实验数据由标定提供'],['config.tofTimeSigmaNs','单个 TOF hit 的独立时间误差，不再包含 T0 项'],['config.tofPathLengthSigmaCm','现有 TOF 动量模型使用的路程误差'],['propagateSharedTimeUncertainty','共同 T0 误差只传播一次到母粒子'],['tofRawToNsScale / tofOffsetNs','Beno 固定 1 / 0，防止重复标定']],[218,CW-218]),
 P('ConfigureInputModes 会处理 Beno 的零时间原点，但不会把真实 T0 的误差清零。真实数据 macro 中的 0.10 ns 等数值仅为配置示例。理想模拟仍使用精确起点和零 T0 误差；若模拟包含 T0 测量抖动，应显式按测量时间契约配置。'),
 src='CeeTpcItofMatcherTask.cxx、CeeTpcEtofMatcherTask.cxx 的 TOF 计算；CeeKFParticleInput.h。')

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
 P('两体完整示例见 macros/runCeeKFParticle.C，通用通道见 runCeeChannels.C；后者末尾可传 inputFormat="benodst"，仍仅支持 tracking。它们通过 OpenCeeInput 选择 Source。在已有 Beno FairTask 循环中，用 adapter 构造 AnalysisEvent，然后 MakeChannels 与 writer.Write。Begin/Finish 放在初始化与结束，不要每事件重建输出文件。'))
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
section('统一输入中的几何、动量和时间',
 table([['字段','含义与维数'],['positionCm / direction','同参考位置的 track 位置与方向；零场不使用输入动量模长'],['positionCovarianceCm2 / sigmaDirectionRad','旧输入的 3×3 位置误差与标量方向误差'],['hasGeometryCovariance / geometryCovariance','可选 (x,y,z,ux,uy,uz) 的 6×6 矩阵；Beno 零场由 4D 几何构造'],['hasTrackingMomentum / trackingMomentumGeV','是否提供物理动量及其三分量'],['hasTrackingCovariance / trackingCovariance','tracking 状态 (r,p) 的 6×6 协方差'],['tofPositionCm / tofTimeNs','选用 TOF 的位置与时间；Beno 使用已扣 T0 的 ns'],['startTimeNs / sigmaStartTimeNs','事件时间原点与共同 T0 误差；Beno 原点为 0']],[235,CW-235]),
 P('geometryCovariance 与 trackingCovariance 含义不同，不能互相复制。旧代码不设置 hasGeometryCovariance 时保持原来的 TOF 误差构造。数据结构变化后应重新编译核心和共享库；顶层 make 会完成所需目标。'))
section('已有两体的 89 个标量目录',
 P('以下表属于原两体 RootWriter。固定 [0]/[1] 仅约束旧 Candidate，不限制新 Channels。DcaProton/Pion 是按 PDG 查找的便利函数；通用角色可取 daughter1、daughter2、he3、deuteron 等任意唯一名字。'))
reuse(*range(14,25))
