page('从一个事件到一行候选：实际运行顺序',
    P('理解输出前，先沿程序真正执行的顺序走一遍。Analyzer 负责事件层的输入检查与匹配；Finder 负责track对的求解和选择。只有 Finder 最终返回的候选，才会进入输出树。'),
    table([['事件步骤','继续执行的条件'],['读取输入事件','ROOT 成功读到该事件'],['选择重建 PV','找到坐标和所用误差有效的顶点'],['确定时间原点','模拟 t₀=0；真实数据有合规时间'],['选择 TPC track','至少两条通过质量与参考状态检查'],['准备动量输入','零场至少两个唯一 TOF 匹配'],['形成核心输入并调用 Finder','至少两条通过 adapter 的track'],['写 Λ 候选及 QA','track对求解成功并通过当前选择']],[225,274]),
    P('中心拟合完成后，候选先经过质量、两track DCA、几何拟合、衰变长度和指向等中心结果的筛选，再做代价较高的完整母粒子协方差计算。这样可以减少对明显不合格组合的重复求解。'),
    P('hEventCounter 按上述过程记录事件通过了哪一步。其中 finder run 表示至少调用了一次 Finder，has Lambda 才表示事件产生了被接受的 Λ。它不能说明每一对失败track的具体原因。'),
    note('真值分母的位置','生成 Λ 的真值 QA 在 PV、时间和track选择之前填入。因此做效率研究时，需要明确分母是否还要求几何接受度、可重建子track等条件，不能仅凭两张直方图相除就认为得到了全部效率。'),
    src='[A] Analyzer::Run、ProcessEvents、FillGeneratedLambdaTruth；[F] BuildCandidate')

page('先跑通普通模拟：编译与最小运行示例',
    P('完整分析程序需要已经配置好的 CEE、ROOT 和 FairRoot 环境。进入 CEErco/CeeV0Analysis 后编译；Makefile 通过 root-config 查找编译选项，并使用 VMCWORKDIR 与 FAIRROOTPATH 查找实验数据类和动态库。'),
    code('cd /path/to/CEErco/CeeV0Analysis\nmake\n./CeeV0Analysis.bin \\\n  --input /path/to/reco.root \\\n  --output CeeV0Lambda.root \\\n  --mode simulation --momentum tof --by-tesla 0 \\\n  --cut-preset legacy --events 10000 --debug 1'),
    P('示例采用当前零场默认路径，处理至多一万个事件。输入可以是单个 ROOT 文件、文件列表或目录；默认输入树名为 ceesim。若文件中的树使用其他名称，需要显式给 --tree。输出使用 RECREATE，同名文件会被本次结果覆盖。'),
    sub('先看是否真的进入了重建'),
    P('第一次运行时，优先检查事件计数、有效 PV、参考状态来源、唯一 TOF 匹配数，以及最终是否出现候选。随后再看质量峰和拓扑分布。若仍卡在输入阶段，直接收紧物理 cut 并不能定位问题。'),
    code("root -l -b -q \\\n  'macros/plotCeeV0.C(\"CeeV0Lambda.root\",\"CeeV0QA.pdf\")'"),
    note('运行宏与二进制使用同一实现','runCeeV0Analysis.C 只是组装参数后调用 CeeV0Analysis.bin。需要完整控制时间、误差和选择时，直接使用命令行更清楚；它没有另一套独立的重建算法。'),
    src='[C] main；CeeV0Analysis/Makefile；[A] ResolveInputFiles、SetupOutput；macros/runCeeV0Analysis.C、plotCeeV0.C')

page('真实数据先确认时间：单位、偏移和分辨率分别设置',
    P('真实数据默认从 CeeT0Hit 取得事件起始时刻。程序要求显式输入 --t0-scale-to-ns，目的是让使用者确认原始时间单位。若原始数值已经是 ns，填 1；若是 ps，才应填 0.001。下例假定已确认原始 T0 为 ns。'),
    code('./CeeV0Analysis.bin \\\n  --input /path/to/real.list --output CeeV0Lambda_real.root \\\n  --mode real --momentum tof --by-tesla 0 \\\n  --t0-source hit --t0-scale-to-ns 1 \\\n  --cut-preset legacy --events 10000 --debug 1'),
    P('时间转换分为两个独立过程：T0 原始值使用 T0 的 scale 与 offset；TOF 原始值使用 TOF 的 scale 与 offset。应先查看 hT0Raw、hT0 和 hTofMinusT0，再确认时间差是否符合事件时标。改变单位与改变分辨率不是同一个操作。'),
    eq(r't_{\rm ns}=t_{\rm raw}\,a+b',size=18),
    P('一个事件有多个有效 T0 hit 时，默认取中位数，也可选 first、earliest 或 mean。配置的 T0 sigma 表示所用事件时间的不确定度，程序不会因 hit 数增多自动除以平方根。'),
    note('固定起始时间仍要有依据','真实数据可用 --t0-source fixed 与 --fixed-t0-ns，但不会因此自动消除公共时间误差。real 默认的 T0 sigma 和额外公共时间 sigma 均为 0.10 ns，它们是占位配置，必须按实际标定结果核对，避免把同一误差算两次。'),
    src='[C] ConfigureDataMode、ValidateDataModeTiming、main；[G] SelectEventT0；[A] EventT0、ProcessEvents')

page('磁场路径怎样使用：输入要求比零场更严格',
    P('v0.6.5 默认先做 track–PV DCA 参考点转换。--magnetic-track-reference stored 仅用于重放旧输入参考约定；--reference-max-path-cm 设置搜索半范围，默认 1000 cm。'),
    P('CEE 磁场沿 y，束流沿 z。有场运行使用 --momentum tracking 和 --by-tesla，数值单位为 T，正负号表示 ±y。程序按 |By| 判断零场：tof 要求 |By|≤10⁻⁶ T，tracking 要求 |By|&gt;10⁻⁶ T。'),
    table([['要求','原因'],['完整 Cartesian 六维状态与协方差','位置、动量和误差必须属于同一参考点'],['有效拟合电荷','螺旋输运需要已知电荷正负'],['不接受旧状态误差回退','方向近似不能替代曲率动量及相关误差'],['当前不绑定 TOF hit','尚未用螺旋 matcher 验证 TOF 对应关系']],[240,259]),
    P('命令行切换到 tracking 时，会关闭默认的零场 DCA 角色分配。Λ 仍按质子与 pion 假设组合，但有效拟合电荷可以排除不相容的正负号。适配层要求拟合电荷有限，且其绝对值至少达到配置阈值 0.25。'),
    P('真实 tracking 模式仍走通用事件时间选择。若使用 --mode real 而不改 T0 来源，默认仍要求 CeeT0Hit；这不意味着该 hit 在替代曲率决定动量。手工固定时间也应有明确的数据约定。'),
    note('运行参数与输入检查','旧 --bz-tesla 仅接受 0 并提示迁移，非零值报错；C++ 字段改为 magneticFieldByTesla，需重编译。当前仓库未见上游生产者实际填写种子与 PV 误差的调用。有场分析仍须先检查完整跟踪状态，再做物理闭合。'),
    src='[C] main；[D] MakeTrackInput、SelectTrackReferenceState；[A] SetupInput、ProcessEvents；docs/CEE_INPUT_PRODUCER.md')

page('三组预设怎样理解：先保留样本，再验证选择',
    P('legacy 便于与旧 CEE 分析做第一轮对照；loose 使用较宽的几何选择；star 借鉴了一组 STAR 数值。它们改变的是选择条件，不会自动把 CEE 的误差模型和实验性能变成 STAR 的实现。'),
    table([['条件','legacy','loose','star'],['质量 / GeV/c²','0.7–2.0','1.05–1.25','1.08–1.18'],['两track DCA 上限 / cm','放开','3','1'],['几何 χ²/NDF 上限','放开','25','10'],['总拓扑 χ²/NDF 上限','放开','25','5'],['衰变长度下限 / cm','0','0.2','5'],['L/σL 下限','0','0','5'],['指向余弦下限','−1','0.90','−1'],['子track到 PV 的原始 χ² 下限','0','0','18.6'],['完整母粒子 Jacobian','要求','要求','要求']],[217,94,94,94]),
    P('表中的“放开”在代码中用 10³⁰ 表示。表内上下限通过条件含等号；DCA 角色中 pion 一侧另要求严格大于分界。即使选择 legacy，仍要求输入和中心解有效、质量和衰变长度为正、默认 TOF 收敛、默认 DCA 角色规则满足，以及完整母粒子 Jacobian 成功。它不能理解为完全不做筛选。'),
    P('总拓扑量使用几何项与 PV 项之和除以总自由度。star 预设中的 18.6 则用于子track到 PV 的原始 χ²，不是 sqrt(χ²)。若把这个数填到旧 deviation 选项，会形成另一套选择。'),
    note('比较样本时记录实际配置','程序先应用模式和预设，再处理逐项覆盖。输出中的 CeeV0AnalysisConfig 保存实际参数。一次运行可用 CLI 覆盖；永久调整预设则修改 CeeV0AnalysisConfig.cxx 的 ApplyCutPreset 并重新编译。'),
    src='[C] ApplyCutPreset、main；[F] ChannelCuts、PassCuts、BuildCandidate')

page('命令行速查：输入、模式与处理规模',
    P('以下选项控制文件和事件层面的行为。路径与树名需要和数据实际结构一致；用较少事件定位输入问题后，再扩大处理规模。布尔开关只写选项名，其余选项后接相应值。'),
    table([['选项','含义与默认'],['--input','必填：.root、.list/.txt 或目录'],['--output','输出路径；默认 CeeV0Lambda.root'],['--tree','输入树名；默认 ceesim'],['--mode','simulation 或 real；默认 simulation'],['--momentum','tof 或 tracking；默认 tof'],['--by-tesla','均匀 By / T；默认 0'],['--cut-preset','legacy、loose、star；默认 legacy'],['--events','最多事件数；−1 表示全部'],['--files','最多输入文件数；−1 表示全部'],['--debug','0、1、2；默认 0'],['--no-production-state','关闭新增八维计算；默认开启'],['--no-truth','模拟模式不绑定 MCTrack'],['--help 或 -h','打印可执行程序自身的帮助']],[214,285]),
    P('debug=1 用于观察部分事件的输入数、可重建track数、唯一 TOF 数和候选数；debug=2 还会报告部分track适配失败。它没有为所有被拒的track对自动保存一行详细诊断。'),
    note('没有候选不等于文件读取失败','程序可以正常读完文件并输出一个空候选树。应结合 hEventCounter 和状态直方图判断停在哪一步，再决定检查输入、时间、匹配还是物理选择。'),
    src='[C] Usage、main；[A] ResolveInputFiles、ProcessEvents')

page('命令行速查：时间与 TOF 匹配',
    P('时间 scale、offset 和 sigma 分别表示单位转换、零点校准和测量不确定度。三者不能互相替代。下表中的真实数据时间误差默认值是当前配置占位，不代表已经测定的探测器性能。'),
    table([['选项','含义与默认'],['--t0-source','真实数据：hit 或 fixed；默认 hit'],['--fixed-t0-ns','固定起始时刻 / ns；默认 0'],['--t0-scale-to-ns','实测 T0 单位转换；real+hit 必须显式给'],['--t0-offset-ns','T0 转换后偏移 / ns；默认 0'],['--t0-sigma-ns','T0 sigma / ns；模拟 0，真实默认 0.10'],['--t0-policy','first、earliest、mean、median；默认 median'],['--tof-time-scale-to-ns','TOF 原始时间转换系数；默认 1'],['--tof-time-offset-ns','TOF 转换后偏移 / ns；默认 0'],['--tof-time-sigma-ns','单条 TOF 时间 sigma / ns；默认 0.10'],['--tof-path-sigma-cm','路程 sigma / cm；默认 0.5'],['--additional-decay-time-sigma-ns','额外公共时间 sigma；模拟 0，真实默认 0.10 ns'],['--max-beta','允许的 beta 上界；默认 0.9995'],['--itof-match-cm','iTOF 匹配距离上限；默认 10 cm'],['--etof-match-cm','eTOF 匹配距离上限；默认 10 cm']],[252,247]),
    P('固定时间源下，T0 的 scale、offset 和多 hit 策略没有作用，CLI 会拒绝这种混用。普通 simulation 还要求固定 t₀=0，并关闭公共时间误差；单条 TOF 的测量误差仍然保留。'),
    src='[C] Usage、ConfigureDataMode、ValidateDataModeTiming、main')

page('命令行速查：track质量、角色与回退误差',
    P('这组参数决定哪些输入track能进入重建，以及信息缺失时允许采用什么误差假设。回退值应通过模拟残差和 pull 分布核对，不能仅为了让更多track被接受而任意缩小。'),
    table([['选项','含义与默认'],['--min-tpc-hits','拟合 hit 数下限；默认 0'],['--max-track-chi2ndf','track χ²/NDF 上限；默认实质放开'],['--max-tracks','每事件最多选取track数；默认 1000'],['--max-candidates','每事件最多接受候选数；默认 100000'],['--require-track-converged','要求存储的track拟合收敛标志'],['--strict-pv-errors','不允许缺失或无效的 PV 误差回退'],['--disable-legacy-track-fallback','零场也必须有有效 vertex seed'],['--disable-dca-role','关闭质子/pion 的 DCA 角色规则'],['--lambda-role-dca-cm','角色分界；默认 1.2 cm'],['--fallback-track-sigma-x-cm','回退track位置 sigmaₓ；默认 0.10 cm'],['--fallback-track-sigma-y-cm','回退track位置 sigmaᵧ；默认 0.10 cm'],['--fallback-track-sigma-z-cm','回退track位置 sigma_z；默认 0.15 cm'],['--fallback-direction-sigma-rad','回退方向 sigma；默认 0.003 rad'],['--fallback-pv-sigma-x-cm','回退 PV sigmaₓ；默认 0.10 cm'],['--fallback-pv-sigma-y-cm','回退 PV sigmaᵧ；默认 0.10 cm'],['--fallback-pv-sigma-z-cm','回退 PV sigma_z；默认 0.15 cm']],[259,240]),
    P('达到track或候选上限时，程序会截断后续处理。若真实事件经常触及上限，应先核对其比例和对样本的影响，不能把上限视为没有物理影响的参数。'),
    src='[C] Usage、main；[D] CeeAdapterConfig；[A] SelectTracks；[F] Finder::Reconstruct')

page('命令行速查：逐项选择与诊断开关',
    P('使用者可以在预设基础上覆盖下表条件。质量、距离和统计量的量纲不同；尤其要区分子track到 PV 的原始 χ²与它的平方根。下限和上限的方向也应与选项名称一起核对。'),
    table([['选项','所控制的量'],['--mass-min / --mass-max','候选不变质量上下限 / GeV/c²'],['--max-daughter-dca','两条子track DCA 上限 / cm'],['--max-geo-chi2ndf','共同顶点几何 χ²/NDF 上限'],['--max-topo-chi2ndf','几何项加 PV 项后的总 χ²/NDF 上限'],['--min-decay-length','|SV−PV| 的下限 / cm'],['--min-ldl','L/σL 的下限'],['--min-pointing-cos','母粒子动量与 PV→SV 的指向余弦下限'],['--min-daughter-pv-chi2','两条子track各自到 PV 的原始 χ²下限'],['--min-daughter-pv-deviation','旧兼容选项：sqrt(χ²) 下限'],['--allow-nonconverged-tof','诊断：允许可用但未达到收敛阈值的迭代'],['--allow-incomplete-mother-jacobian','诊断：只保留特定外层扰动失败的中心解']],[267,232]),
    P('最后两个开关分别放宽不同条件。它们不会让中心 POCA、共同顶点或无效速度自动变成有效解；也不应作为“增加统计量”的普通物理配置。诊断结果中有些误差和拓扑量会是 NaN，必须结合状态标志解释。'),
    note('当前没有对应 CLI 的内部参数','TOF 最大迭代次数、顶点与时间收敛阈值、Jacobian 步长等在配置结构中设置并写入运行元数据，但没有同名命令行选项。不要照着成员名自行拼一个 --参数。'),
    src='[C] Usage、main；[F] BuildCandidate、PassCuts、FinderConfig')

page('ROOT 输出怎样读：一行保存的是被接受的 Λ',
    P('当前候选树名为 LambdaCandidates。每行保存一个通过当前配置的 Λ 假设，daughter 0 是质子假设，daughter 1 是 pion 假设。一个事件可以有多行；没有候选的事件只反映在事件与 QA 统计中，不会自动生成一行空记录。'),
    table([['字段组','读图时应回答的问题'],['event、trackId0/1','候选来自哪个事件、哪两条重建track？'],['pvX/Y/Z、svX/Y/Z','使用的主顶点和重建衰变顶点在哪里？'],['px/py/pz、energyGeV、massGeV','SV 母粒子的动量、能量和不变质量是多少？'],['geometricChi2 与 NDF','两条track形成共同顶点是否相容？'],['primaryVertexConstraintChi2','母track相对 PV 的单独统计项有多大？'],['topologicalChi2 与 NDF','加入 PV 项后，总统计量是多少？'],['decayLengthCm、pointingCosine','两个顶点相距多远，母粒子动量是否指回 PV？'],['decayTimeNs、tofConverged','TOF 耦合求解给出什么衰变时刻，是否收敛？']],[231,268]),
    P('decayLengthCm 是几何距离 |SV−PV|，不是沿动量方向的有符号投影；decayTimeNs 是事件时标上的衰变时刻，不是固有寿命 cτ。v0.6.3 的输出格式为 6，另有 motherTopo* 分支，保存八维约束副本、S 和 cτ；字段含义见八维状态一章。'),
    note('真值标签在重建之后附加','truthAvailable=1 才表示两条track的 MC 关系可解析；truthMatched=1 还要求来自同一个真 Λ，且质子/pion 角色正确。缺少真值不能直接归入已知组合背景。'),
    src='[A] CreateCandidateBranches、FillRow、MatchTruth；[F] BuildCandidate')

page('输出协方差怎样还原：七维状态与 28 个元素',
    P('原有输出 motherCovariance 不是一张独立的新测量表，而是 SV 母粒子七维状态的不确定度。其行列依次对应位置、动量和能量。这一分支仍为 28 项；新增的 motherTopoCovariance 才是 36 项，二者不能混读。'),
    eq(r'\mathbf y=(x,y,z,p_x,p_y,p_z,E)^T,\qquad C_{ij}=\operatorname{Cov}(y_i,y_j)',size=17),
    P('由于协方差对称，树中只保存逐行下三角，共 7×8/2=28 个数。索引从 0 开始；第 i 行依次保存第 0 到 i 列。'),
    eq(r'k=\frac{i(i+1)}{2}+j,\qquad 0\leq j\leq i\leq6',size=18),
    code('int k = 0;\ndouble C[7][7] = {};\nfor (int i = 0; i < 7; ++i) {\n  for (int j = 0; j <= i; ++j) {\n    C[i][j] = motherCovariance[k++];\n    C[j][i] = C[i][j];\n  }\n}'),
    P('读取前先检查 motherCovarianceValid。位置块给出 SV 误差；动量与能量块共同给出质量误差；位置—动量相关则影响母粒子相对 PV 的统计检验。保留非对角项，是为了把这些来自同一次拟合的相关变化一起传播。'),
    note('单位随行列变化','位置—位置项用 cm²；位置—动量项用 cm·GeV/c；动量—动量项用 (GeV/c)²。不能把整个矩阵都当作“位置误差”，也不能把某个协方差元素直接当标准差。'),
    src='[A] FillRow、SetupOutput；[Y] MotherState、MotherMassAndError')

page('0.6.1 的诊断保留：哪些失败可以留下来',
    P('默认物理配置要求中心解可用，并要求完整母粒子 Jacobian 成功。外层 Jacobian 对两条六维输入共 12 个分量做正负扰动，因此通常有 24 次扰动重拟合。某一侧失败时，不能再把算不完整的协方差当成有效结果。'),
    P('--allow-incomplete-mother-jacobian 只保留一类情况：中心 POCA、顶点和可用中心解已经获得，但外层扰动重拟合失败。质量、SV、母粒子动量、DCA 和几何 χ²等仍可用于定位问题；最终母粒子协方差及依赖它的量被明确置为 NaN。'),
    table([['诊断字段','解释'],['motherJacobianComplete','外层母粒子 Jacobian 是否完整'],['motherJacobianPerturbationFailed','是否属于可保留的扰动失败类型'],['motherJacobianFitsAttempted / Succeeded','尝试数与成功数；普通完整解为 24/24'],['motherJacobianFailureCount','失败的正负侧总数；不是失败列数'],['motherJacobianPlusFailureMask / MinusFailureMask','12 列的正侧、负侧失败位图'],['motherJacobianFailedColumn / FailedSide','第一次失败的列；正侧 +1、负侧 −1'],['motherJacobianFirstFailureFitStatus','第一次失败时的拟合状态'],['covarianceDependentCutsEvaluated','是否真正计算并应用最终协方差相关 cut']],[285,214]),
    P('保留行的质量误差、SV 误差、σL、L/σL、PV 统计项和总拓扑项无效；NaN 不可改成 0。几何 χ²与子track到 PV 的 χ²仍使用输入协方差，所以也不能把它们称作“完全不依赖误差”。'),
    note('从诊断文件恢复严格样本','通常先要求 motherJacobianComplete==1。若同一文件还用了 --allow-nonconverged-tof，则须同时要求 tofConverged==1；仅凭 Jacobian 标志不能恢复全部默认条件。'),
    src='[F] BuildCandidate、InvalidateCovarianceDependentObservables、PassCuts；[A] CreateCandidateBranches；MIGRATION_0_6_1.md')

page('先做哪些 QA：让每张图回答一个明确问题',
    P('第一次检查结果时，建议沿输入、求解、物理量的顺序读图。只看到 Λ 质量峰，能够说明某些组合形成了合理中心值，却不能证明位置误差、质量误差和所有拓扑 cut 已经标定正确。'),
    table([['查看对象','首先要确认的事情'],['hEventCounter','事件主要停在哪一阶段？'],['hTrackReferenceSource','多少track用了真实 seed，多少用了回退？'],['hReferenceAdapterStatus / hTrackAdapterStatus','参考状态或最终核心输入为何被拒绝？'],['hPvStoredErrors','PV 误差是否来自文件，还是回退配置？'],['hT0Raw、hT0、hTofMinusT0','时间单位、零点和有效时间差是否一致？'],['hTofMatches、hTofMatchDistance','匹配数量与残差是否合理？'],['hLambdaMass、真值与背景分量','质量峰和组合背景分别发生了什么变化？'],['几何、PV-only、总拓扑 χ²','是哪一项改变了候选相容性？']],[260,239]),
    P('若要验证误差，可在有匹配真值的样本上比较重建值与真值的残差，并研究残差除以重建标准差后的 pull。这样的检查应同时覆盖不同动量、衰变位置和输入状态来源。对近似高斯、无偏且未受强选择扭曲的验证样本，pull 的均值应接近 0、宽度应接近 1；偏离时还要排查错配尾与选择效应。'),
    P('还可将 seedParentBeta 改为 0.5、0.7、0.9，比较共同收敛样本的最终状态和各自收敛率。内置 plotCeeV0.C 提供质量及若干基础分布。开启诊断保留后，要另外按有效性标志分组查看；不能默认每张内置直方图已经替你恢复了严格物理样本。'),
    note('知道这些图没有记录什么','事件计数不是候选级完整 cutflow；adapter 状态图也不包含所有早期track质量拒绝。若要量化各项选择效率，需为明确的分母和失败类别另做统计。'),
    src='[A] SetupOutput、SelectTracks、ProcessEvents；macros/plotCeeV0.C')

page('数据类安全修复的边界：初始化不等于真实测量',
    P('0.6.1 配套代码把许多未填字段改为确定的初始值，使数据对象不再依赖未初始化内存。CeeTpcTrack 的构造函数初始化拟合参数、协方差数组和 seed 数组；未填 MC 关系为 −1，未保存 seed 时标志为 false。'),
    P('这些改动提高了数据类的可解释性，但不会替上游补做测量。比如默认拟合收敛标志为 true，不能由这个默认值反推拟合确实成功；是否可信，仍取决于生产者是否填写。dE/dx、NDF 等字段同样需要有明确来源。'),
    sub('保存了一个 seed，还要检查它是否完整有效'),
    P('SetVertexSeed 会复制给定的位置、动量和矩阵；若输入矩阵不够六维，缺项补 0，然后设置存在标志。setter 本身没有保证物理有效性，适配层随后仍需检查数值和协方差。另一个 5×5 track参数协方差使用 DCA、角度和 q/p，不能直接当成六维 Cartesian 矩阵。'),
    P('CeeTpcVertex 的 Clear 会清除坐标、误差、rank 和track DCA 向量。当前数据类版本为 Track 8、Vertex 2；这说明当前 schema 的版本，并不证明某个旧 ROOT 文件已保存全部新增信息。'),
    note('更新代码后还必须更新生产链','修改持久数据类后，应重新编译对应 CEE 库与 ROOT 字典，并用实际文件验证读写。当前仓库未显示生产者填写完整 seed 与 PV 误差的调用；普通零场可按文档使用回退，磁场与精确误差研究则需要真正写入这些测量信息。'),
    src='CeeTpcTrack.cxx：构造函数、SetVertexSeed、GetTrackParameterCovariance；CeeTpcVertex.h：Clear；[D] CopyVertexSeed；docs/CEE_INPUT_PRODUCER.md')

page('失败状态速查：先看是哪一层没有完成',
    P('motherCovarianceStatus 使用下列整数编码。除有效结果和显式保留的扰动失败 QA 外，其他失败通常已经被 Finder 拒绝，不能假设在输出树中一定能找到相应行。调试时应结合函数返回与诊断信息定位。'),
    table([['值','MotherCovarianceStatus 的含义'],['0','尚未计算母粒子协方差'],['1','协方差有效'],['2','共同顶点状态更新失败'],['3','数值误差传播无法完成'],['4','公共时间导数求解或稳定性检查失败'],['5','协方差含非有限数值'],['6','协方差未通过半正定检查'],['7','最终母粒子传播不支持所配置的过程噪声'],['8','中心解成功，但 24 次正负扰动中至少一侧拟合失败']],[46,453]),
    Spacer(1,10),
    P('motherJacobianFirstFailureFitStatus 另使用 VertexFitStatus 编码：0 为未运行，1 为收敛，2 为几何退化，3 为输运失败，4 为数值失败，5 为达到最大迭代次数。它描述首个失败扰动的顶点求解，与上表不是同一套编号。'),
    P('正负失败掩码的 bit j 对应第 j 个输入分量，0-5 是第一条track的 x,y,z,pₓ,pᵧ,p_z，6-11 是第二条。一个分量的正负两侧均失败时，失败数增加 2；FailedColumn 与 FailedSide 仅指出最先遇到的一侧。'),
    P('若大量失败集中在同一动量或开角区间，应先检查输入、几何退化、最近点分支和差分步长。放宽物理 cut 不会修复一次没有完成的数值求解。'),
    src='[Y] MotherCovarianceStatus；CeeV0LiteVertexFitter.h：VertexFitStatus；[V] RecalculateMotherCovariance。')
