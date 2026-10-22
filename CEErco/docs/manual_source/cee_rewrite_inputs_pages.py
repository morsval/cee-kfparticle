page('先认清输入：探测器量到什么，程序要重建什么',
    P('当前 CEE 分析的目标，是从两条 TPC 重建track中寻找 Λ→p+π⁻。当前时间模型按 Λ 从主顶点产生、飞行一段距离后衰变来求解；TPC 记录的是两条带电的子粒子track，TOF 记录的是它们到达探测器的位置和时刻。程序要把这些测量重新拼成一个共同衰变顶点，并求出母粒子的动量、质量及误差。'),
    P('零磁场下，TPC track近似为直线。方向可以由track拟合得到，动量大小却不能由曲率确定。因此这条分析路径取 TPC 的方向，再在给定子粒子质量假设后，结合 TOF 确定动量大小。'),
    table([['输入对象','当前用途'],['CeeTpcTrack','两条track的参考状态、方向和拟合质量指标'],['CeeTpcVertex','事件重建主顶点及其位置误差'],['CeeiTOFHit / CeeEtofHit','零场路径的到达位置与时间'],['CeeT0Hit','真实数据可选的事件时间原点'],['MCTrack','模拟真值对照；不参与重建选择']],[210,289]),
    note('核心与分析程序的范围不同','CeeV0Lite 核心还包含 Kₛ⁰ 和反 Λ 的通道。当前 CeeV0Analysis.bin 默认只重建并输出 Λ，且两条子track都来自 TPC；本手册不把这些接口能力写成已完成的 Ξ 级联分析。'),
    src='[A] SetupInput、ProcessEvents；[C] AnalysisConfig；[F] Finder::Reconstruct')

page('主顶点从哪里来：先选位置，再解释误差',
    P('主顶点 PV 来自输入文件中的 CeeTpcVertex。即使分析的是模拟数据，程序也使用这个重建顶点，而不会从 MCTrack 取一个更准确的真值顶点。这样，质量和拓扑结果才反映实际重建链的表现。'),
    P('一个事件可以保存多个顶点。适配层先在每个对象内检查坐标，去掉非有限的条目，再选 rank 最大者；分析层随后在不同对象之间继续比较 rank。若 rank 相同，保留先遇到的顶点。默认允许缺少 rank，但它的优先级低于有 rank 的有效条目。'),
    sub('坐标有效，并不表示误差也已保存完整'),
    P('GetXErr、GetYErr、GetZErr 在这里解释为位置标准差，而不是方差。程序将它们平方，组成对角的 PV 协方差。当前数据类没有提供完整的顶点非对角相关。'),
    eq(r'C_{\rm PV}=\mathrm{diag}(\sigma_x^2,\sigma_y^2,\sigma_z^2)',size=17),
    P('某一项误差缺失、非正或非有限时，默认只替换该项。当前回退标准差为 x、y 各 0.10 cm，z 为 0.15 cm。这些数值是配置假设，不能直接当作 CEE 的实测顶点分辨率。'),
    note('如何在输出中识别回退','pvUsedStoredErrors=1 表示三项误差均来自有效的存储值。只要一项使用回退，该标志就为 0。若希望缺误差的事件直接退出，可用 --strict-pv-errors。'),
    src='[A] SelectPrimaryVertex；[D] MakePrimaryVertex、CeeAdapterConfig')

page('track参考状态：位置、方向和误差必须属于同一处',
    P('一条track可以在不同位置表示。把某一处的位置和动量，与另一处的协方差拼在一起，会改变误差的含义。因此，进入 TOF 匹配和进入重建核心时，程序都调用同一个参考状态选择函数。'),
    P('优先使用文件中保存的 vertex seed：其中包含同一参考点的三维位置、三维动量和 6×6 协方差。程序不只检查 HasVertexSeed 标志，还检查数值有限性、动量是否非零以及协方差是否有效。'),
    sub('没有完整 seed 时，零场和磁场的处理不同'),
    P('在零场模式下，默认允许整组回退到 GetX/Y/Z 与 GetPx/Py/Pz，并配上配置中的位置与方向误差。回退位置标准差为 0.10、0.10、0.15 cm，方向标准差为 0.003 rad。若希望只接受完整 seed，可关闭旧状态回退。'),
    P('磁场模式则必须有有效 seed。因为它真正使用拟合动量及完整六维误差，不能用一个任意的方向向量加回退误差来替代。零场即使读到有效 seed，也只使用其中的位置误差和概括后的方向误差，TOF 再提供动量大小的不确定度。'),
    note('有磁默认参考点与 PV 约束是两回事','v0.6.5 在配对前将有效输入沿测量螺旋线移到当前重建 PV 的 DCA 点，同时输运完整状态与协方差。SetVertexSeed 只保存传入值；本次转换在内存中完成，不强迫子 track 穿过 PV。'),
    src='[D] SelectTrackReferenceState、CopyVertexSeed、DirectionSigmaFromMomentumCovariance')

page('TPC 与 TOF 怎样对应：先建立候选，再做唯一匹配',
    P('重建 Λ 之前，必须知道每条 TPC track对应哪一个 TOF hit。当前零场分析按track直线投影到 hit 所在的探测器平面，再计算平面内的距离残差。这个距离只用于判断匹配相容性，不是随后计算速度使用的完整飞行路程。'),
    table([['TOF 区域','当前投影与残差'],['iTOF 模块 0–8、15–24','投影到 hit 的 x 平面，比较 y、z'],['iTOF 模块 9–14；eTOF','投影到 hit 的 z 平面，比较 x、y'],['未识别的 iTOF 模块','沿用旧分析的 xy 平面直线距离']],[210,289]),
    P('距离不超过配置上限的“track—hit”组合进入候选匹配集合。iTOF 与 eTOF 默认上限均为 10 cm。程序把二者放在同一个分配问题中，要求每条track最多用一个 hit，每个 hit 也最多分给一条track。'),
    P('选择分两层：首先尽可能增加成功匹配数；在这个数目下，再减小按探测器距离上限归一化后的总残差。因此它不是按输入顺序逐条抢最近 hit，也没有固定规定先用 iTOF、再用 eTOF。'),
    note('适用范围','这套直线投影用于零场路径。当前磁场分析不借用它匹配 TOF；若将来要用 TOF 校验磁场track，需要另行实现与螺旋track一致的匹配。'),
    src='[G] ProjectionDistanceCm、AssignUniqueTofHits；[A] CollectTofHits、ProcessEvents')

page('质子和 pion 的身份：先是假设，再由候选接受检验',
    P('当前程序并没有先通过一套完整的 dE/dx PID，把每条track确定为质子或 pion。应用层给每条track同样的两种质量假设；dE/dx 会写入输出树，供检查使用，但没有用于当前在线筛选。'),
    P('对一对track，程序分别尝试“第一条是质子、第二条是 pion”和相反分配。每种假设有自己的子粒子质量，因此同样的速度也会得到不同的动量大小：'),
    eq(r'\beta=\frac{\ell}{c\,\Delta t},\qquad p=\frac{m\beta}{\sqrt{1-\beta^2}}',size=18),
    P('这里 TOF 帮助确定的是速度。要从速度得到动量，仍需输入质量假设。不能把这一步写成“TOF 独立测出动量，再唯一认出粒种”。最终还要看共同顶点、母粒子质量及拓扑是否相容。'),
    sub('当前默认还保留旧分析的 DCA 角色规则'),
    P('质子假设要求到 PV 的 DCA≤1.2 cm，pion 假设要求 DCA&gt;1.2 cm。这是当前采用的经验角色分配，会改变候选样本；它不是每个真实 Λ 都必须满足的物理定律。可用 --disable-dca-role 做关闭后的对照。'),
    note('零场电荷标志的含义','零场不能从曲率测出电荷正负，适配层将其标记为未知。输出 daughterPdg0=2212、daughterPdg1=-211 表示本次重建假设，是否与真值一致要另看真值匹配列。'),
    src='[A] kLambdaDaughterHypotheses、FillRow；[D] MakeTrackInput；[F] BuildCandidate、Reconstruct；[Y] EstimateTofMomentum')

page('普通模拟的时间原点：t₀=0 是输入约定',
    P('当前 simulation 模式有一个明确前提：TOF 时间已经以发生器的碰撞时刻为原点。因此事件起始时刻设为精确的 0 ns，不读取 CeeT0Hit，也不额外加入事件公共起始时间误差。'),
    P('这与“真实数据没有测到 T0，于是暂时拿 0 代替”不同。后者仍存在时间标定和分辨率问题，不能借用普通模拟的精确时间假设。输入若不满足发生器 t=0 的约定，应先检查时间定义，而不是通过添加一个随意的 T0 误差掩盖差异。'),
    table([['simulation 的时间设置','当前值'],['事件时间原点','fixed，0 ns'],['sigmaStartTimeNs','0'],['additionalDecayTimeSigmaNs','0'],['公共时间误差传播','关闭'],['单条 TOF 时间与路程误差','仍保留；默认 0.10 ns、0.5 cm']],[275,224]),
    P('没有公共 T0 误差，不等于所有时间都没有误差。每个 TOF hit 的时间分辨率和路程误差仍会影响速度与动量；母粒子的衰变时刻也仍需由求解过程确定。程序会检查配置一致性，拒绝在 simulation 模式下混入实测 T0 的单位转换选项。'),
    note('真值只用来核对结果','模拟可读取 MCTrack 标记真 Λ 与组合背景；--no-truth 可关闭。MCTrack 缺失时仍可重建。真实数据路径不绑定真值分支，两种模式都使用重建 PV 和重建track。'),
    src='[C] ConfigureDataMode、ValidateDataModeTiming、main；[A] SetupInput、EventT0')
