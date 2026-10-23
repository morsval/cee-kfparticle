page('为什么再建立一个八维母粒子副本',
    P('前面的七维拟合回答了一个问题：这两条track能否来自同一个衰变点，它们合成的母粒子状态是什么。进一步假设这个母粒子来自碰撞主顶点时，还能用 PV 的位置和误差信息调整它的状态。'),
    P('v0.6.2 为这一步建立独立的 motherTopo。原 mother 继续保存未加 PV 约束的结果；新副本保存增加“来自该 PV”条件后的估计。这样可以分别查看原质量谱与约束带来的变化。'),
    table([['对象','包含什么','当前用途'],['Candidate.mother','原始 SV、P、E 与 7×7 协方差','原质量谱和原有选择'],['motherTopo','约束后 r、p、E、S 与 8×8 协方差','端点输运、飞行量和附加分析'],['topoFlight','由新副本计算的质量、L、cτ 与误差','读取时须检查 valid']],[112,209,178]),
    Spacer(1,12),
    P('是否加入 S，与磁场是否沿 y 没有先后依赖。当前 CEE 重建的 Λ、反 Λ 和 K⁰ₛ 都是中性母粒子，在这个输运模型中沿直线飞行；它们的带电子粒子仍使用前文的 By 螺旋线。'),
    note('本次升级的实际价值','S 与相关协方差把产生端、衰变端和动量联系起来。只存一个 cτ 中心值无需八维状态；要正确传播端点变化和飞行量误差，则需要同时保存这些关联。'),
    src='[Y] ProductionState / FlightObservables；[F] BuildProductionCopy；[PV8] FitProductionVertex。')

page('S 描述同一粒子的哪一段飞行',
    P('先讨论中性分支，产生端到衰变端的关系可以写成下式。两端都是同一个约束副本中的估计，不能把原七维 SV 与新副本的动量混在一起使用。'),
    eq(r'\mathbf r_D=\mathbf r_P+S\mathbf p'),
    eq(r'S=\frac{\mathbf p\cdot(\mathbf r_D-\mathbf v_{\rm fit})}{|\mathbf p|^2},\qquad L=S|\mathbf p|'),
    P('这里 v_fit 是联合拟合中的产生顶点估计。约束收敛后，两端连线与 p 平行；程序将衰变状态沿 −S 输运到产生端。有限数值容差内，返回的产生端与 v_fit 重合。输入事件中的原 PV 对象不会被改写。'),
    table([['量','数值单位','含义'],['S','cm/(GeV/c)','完整产生→衰变区间的有符号长度／动量'],['L','cm','沿母粒子飞行方向的有符号长度'],['cτ','cm','固有衰变长度 Sm'],['decayTimeNs','ns','原算法的实验室系衰变时刻']],[95,125,279]),
    Spacer(1,10),
    P('例如 p=1 GeV/c、L=5 cm，则 S=5 cm/(GeV/c)。若重建质量为 1.115 GeV/c²，采用 c=1 的数值约定得到 cτ=5.575 cm。固有时间还需除以光速，约为 0.186 ns。'),
    P('S 可以为负，用于表达方向不相容的几何情况。本接口不把负值截成零；现有候选选择通常会先排除明显反向的组合。无效状态的 S 则为 NaN，不能解释成零寿命。'),
    src='[PV8] FitProductionVertex / GetFlightObservables；示例数值仅用于单位说明。')

page('新副本的完整八维协方差',
    P('行列顺序固定为 (x,y,z,px,py,pz,E,S)。刚完成拟合时，x、y、z 指产生端位置；调用 TransportToDecayVertex 后，它们才表示新副本的衰变端。'),
    Covariance(8),
    P('最后一行不是附属信息。C_SS 决定 S 的误差；C_S,p 决定动量误差与飞行参数误差会叠加还是抵消；C_S,E 与质量相关项一起进入 cτ 的误差。'),
    P('例如产生端沿前进方向偏移一点，估计的飞行距离通常会减少。因此产生端位置与 S 可以负相关。若忽略这一项，输运回衰变端时就会得到错误的位置误差。'),
    note('数值维度不等于独立测量的数目','第八维由已有位置、动量及其相关误差导出。没有增加一台测量 S 的探测器，也没有把已知的标准寿命当作输入。'),
    src='[M] Mat8；[Y] ProductionState；[PV8] 由联合后验生成 S 行列。')

page('为什么拟合时先用十维联合输入',
    P('要处理 PV 与母粒子的关联，不能把 PV 只当成一个外部坐标。程序先把七维母粒子与三维 PV 并列成一个十维变量，再对它们之间的几何关系施加约束。'),
    eq(r'\mathbf X=(\mathbf r_D,\mathbf p,E,\mathbf v)^T'),
    table([['协方差块','母粒子 7 列','PV 3 列'],['母粒子 7 行','C_m：母粒子边缘协方差','B：Cov(mother,PV)'],['PV 3 行','Bᵀ','V：PV 协方差']],[125,190,184]),
    Spacer(1,12),
    P('边缘协方差表示已经考虑到 PV 不确定度以后，母粒子自身的总误差。交叉块 B 记录两者一起变化的部分。若母粒子会随 PV 同向平移，差值 SV−PV 的误差会部分抵消；只保留两个对角块就看不到这种抵消。'),
    P('独立接口 ProductionVertexInput 允许调用者给出 B。B 全为零仅对应独立输入假设，不是无条件正确的默认物理事实。整个十维矩阵必须对称、有限且半正定，否则拟合返回明确的失败状态。'),
    note('输入 PV 没有变成无误差的真值','几何约束可以要求track通过同一个产生点，但这个点的位置仍有测量误差。联合拟合允许 PV 估计和母粒子估计在各自误差范围内同时调整。'),
    src='[PV8] ProductionVertexInput、ValidCovariance、FitProductionVertex。')

page('TOF 已经用过 PV，怎样避免漏算相关性',
    P('在 TOF 模式中，PV 先参与母粒子飞行时间的估计，再影响子粒子的 TOF 动量，最后影响拟合出的母粒子。因此，即使原始track测量与 PV 独立，算法输出也会与 PV 相关。'),
    eq(r'D=\frac{\partial\mathbf y_7}{\partial\mathbf v}'),
    P('这个 7×3 矩阵描述母粒子对 PV 的响应。BuildProductionCopy 对每个有误差的 PV 坐标做正负扰动，重新执行完整 SolvePair。用步长 h 与 h/2 分别求导并检查稳定性；三个坐标最多增加 12 次 TOF 求解。失败时会标记新副本无效，不填一组看似正常的零误差。'),
    eq(r'C_m=C_{\rm cond}+DVD^T,\qquad B=DV'),
    P('C_cond 沿用原七维母粒子误差，包括按原配置加入的公共 T0 项；新增 DVDᵀ 表示 PV 不确定度传入母粒子的部分。交叉块 DV 必须同时加入，否则两端之间的共同变化仍会被算错。'),
    sub('这个处理仍有输入假设'),
    P('当前 Finder 假定原始track、PV 和 T0 测量之间没有额外的上游相关性。若同一对子track参与过 PV 拟合，或者上游输出已受 PV 约束，这个假设可能不成立；需要另行提供相应联合协方差。'),
    P('本次补入的是 TOF 算法对 PV 的响应，没有重做全部原始 TOF hit、路程和方向误差的完整联合模型。新误差仍是给定输入模型下的一阶近似。'),
    src='[F] BuildProductionCopy / SolvePair；步长和稳定性阈值见 ProductionVertexConfig。')

page('两个产生顶点约束怎样写成方程',
    P('三维中，一条方向已知的直线能够沿自身方向移动；要求它通过某个点，只提供两个独立的横向条件。本版用平滑的解析表达式构造这两个约束，避免每轮重新选横向基造成导数跳变。'),
    P('先由初始动量选定单位方向 n，以及与它垂直的两个单位向量 e₁、e₂。它们在本次拟合中保持固定。记 d=SV−PV，对 a=1、2 定义：'),
    eq(r's_{\rm plane}=\frac{\mathbf n\cdot\mathbf d}{\mathbf n\cdot\mathbf p}'),
    eq(r'g_a=\mathbf e_a\cdot\mathbf d-s_{\rm plane}\,\mathbf e_a\cdot\mathbf p=0'),
    P('g₁、g₂ 同时为零时，d 与 p 平行。s_plane 是约束方程里的辅助参数；在约束面上它等于最终 S，计算过程中不应把它与保存的 S 混为一谈。'),
    eq(r'\mathbf w_a=\mathbf e_a-\frac{\mathbf e_a\cdot\mathbf p}{\mathbf n\cdot\mathbf p}\mathbf n'),
    eq(r'H_a=(\mathbf w_a^T,\ -s_{\rm plane}\mathbf w_a^T,\ 0,\ -\mathbf w_a^T)',14),
    P('H 的两行分别对应两个约束，十列依次对应 SV、p、E、PV。对 E 的直接导数为零；但如果 E 与位置、动量相关，它仍可以通过增益矩阵获得更新。'),
    src='[PV8] Constraint；方向接近本地参数化奇点时返回失败。')

page('迭代重线性化，为什么不能反复加入同一个 PV',
    P('位置与动量之间的平行关系是非线性的，所以程序在当前估计附近求 H，再更新估计。这里可以理解为迭代扩展卡尔曼式约束拟合；它没有新增逐时刻探测器测量序列。'),
    P('记原联合先验为 X₀、C₀，第 k 轮状态为 Xₖ。每轮都使用同一个 C₀，而不是把上一轮收缩后的协方差当成新的独立先验。'),
    eq(r'W_k=H_kC_0H_k^T,\qquad K_k=C_0H_k^TW_k^{-1}'),
    eq(r'\boldsymbol{\nu}_k=\mathbf g(\mathbf X_k)+H_k(\mathbf X_0-\mathbf X_k)',15),
    eq(r'\mathbf X_{k+1}=\mathbf X_0-K_k\boldsymbol{\nu}_k'),
    P('这样做是在反复改进同一个问题的线性化。如果每轮都把同一个 PV 当作一份新测量，误差会随着迭代次数被人为压小。'),
    eq(r'C_{\rm post}=(I-KH)C_0(I-KH)^T'),
    P('几何条件采用精确约束，所以此处额外的 R 为零；PV 本身的测量误差已经包含在 C₀ 中。收敛要求状态变化和几何残差同时足够小，默认最多 30 轮。'),
    P('在收敛线性化处，新增 χ²=νᵀW⁻¹ν，自由度为 2；新副本的总 χ² 再加原几何拟合 χ²。它与旧 topologicalChi2Ndf 来自不同计算，不能仅凭同名“拓扑”直接交换。'),
    src='[PV8] FitProductionVertex / Gain；协方差使用 Joseph 形式，原输入保持不变。')

page('从十维联合后验得到 S 及其误差',
    P('约束完成后，程序先在新的衰变端组成 (SV,p,E,S)。前三维取自拟合后的 SV，动量和能量也取自同一份后验；S 用两端差值沿动量的投影计算。'),
    eq(r'S=\frac{\mathbf p\cdot\mathbf d}{p^2},\qquad \mathbf d=\mathbf r_D-\mathbf v_{\rm fit}'),
    eq(r'\frac{\partial S}{\partial\mathbf r_D}=\frac{\mathbf p^T}{p^2},\qquad\frac{\partial S}{\partial\mathbf v}=-\frac{\mathbf p^T}{p^2}',15),
    eq(r'\frac{\partial S}{\partial\mathbf p}=\frac{\mathbf d^T-2S\mathbf p^T}{p^2},\qquad\frac{\partial S}{\partial E}=0',15),
    P('将这些导数放在 8×10 变换矩阵 T 的最后一行；前七行取出相应的母粒子变量。然后统一传播：'),
    eq(r'C_8=TC_{\rm post}T^T'),
    P('这一乘法同时给出 S 的方差及它与全部七个变量的协方差。它使用的 PV 与 SV 已经在同一个联合后验中，因此不能再在最后简单加一次“PV 误差”，否则会重复计入。'),
    note('约束后的衰变点可以改变','新副本的衰变端不要求逐坐标等于原七维 SV。那是吸收 PV 条件后得到的另一份估计；原 SV 保存在原 mother 中，仍可用于比较。'),
    src='[PV8] FitProductionVertex 中的 8×10 矩阵 t；下标 7 对应 S。')

page('在产生端和衰变端之间往返输运',
    P('得到衰变端八维状态后，程序先向后输运到产生端再保存。对中性粒子，没有材料效应时，两个方向的中心值变换都很直接：'),
    eq(r'\mathbf r_P=\mathbf r_D-S\mathbf p,\qquad\mathbf r_D=\mathbf r_P+S\mathbf p'),
    P('p、E、S 的中心值在这次端点变换中保持不变。位置行的导数则必须同时包含对动量和 S 的依赖：'),
    eq(r'\frac{\partial\mathbf r_{\rm new}}{\partial\mathbf r}=I,\quad\frac{\partial\mathbf r_{\rm new}}{\partial\mathbf p}=\pm SI,\quad\frac{\partial\mathbf r_{\rm new}}{\partial S}=\pm\mathbf p',14),
    eq(r'C_{\rm new}=JC_{\rm old}J^T'),
    P('若把 S 当作一个没有误差的常数，就会漏掉最后一列及 C_S,p 等交叉项。中心位置可能仍然算对，返回的协方差却已不再对应同一份信息。'),
    code('ProductionState copy = candidate.motherTopo;\nif (copy.valid) {\n  TransportToDecayVertex(copy);\n  TransportToProductionVertex(copy);\n}'),
    P('状态中的 atProductionVertex 标明当前端点。重复请求已经所在的端点不会再次移动。往返检查应同时比较状态和完整矩阵；只检查位置返回原处还不够。'),
    src='[PV8] TransportEndpoint / TransportToDecayVertex / TransportToProductionVertex。')

page('L 与 cτ 的误差为什么都要读第八行',
    P('当前副本的质量由 E 和 p 计算；没有输入标准寿命。L 是有符号飞行长度，cτ 是固有衰变长度。它们的误差按同一份八维协方差传播。'),
    eq(r'm=\sqrt{E^2-p^2},\qquad L=Sp,\qquad c\tau=Sm'),
    eq(r'\sigma_L^2=p^2C_{SS}+S^2\sigma_p^2+2Sp\,\operatorname{Cov}(S,p)',15),
    eq(r'\sigma_{c\tau}^2=m^2C_{SS}+S^2\sigma_m^2+2Sm\,\operatorname{Cov}(S,m)',15),
    P('其中 p=|p| 是模长。Cov(S,p) 由三个 C_S,p_i 按方向投影；质量相关项同时包含动量与能量：'),
    eq(r'\operatorname{Cov}(S,m)=\frac{EC_{S,E}-\sum_i p_iC_{S,p_i}}{m}',16),
    P('代码为 m、L、cτ 分别构造八维梯度，再计算梯度与协方差的二次型。这样同一条规则就包含了所有交叉项；没有把各项误差当成独立误差随意相加。'),
    note('使用时先看有效标志','质量平方非正、矩阵无效或方差计算失败时，派生量不输出有效数值。非常小的负方差只有在明确的浮点舍入容差内才会按零处理。'),
    P('这些误差是局部一阶结果。位移很大、输入误差很宽或估计器明显非线性时，需要额外的模拟检验；维度增加本身不能保证 pull 分布正确。'),
    src='[PV8] GetFlightObservables / Error；ctauCm、sigmaCtauCm 均以 cm 保存。')

page('Finder 实际何时调用新拟合',
    P('v0.6.2 在候选通过原有选择后才调用 BuildProductionCopy。新计算默认开启，但它是否成功不会改变原七维候选是否被接受；因此可以在相同样本中比较约束前后结果。'),
    table([['顺序','本版本执行内容'],['1','SolvePair：得到原始共同顶点与母粒子'],['2','计算七维协方差、共享 T0 项及原拓扑量'],['3','执行原 mass、DCA、χ²、长度与 pointing 等选择'],['4','准备含 PV 相关性的联合输入'],['5','FitProductionVertex：生成 motherTopo'],['6','GetFlightObservables：生成 topoFlight，写入输出']],[60,439]),
    Spacer(1,12),
    P('若只想复现原算法的计算开销，可在运行时加 --no-production-state。C++ 对应 FinderConfig::buildProductionState=false。关闭后新状态标记为 Disabled，不能把空结果当作物理零值。'),
    P('原母粒子协方差无效的 QA 行不会得到有效八维副本。TOF 外层未收敛时，即使原配置允许保留该候选，新副本仍返回 PvResponseFailed。'),
    note('新旧选择的对应关系','现有 --max-topo-chi2ndf 仍作用于原 topologicalChi2Ndf。新 motherTopoChi2Ndf、S、cτ 和新 L/σL 暂未加入自动选择；以后改用它们前应重新研究效率和背景。'),
    src='[F] BuildCandidate / BuildProductionCopy；[C] --no-production-state；[A] FillRow。')

page('新增 ROOT 分支怎样读取',
    P('本版输出格式为 6，保留原分支，并新增以下字段。原 motherCovariance[28] 与新的 motherTopoCovariance[36] 表示不同对象，不能相互覆盖。'),
    table([['新分支','含义'],['motherTopoValid / Status','八维副本是否有效；失败原因'],['motherTopoCorrelationModel','0 外部联合输入；1 tracking 独立 PV；2 TOF PV 响应'],['motherTopoParameters[8]','产生端的 x,y,z,px,py,pz,E,S'],['motherTopoCovariance[36]','对应八维协方差，下三角逐行'],['motherTopoS / SigmaS','S 与标准差'],['motherTopoPvChi2 / PvNdf','新 PV 约束的 χ² 与自由度'],['motherTopoChi2 / Ndf / Chi2Ndf','加上原几何拟合后的统计量'],['motherTopoFlightValid','派生物理量是否有效'],['motherTopoMassGeV / SigmaMassGeV','约束后的质量与误差'],['motherTopoDecayLengthCm / SigmaDecayLengthCm','有符号 L 与误差'],['motherTopoCtauCm / SigmaCtauCm','固有衰变长度及误差']],[284,215]),
    Spacer(1,8),
    P('表中斜线后的名称沿用 motherTopo 前缀，例如 Status 指 motherTopoStatus。','small'),
    P('另有端点标志、迭代次数及 motherTopoDecayLengthSignificance。36 项的索引为 k=i(i+1)/2+j（i≥j）；S 行是 28–35，35 为 C_SS。完整字段列表见 MIGRATION_0_6_2.md。','small'),
    src='[A] CandidateRow / CreateCandidateBranches / FillRow；参数顺序、单位和模型写入运行摘要。')

page('新状态无效时，怎样判断停在哪一步',
    P('motherTopoValid 为 0 时，先看 motherTopoStatus。原候选可能仍然有效；新状态的参数、矩阵和飞行量则保持 NaN。不要只按原 motherCovarianceValid 来筛选新寿命变量。'),
    table([['值','状态','含义'],['0','NotRun','还没有调用'],['1','Valid','八维状态有效'],['2','Disabled','配置关闭新计算'],['3','InvalidMother','原母粒子或协方差无效'],['4','PvResponseFailed','TOF 未收敛或 PV 导数失败／不稳定'],['5','InvalidJointCovariance','联合矩阵非法'],['6','SingularConstraint','约束方差不可逆或方向参数化失效'],['7','NotConverged','约束迭代未满足终止条件'],['8','UnphysicalState','质量、误差或端点输运无效'],['9','UnsupportedCharge','旧版拒绝带电输入；本版保留编号'],['10','InvalidConfiguration','配置或几何统计量非法']],[38,166,295]),
    Spacer(1,12),
    P('v0.6.3 新增状态 11（带电调用漏传磁场）和 12（超出局部飞行范围），见带电扩展一章。','small'),
    P('没有人为保证每个原候选都能得到有效的新状态。分析时应先统计这些状态的比例，再检查它们随动量、衰变长度及输入误差的变化，避免只画“成功候选”的分布而忽略选择效应。'),
    src='[Y] ProductionFitStatus；[F] BuildProductionCopy；[PV8] FitProductionVertex。')

page('这次验证了什么，还不能据此推断什么',
    P('新增代码已经经过可手算情形和独立算法的双重核对。已知动量时，横向位置应成为按方差加权的共同估计；PV 和 SV 共同平移时，误差应在两端差值中抵消。这两类检查直接检验相关项的符号。'),
    table([['检查','解决的问题'],['test_production_state.cxx','解析 χ²、完整协方差、负 S、寿命误差、端点往返及异常状态'],['verify_production_reference.py','以独立八参数广义最小二乘，核对十维约束算法的状态、矩阵与 χ²'],['test_end_to_end.cxx','新状态接入 TOF / tracking；原候选和选择不被覆盖；失败标志明确'],['分析层 check','配置、语法和真实 CLI 参数解析']],[240,259]),
    Spacer(1,12),
    P('独立广义最小二乘程序直接拟合 (r_P,p,E,S)，预测十维测量 (r_P+Sp,p,E,r_P)。它使用另一种参数化和求解器，能够发现把 S 或相关协方差拼错的实现问题。'),
    P('这些检查还不能替代真实 CEE/ROOT 环境的链接、文件读写、探测器模拟 pull 与效率研究。上游遗漏的测量相关性、材料效应和非均匀磁场图也不会因为使用八维状态而自动补齐。'),
    note('与 Ξ⁻ 两级衰变的关系','v0.6.3 已补入均匀 By 下带电母粒子的端点输运，但 Λ 与伴随 π 的组合、质量约束和整条级联链反向重拟合仍需分别实现。将来把 Ξ 顶点再用于 Λ 时，必须处理两者共用测量造成的相关性。'),
    src='CeeV0Lite/tests 中上述三个文件；CeeV0Analysis/Makefile 的 check；本版不改 STAR 参考代码。')
