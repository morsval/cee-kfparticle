class MakerArchitecture(Flowable):
    def __init__(self): self.width=CW; self.height=253
    def draw(self):
        c=self.canv
        rows=[('macro / 用户配置','选择输入、基础 cut、PID 与输出变量',BLUE),
              ('AnalysisMaker','事件 / track 选择 → PID → 选中 track 列表',TEAL),
              ('Interface → Finder','共同 SV、Kalman 更新、7D 母粒子与 8D PV 副本',BLUE),
              ('AnalysisMaker → Writer','用户母粒子 cut → 选定分支或完全自定义 ROOT tree',TEAL)]
        for i,(title,desc,col) in enumerate(rows):
            y=205-i*61
            c.setFillColor(PALE);c.setStrokeColor(GRID);c.roundRect(10,y,CW-20,48,5,stroke=1,fill=1)
            txt(c,title,23,y+30,10.5,col,bold=True);txt(c,desc,23,y+11,9.2,INK)
            if i<3:arrow(c,CW/2,y,CW/2,y-13,GRAY,1.1)

page('v0.7.0：从用户分析类控制重建',
 P('原来有 CeeV0Analysis::Analyzer，能够读 ROOT、选 track、运行 Finder 并写树，但这些职责集中在一起。新版本新增简化的 CeeKFParticleAnalysisMaker 类，使 PID、分析 cut、输入格式和输出结构可以分别修改。'),
 MakerArchitecture(),
 P('设计参照 STAR 的 StKFParticleAnalysisMaker / StKFParticleInterface，但不依赖 StMaker、singleton 或 STAR 的数据格式。macro 直接调用类；共同顶点与 Kalman 拟合仍使用现有 CeeV0Lite 核心。'),
 note('物理算法基线不变','保留 v0.6.5 的有磁 track–PV DCA 默认参考、零场 TOF 迭代、原 7D 母粒子和附加 8D PV 约束副本。本次主要变化是用户分析接口与输出自由度。'),
 src='CeeV0Analysis/include、src/CeeKFParticleAnalysisMaker；Interface；RootSource；RootWriter。')

page('分清五层职责：输入格式只在 reader 中出现',
 table([['组件','负责的事'],['CeeKFParticleRootSource','读 CeeTpcTrack / TOF / PV，转成 AnalysisEvent 值对象'],['CeeKFParticleAnalysisMaker','事件 cut、track cut、PID、选中 track 列表、用户母粒子 cut'],['CeeKFParticleInterface','把选中输入交给 Finder，提供候选与参考点转换诊断'],['CeeV0Lite::Finder','组对子、共同顶点拟合、Kalman 更新、协方差及 PV 副本'],['RootWriter / CandidateSink','选择输出分支，扩展用户变量，或完全自己写 ROOT tree']],[180,CW-180]),
 P('AnalysisEvent 保存 PV、时间信息与一组 AnalysisTrack；不持有 ROOT tree 内部对象指针。用户以后换成 picoDst 类输入，只需要写新的 CeeKFParticleEventSource，填好这些结构。Maker、Finder 和 Writer 不需要知道 picoDst 的类名。'),
 P('每条有磁输入 track 的位置、物理动量与完整 6×6 协方差必须位于同一参考点，并带有实测电荷。Finder 在 PID 后、组对前统一转换到 track–PV DCA 参考；它不是对子 track 加 PV 拟合约束。'),
 src='CeeKFParticleEvent.h；CeeKFParticleRootSource.h/.cxx；CeeV0LiteFinder.h/.cxx。')

page('运行新 macro：直接构造并调用 Maker',
 code('cd CEErco/CeeV0Analysis\nmake shared\nroot -l -b -q \\\n  \'macros/runCeeKFParticle.C("input.list","mothers.root",0.5)\''),
 table([['参数','意义'],['input / output','一个 ROOT 文件或 .list/.txt；母粒子输出文件'],['ByTesla','带符号均匀 By；0 为零场 TOF，非零用 tracking 动量'],['maximumEvents','读取输入条目上限；-1 为全部'],['realData','false：模拟精确 t=0；true：使用明确的实际时间配置'],['overwrite','默认 false；显式 true 才覆盖同名输出']],[140,CW-140]),
 code('CeeKFParticleRootSource source(config);\nUserCeeAnalysis maker(config.finder);\nCeeKFParticleRootWriter writer("mothers.root");\nstd::string error;\nif (!maker.Run(source, writer, error))\n  throw std::runtime_error(error);'),
 P('共享库需要已配置的 ROOT / FairRoot / CEE 环境。完整例子在 macros/runCeeKFParticle.C；上面短片段省略了输入、磁场和 cut 配置。列表每行一个 .root 路径，支持空行及 # 注释；相对路径按运行目录解释。真实零场数据的时间单位、偏置和误差需要使用实验标定。','small'),
 src='CeeV0Analysis/macros/runCeeKFParticle.C；Makefile；AnalysisMaker_v0.7.0.md。')

page('PID 与 track cut：用户来定义粒子假设',
 P('AnalysisTrack 提供 nHits、chi2Ndf、fitConverged、原始 dEdx，以及可选的 nSigmaProton、nSigmaPion、tofMass2GeV2。当前 reader 没有相应 PID 标定时，后三者为 NaN；不能把它们理解为已测得的零值。'),
 code('unsigned IdentifyTrack(const AnalysisTrack& t,\n                       const AnalysisEvent&) const {\n  unsigned mask = 0;\n  if (std::isfinite(t.nSigmaProton) &&\n      std::fabs(t.nSigmaProton) < 3.)\n    mask |= kProtonHypothesis;\n  if (std::isfinite(t.nSigmaPion) &&\n      std::fabs(t.nSigmaPion) < 3.)\n    mask |= kPionHypothesis;\n  return mask;\n}'),
 table([['返回值','行为'],['0','不用于重建'],['一个 hypothesis bit','仅保留该质量假设'],['两个 bit','保留 PID 歧义，允许两种质量假设']],[168,CW-168]),
 P('该例要求你先提供标定后的 nSigma，否则会拒绝所有 track。默认宏保留输入 mask，并未替你完成 PID。也可使用 SetPIDSelector callback；MakerTrackCuts 控制基础 hits / chi2 / convergence，AcceptTrack 可写更完整选择。实测电荷仍由 Finder 检查，不使用 MC 真值决定 PID。','small'),
 src='CeeKFParticleAnalysisMaker::AcceptTrack、IdentifyTrack；CeeKFParticleEvent.h。')

page('预先选择 track：减少后续组合与拟合',
 table([['选择位置','索引含义与作用'],['source.SetTrackIndexSelector','原 ROOT TClonesArray 索引；在适配、TOF 唯一匹配之前筛选'],['maker.SetTrackListSelector','AnalysisEvent.tracks 的索引；只对这些对象运行 PID 并送入 Finder'],['直接填 AnalysisEvent.tracks','未来 reader 或上游已有选中列表时，不需要再提供所有 track']],[190,CW-190]),
 P('两套索引不要混用。原输入 ID 存在 TrackInput.id 中，选择后不重新编号。空列表表示零条 track；重复、越界索引或重复/负的选中 track ID 会报错。每次 Make 都清空上一事件的候选和诊断。'),
 code('maker.SetTrackListSelector(\n  [](const AnalysisEvent& e) {\n    return std::vector<std::size_t>{1, 5};\n  });'),
 P('示例数字只是说明索引写法，实际列表必须来自当前事件。全事件输入通常仍需首次线性遍历才能做 PID；节省的是后面的 pair 和拟合。未在本版宣称已有实测加速倍数。'),
 note('零场 TOF 的选择时机','唯一 hit 分配在 reader 内完成。若某些 track 不应参与分配，应在 source 索引选择器中提前排除；Maker 的后续 PID 不会自动重新匹配。新 reader 超过 maximumTracksPerEvent 会报错，不按顺序静默截断。'),
 src='RootSource::Next；AnalysisMaker::Make；Interface::Reconstruct。')

page('母粒子 cut：基础筛选与分析筛选分开',
 P('Finder 的 lambdaCuts / k0Cuts 在重建内部作用于原始候选，提前拒绝的对象不会进入用户分析层。附加 PV 副本建立之后，Maker 的 AcceptMother 才能对完整候选进行分析选择。'),
 code('bool AcceptMother(const Candidate& c,\n                  const AnalysisEvent&) const {\n  if (!(c.massGeV > 1.05 && c.massGeV < 1.25))\n    return false;\n  if (!(DcaProton(c) > 0.1 && DcaPion(c) > 0.1))\n    return false;\n  if (!c.motherTopo.valid) return false;\n  const auto& p = c.motherTopo;\n  return p.pvConstraintNdf > 0 &&\n         p.pvConstraintChi2 / p.pvConstraintNdf < 5.;\n}'),
 P('这些数值只展示写法，不是经过验证的 CEE 工作点。需要放宽选择时，先检查是否已被 Finder 的基础 cut 拒绝。新示例宏默认不加 L/sigmaL 或 S/sigmaS cut。'),
 table([['接口','得到什么'],['maker.GetCandidates()','经过用户母粒子选择的当前事件候选'],['maker.GetInterface().GetCandidates()','通过 Finder 基础 cut、尚未经过用户母粒子选择的候选']],[246,CW-246]),
 P('返回的是只读引用，下次 Make 会更新；跨事件保留时请复制。可能无效的量先用 isfinite 或 !(x &gt; cut) 检查，防止 NaN 从只检查 x &lt; cut 的写法中漏过。','small'),
 src='Finder::BuildCandidate；AnalysisMaker::AcceptMother、GetCandidates。')

page('DcaProton 和 DcaPion：原来就有底层量',
 table([['你要的名字','原始存储与当前查询'],['DcaProton(c)','daughterDcaToPvCm 中 |PDG|=2212 的一项'],['DcaPion(c)','daughterDcaToPvCm 中唯一 |PDG|=211 的一项'],['旧 ROOT 分支','daughterDcaToPv0Cm / daughterDcaToPv1Cm'],['Lambda 默认角色','daughter 0 是 proton，daughter 1 是 pion']],[170,CW-170]),
 eq(r'DCA_i=|\mathbf r_{i,\rm POCA(PV)}-\mathbf v_{\rm PV}|',18),
 P('核心使用共同顶点更新前、该质量假设下的 inputStates，输运到 PV 的 POCA 再算距离。有磁场沿螺旋线，零场沿直线；这是无符号三维距离，单位 cm。零场动量来自最终 TOF 自洽解。'),
 P('它不是两条 daughter 相互的 DCA，也不是共同顶点 Kalman 更新后的 daughter 再加 PV 约束。新增角色查询适用于反 Lambda；K0S 有两条 pion，DcaPion 因角色不唯一返回 NaN，请使用显式 0/1 变量。'),
 note('与 chi2Primary 的区别','Chi2PrimaryProton/Pion 使用 daughterChi2FromPv，是带 track / PV 协方差的原始 chi2；daughterPvSignificance 才是平方根。类似 STAR 的名字不代表误差模型或阈值可以直接照搬。'),
 src='Finder::MeasureDaughterAtPv、BuildCandidate；CeeKFParticleObservables.cxx。')

page('输出树由你决定：选列、加列或完全自写',
 code('writer.SetVariables({"massGeV", "DcaProton",\n  "DcaPion", "motherTopoValid", "topoMassGeV"});\nwriter.AddColumn("massMeV",\n  [](const Candidate& c,const AnalysisEvent&) {\n    return 1000. * c.massGeV;\n  }, "MeV/c^2");\nwriter.SaveCovariances(true);'),
 P('SetVariables 替换整个列表，之后再 AddColumn。未知或重复名字报错，Begin 后锁定。固定 eventId 为 64 位整数，candidateIndex 为整数；可选标量和用户列为 double。元数据保存版本、单位、阶段、列含义和用户配置文本。'),
 table([['可选矩阵分支','顺序'],['motherCovariance[49]','7×7，x,y,z,px,py,pz,E'],['motherTopoCovariance[64]','8×8，在上述顺序后增加 S'],['daughter0/1Covariance[36]','各 6×6，x,y,z,px,py,pz']],[225,CW-225]),
 P('矩阵按 row-major 展平。无效母粒子/PV 副本的矩阵及相应派生量写 NaN，不沿用上一候选的值。完整 Candidate 仍可直接取得。'),
 note('完全自定义树','实现 CandidateSink 的 Begin / Write / Finish；examples/UserTreeSink.h 演示自己定义 branch 和 Fill。默认输出是一候选一行；接受事件即使没有候选也调用 Sink::Write，方便你另存事件计数或 event 数组。'),
 src='CeeKFParticleRootWriter.h/.cxx；examples/UserTreeSink.h。')

page('换输入、迁移旧分析与验证边界',
 P('examples/PicoDstAdapterExample.h 提供一个可编译的 callback reader：由用户读取未来的数据格式，填 AnalysisEvent 和选中的 AnalysisTrack。它没有假定一个目前尚不存在的 CEE picoDst schema。所有 r、p、C 必须位于同一参考点。'),
 table([['保留 / 新增','范围'],['旧 Analyzer / CLI','保留固定 Lambda 树、QA、truth 匹配与旧 macro'],['新 Maker 入口','Lambda / 反 Lambda / K0S 开关；默认 macro 仅 Lambda'],['超氚','仍为独立 toy/拟合示例；本次没有加入 helium PID 通用通道'],['新 Reader / Writer','不自动复制旧 truth/QA；自定义 Sink 可自行添加']],[145,CW-145]),
 P('验证包含 9 组核心回归与旧分析检查；零场/有磁下 Maker 与直接 Finder 的同输入质量及 7D 协方差完全一致；PID、子集、重复 ID、事件清空、用户 cut 和错误关闭均有覆盖。'),
 P('真实 ROOT 测试验证写出后读回、自定义列、int64 event ID、NaN 清空及协方差；多文件 reader 测试使用 ROOT 加合成 CEE 形状字典。它验证读取路径与接口契约，不能代替真实 FairRoot/CEE 字典及真实实验文件的链接测试。'),
 note('失败文件与输出保护','默认 CREATE 保护已有输出。读写失败时 Run 返回 false 并尝试关闭文件；留下的文件可能是不完整结果。线程间不要共享同一个 Maker / Source / Writer。'),
 src='tests/test_maker.cxx、test_root_writer.cxx、test_root_source.cxx；AnalysisMaker_v0.7.0.md。')

# The machine-readable table and this chapter share the same catalog.
variable_rows=json.loads((BASE.parent/'CandidateVariables_v0.7.0.json').read_text())
for start in range(0,len(variable_rows),9):
    group=variable_rows[start:start+9]
    rows=[['查询 / 输出名','单位','含义与阶段']]
    for r in group:
        name=re.sub(r'([a-z0-9])([A-Z])',r'\1<wbr/>\2',r['name'])
        # ReportLab CJK wrapping already handles long identifiers; no HTML wbr.
        name=name.replace('<wbr/>','')
        stage={'identity':'标识','original':'原候选 / 输入','event':'事件','status':'状态','PV refit':'PV 约束后'}[r['stage']]
        rows.append([name,r['unit'],r['中文说明'].replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')+'<br/><font color="#51626F">'+stage+'</font>'])
    page(f'母粒子变量表 {start+1}-{start+len(group)} / {len(variable_rows)}',
      P('以下名字可直接用于 SetVariables，或 CandidateValue(name,candidate,event)。区分原 7D 候选与附加 PV 约束 8D 副本；PV 副本无效时对应查询值为 NaN。','small'),
      table(rows,[161,75,CW-236]),
      P('S 是飞行参数，L/sigmaL 与 S/sigmaS 是候选飞行显著性，均不是质量谱信号显著度。完整 C++ 表达式见 docs/CandidateVariables_v0.7.0.csv。','small'),
      src='CeeKFParticleObservables::CandidateObservables；CandidateVariables_v0.7.0.csv。')

page('变量表补充：完整状态与诊断仍能直接取',
 table([['Candidate 中的位置','内容'],['mother / motherTopo','7D / 8D 位置、动量、能量、S、完整协方差及状态'],['decayVertexCovarianceCm2','母粒子 PV 约束前 SV 的 3×3 协方差'],['daughterAtVertex[0/1]','共同顶点更新后的 daughter 位置、动量、质量、PDG、电荷和 6×6 协方差'],['mother.jacobian*','扰动拟合次数、成功数、失败列/方向、bitmask 和失败状态'],['daughterDeviationFromPv[0/1]','旧兼容别名，等于 sqrt(chi2)；新代码用明确命名'],['GetSelectedTracks()','送入 Finder 前的用户选中输入与 PID mask'],['GetReferenceDiagnostics()','Finder 内部参考转换状态、路径、DCA 和收敛信息']],[211,CW-211]),
 P('89 个可选标量覆盖常用 cut 与输出。矩阵、完整对象和诊断结构不需要拆成固定标量才能读取；用户可通过 AddColumn 暴露新量，也可在自己的 Sink 中保存数组或对象。'),
 P('对 PV 约束后的状态，先检查 motherTopo.valid；使用长度与误差时再检查 topoFlight.valid。motherTopo.positionCm 通常位于 PV 端点，不能直接把它当成新的 SV。'),
 src='CeeV0LiteTypes.h；CeeKFParticleAnalysisMaker.h；CeeKFParticleInterface.h。')
