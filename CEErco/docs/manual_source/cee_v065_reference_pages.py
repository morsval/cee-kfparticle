class ReferenceDiagram(Flowable):
    def __init__(self): self.width=CW;self.height=185
    def draw(self):
        c=self.canv
        # Local circle segment in the plane transverse to By (schematic).
        center=(205,205);radius=145
        def pt(a):return (center[0]+radius*math.cos(a),center[1]+radius*math.sin(a))
        path=c.beginPath()
        for i in range(101):
            x,y=pt(-2.55+i*.021)
            if i==0:path.moveTo(x,y)
            else:path.lineTo(x,y)
        c.setStrokeColor(TEAL);c.setLineWidth(2);c.drawPath(path)
        r=pt(-.65);a=pt(-1.75);pv=(a[0]-.18*(center[0]-a[0]),a[1]-.18*(center[1]-a[1]))
        c.setDash(3,3);c.setStrokeColor(GOLD);c.line(*pv,*a);c.setDash()
        for p,col in [(r,BLUE),(a,TEAL),(pv,GOLD)]:c.setFillColor(col);c.circle(*p,4,fill=1,stroke=0)
        txt(c,'输入参考点 r',r[0]+8,r[1]+5,10,BLUE)
        txt(c,'track–PV DCA 点',a[0]+10,a[1]-28,10,TEAL)
        txt(c,'重建 PV',pv[0]-67,pv[1]-3,10,GOLD)
        tx,ty=(-math.sin(-1.75),math.cos(-1.75));arrow(c,a[0],a[1],a[0]+45*tx,a[1]+45*ty,TEAL,1.7)
        txt(c,'该点的动量方向',a[0]+45,a[1]+22,9,TEAL)
        txt(c,'参考点沿同一测量螺旋线移动；DCA 通常非零。',8,8,10,GRAY)

page('有磁默认参考：每条 track 先移到 track–PV DCA 点',
 P('PV（Primary Vertex）是当前事件选出的重建主顶点。DCA（Distance of Closest Approach）是 track 到它的最近距离；对应位置是螺旋线上的最近点 POCA（Point of Closest Approach）。新的参考位置位于 track 上，通常不等于 PV。'),
 ReferenceDiagram(),
 table([['顺序','v0.6.5 实际执行'],['1　准备输入','读取同一参考点的 r、p、完整 6×6 协方差及拟合电荷'],['2　统一参考','每条 track 求到当前 PV 的名义 POCA，转换完整状态'],['3　进入配对','使用新参考点处的切线给 SV 初值，继续原共同顶点迭代'],['4　保存候选','原 7D 结果和附加母粒子 PV 约束结果仍分别保存']],[100,399]),
 P('转换放在 Finder::Reconstruct 的配对循环之前，每个输入 track 只做一次。直接调用 Finder、CeeV0Analysis.bin 与配套 .C macro 都使用这个默认行为。零场 TOF 模式跳过这一步。'),
 src='[F] PrepareTrackReference、Finder::Reconstruct；[T] ReanchorAtClosestApproach；CeeTpcTrack.h。')

page('参考点转换怎样传播误差：固定名义弧长',
 P('先在当前测量螺旋线上求出名义弧长 ℓ₀，使到 PV 的距离最小。然后固定这个参考坐标，用固定弧长输运把整组状态及误差移过去：'),
 eq(r'\ell_0=\arg\min_{|\ell|\leq\ell_{\max}}|\mathbf r(\ell)-\mathbf v_{\rm PV}|^2',18),
 eq(r'\mathbf a^{\prime}=f(\mathbf a,\ell_0),\quad F=\left.\frac{\partial f}{\partial\mathbf a}\right|_{\ell_0},\quad C^{\prime}=FCF^T',18),
 P('a=(r,p) 为六维 Cartesian 状态；完整位置–动量交叉块保留。PV 在这里仅用于选择参考点，不作为对子 track 的测量约束，所以不把 PV 误差加进 track 协方差，也不把 DCA 改成零。纯参考转换关闭过程噪声，避免把重复参数化当成新的散射。'),
 table([['操作','Jacobian 的含义'],['新增：输入参考点转换','名义弧长固定；仅改变同一测量 track 的参数表示'],['原有：输运到当前 SV 的 POCA','每个输入扰动都重新求 POCA；包含路径对状态的导数'],['母粒子 PV 约束','使用 PV 的位置及误差更新母粒子；是独立的拟合步骤']],[165,334]),
 note('不要只改 SetVertexSeed 的位置','必须同时更新该位置处的动量和完整协方差。SetVertexSeed 只保存调用者传入值；本版在重建内存副本中统一参考，不改写输入 ROOT 文件。'),
 src='[T] FindClosestApproachParameters、PropagateState、TransportToClosestApproach；[F] PrepareTrackReference。')

page('默认值、失败处理与输出检查',
 table([['项目','默认／含义'],['usePrimaryVertexDcaReference','true；非零 By 时启用'],['弧长搜索范围','相对输入参考点 ±1000 cm；不是相对 PV 的距离'],['求根设置','最多 64 次；驻点残差容差 10⁻⁸ cm'],['必需输入','已知非零电荷、拟合动量、完整有效 6×6 协方差'],['失败处理','不收敛、边界截断、无效输入／协方差：跳过该 track'],['stored 选项','显式使用旧存储参考点，供同 cut 对照']],[175,324]),
 sub('三个新增 QA 直方图'),
 code('hPvDcaReferenceStatus\nhPvDcaReferencePath\nhPvDcaReferenceDca'),
 P('状态 0=未转换，1=DCA 参考成功，2=输入无效，3=输运失败。路径图记录转换的有符号弧长，DCA 图记录成功 track 到 PV 的距离。debug=2 会报告失败 track ID。原 hTrackReferenceSource 仍表示适配层读取 seed／回退的来源。'),
 P('元数据保存参考策略、最大路径、容差、迭代数、固定弧长协方差约定以及 ceeV0LiteVersion=0.6.5。候选树分支布局未变，outputFormat 仍为 6；numberFinderTracks 是送入 Finder 前的数量，参考失败数另看 QA。'),
 note('有界螺旋线的适用范围','算法沿用均匀 By 模型，不包含材料或实际可用 track 段限制。低动量多圈 track 可能选到较远的一圈；搜索半范围可调，必须结合实验接受度验证。'),
 src='[F] FinderConfig、TrackReferenceDiagnostic；[A] SetupOutput、ProcessEvents；main.cxx。')

page('同一批超氚信号重跑：S 长尾减少，但未完全消失',
 P('固定生成种子 20260929，从原 6000 个信号中保留同一批 5444 个接受度／PID 后测量输入。By=+0.5 T，子粒子电荷 +2、−1。只切几何 χ²&lt;10、子 track DCA&lt;1 cm 并要求 PV 拟合有效；未施加质量窗、PV χ²、飞行显著性或 truth 质量 cut。'),
 eq(r'S=L_{\rm arc}/|\mathbf p|,\qquad {\rm pull}_S=(S_{\rm rec}-S_{\rm true})/\sigma_S',18),
 table([['共同 5410 个信号','原参考','默认 PV DCA'],['S pull 标准差','1.192124','0.999836'],['|pullS| > 5','9','0'],['最大 |pullS|','25.4803','4.6893'],['最大 |Srec−Strue|','11.2480','4.0428'],['PV 约束前 SV 误差 > 10 cm','10','1']],[245,127,127]),
 P('S 的单位为 cm/(GeV/c)。默认路径共得到 5415 个有效信号，整体 pull 标准差 1.001008。与先前 toy 的 DCA 参考对比，位置、动量与协方差逐项相同，有效候选 ID 完全相同。MC truth 只用于本页诊断。'),
 note('仍然存在错误局部解','候选 2720：Strue=4.33552，默认重建 S=0.292714，σS=0.992641，SV 偏差 21.10 cm。无 >5σ 尾部不等于每个 S 都准确；参考位置改变了切线初值，但没有增加多 seed 搜索或 truth 引导。'),
 P('新增测试覆盖电荷／磁场反号、驻点、完整协方差往返、失败拒绝、自动与手动参考一致、零场旁路。真实 CEE/ROOT I/O 和探测器级闭合尚需在实验环境验证。','small'),
 src='v0.6.5_review/replay_signal.py、validation/signal_replay_summary.json；tests/test_pv_dca_reference.cxx。')

page('运行入口：ROOT 文件或 list 与可调 cut',
 P('配套 runCeeV0Reco.C 放在 CeeV0Analysis/macros。它调用正式分析程序，支持一个 ROOT 文件或每行一个 ROOT 路径的 .list/.txt；输入必须已包含重建 track 和 PV。零场另需 TOF hit。默认树为 ceesim；当前分析输出 Λ→p+π⁻。'),
 sub('有磁模拟示例：在 CeeV0Analysis 目录运行'),
 code('make clean\nmake -j4\nroot -l -b -q \'macros/runCeeV0Reco.C(\n  "/data/files.list","lambda_By.root",-1,\n  "tracking",0.5)\''),
 P('上面为便于阅读分行；实际 shell 命令应写在一行。默认 reference 为 pv-dca，使用测量电荷、动量及完整协方差；不需要在宏里逐个调用 SetVertexSeed。'),
 sub('在宏的配置区调整'),
 code('const char* magneticTrackReference = "pv-dca";\nconst double referenceMaximumAbsPathCm = 1000.;\nconst double maximumDaughterDcaCm = 1.e30;\nconst double maximumGeometricChi2Ndf = 1.e30;\nconst double minimumDecayLengthSignificance = 0.;'),
 P('上述是本宏实际默认值，沿用宽松选择，不因本次参考点更新自动增加 L/σL cut。需要更紧选择时改相应变量；stored 只用于参考点对照。直接运行程序也可用 --magnetic-track-reference pv-dca 与 --reference-max-path-cm 1000。'),
 note('重编译范围','FinderConfig 与 Reconstruct 接口新增字段／参数，核心及其调用程序必须一起清理重编；不能只替换头文件。此次 CeeTpcTrack 只改注释，持久 schema 无新增成员。'),
 src='CeeV0Analysis/macros/runCeeV0Reco.C、README_CN.txt；MIGRATION_0_6_5.md；main.cxx。')
