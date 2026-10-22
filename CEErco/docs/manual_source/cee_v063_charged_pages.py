class ChargedFlightFigure(Flowable):
    def __init__(self): self.width=CW; self.height=238
    def draw(self):
        import math
        c=self.canv; px=53; py=194; radius=485; angle=.65
        def point(t): return px+radius*math.sin(t), py+radius*(math.cos(t)-1)
        dx,dy=point(angle)
        c.setStrokeColor(GRAY);c.setLineWidth(1);c.setDash(4,3);c.line(px,py,dx,dy);c.setDash()
        path=c.beginPath();path.moveTo(px,py)
        for i in range(1,81): path.lineTo(*point(angle*i/80))
        c.setStrokeColor(BLUE);c.setLineWidth(2.6);c.drawPath(path)
        arrow(c,px,py,px+87,py,TEAL,2)
        arrow(c,dx,dy,dx+72*math.cos(angle),dy-72*math.sin(angle),TEAL,2)
        for x,y in [(px,py),(dx,dy)]:c.setFillColor(BLUE);c.circle(x,y,4.5,fill=1,stroke=0)
        txt(c,'PV：产生端',px-18,py+20,11,BLUE,bold=True)
        txt(c,'SV：衰变端',dx-8,dy+18,11,BLUE,bold=True)
        txt(c,'产生时动量',px+34,py+7,9,TEAL)
        txt(c,'衰变时动量',dx+20,dy-49,9,TEAL)
        txt(c,'超氚的弧线',177,181,10,BLUE,bold=True)
        txt(c,'两端直线连线',153,135,10,GRAY)
        txt(c,'By > 0，q = +1；z 向右、x 向上，弯转作放大示意。',15,12,9,GRAY)

class ChargedArcMatrix(Flowable):
    def __init__(self):self.width=CW;self.height=180
    def draw(self):
        c=self.canv;left=54;cw=(CW-65)/3;top=145;rh=38
        entries=[[r'\sin\theta/\kappa','0',r'(1-\cos\theta)/\kappa'],['0','S','0'],[r'-(1-\cos\theta)/\kappa','0',r'\sin\theta/\kappa']]
        for j,l in enumerate(['p_{Dx}','p_{Dy}','p_{Dz}']):drawmath(c,l,left+(j+.5)*cw,top+12,12,anchor='center')
        for i,l in enumerate(['x','y','z']):
            y=top-(i+1)*rh;drawmath(c,l,29,y+14,12,anchor='center')
            for j in range(3):
                c.setFillColor(PALE if i==j else WHITE);c.setStrokeColor(GRID);c.rect(left+j*cw,y,cw,rh,fill=1,stroke=1)
                drawmath(c,entries[i][j],left+(j+.5)*cw,y+13,13,anchor='center',maxw=cw-9)
        txt(c,'A(S)：乘以衰变端动量，得到产生端到衰变端的位移。',17,5,9,GRAY)

page('超氚也带电，母粒子的飞行track必须弯转',
    eq(r'{}^3_{\Lambda}{\rm H}^{+}\rightarrow{}^3{\rm He}^{++}+\pi^-'),
    P('这个衰变的电荷关系是 +1=+2−1。两条子track在磁场中弯转，衰变以前的超氚也会弯转。用子track找到 SV 以后，必须用母粒子自己的电荷和动量，把它的螺旋线追溯到 PV。'),
    ChargedFlightFigure(),
    P('图中产生时与衰变时的动量方向不同。两端直线连线通常既不平行于产生端切线，也不平行于衰变端切线。因此，对带电超氚直接要求 (SV−PV) 与 p 平行，会施加错误的几何条件。'),
    note('此前两个版本的区别','正式 v0.6.2 明确拒绝带电母粒子的 PV 拟合；2026-09-29 的独立超氚测试另用 Python 螺旋线约束处理这一步。v0.6.3 将相应能力合入正式 C++ 核心，并与原独立测试逐候选核对。'),
    src='[PV8] FitProductionVertex 新重载；examples/hypertriton_by/core_bridge.cxx。')

page('带电粒子的 S 必须按弧长理解',
    P('在均匀磁场、忽略材料能损的模型中，磁场只改变动量方向，不改变模长。因此 S 仍等于有符号三维弧长除以动量模长。记衰变端为 D，产生端为 P：'),
    eq(r'S=L_{\rm arc}/p,\qquad \kappa=0.00299792458\,qB_y,\qquad\theta=\kappa S',15),
    eq(r'\mathbf r_P=\mathbf r_D-A(S)\mathbf p_D,\qquad\mathbf p_P=R(S)\mathbf p_D',15),
    ChargedArcMatrix(),
    eq(r'R(S)=\frac{\partial A}{\partial S}'),
    P('S 的单位是 cm/(GeV/c)。这里的 θ 已经等于 κS，不应再除一次 p；前文以路程 s 为自变量的子track输运才写成 θ=κs/p。'),
    P('κ→0 时，A(S)→SI、R(S)→I，连续回到中性／零场结果。实现使用小角度展开计算 sinθ/θ 及 (1−cosθ)/θ，避免小磁场下两个近似相等数相减造成精度损失。'),
    src='[PV8] HelixMatrices；符号遵循 q(v×B)，正电荷向 +z 运动且 By>0 时向 −x 弯转。')

page('带电 PV 约束怎样同时求解 S',
    P('输入仍是十维联合状态 X=(SV,pD,E,PV) 和完整协方差 C₀。S 没有独立测量误差，也没有先填一个标准寿命；它是与其他参数一起求出的自由飞行参数。'),
    eq(r'\mathbf g(\mathbf X,S)=\mathbf r_D-\mathbf r_{\rm PV}-A(S)\mathbf p_D=0'),
    P('每一轮在当前估计附近线性化。记 H 为 g 对十维 X 的导数，v 为螺旋线位移对 S 的导数，则：'),
    eq(r'H=(I,-A,0,-I),\qquad\mathbf v=R(S)\mathbf p_D'),
    eq(r'\mathbf g+H\,\Delta\mathbf X-\mathbf v\,\Delta S=0'),
    P('沿 v 的变化可以通过调整 ΔS 吸收；垂直 v 的两个方向才约束母粒子与 PV 是否相容。程序将这一轮线性方程投影到两个横向方向，求出状态更新，再由剩下的一条方程求 ΔS。'),
    P('因此三个位置方程只增加两个自由度。这个做法不需要给 S 一个虚构的初始方差，也不要求整个十维协方差可逆；例如精确已知的 PV 可以采用零协方差。'),
    note('迭代不等于重复测量','每轮都重新使用同一份 X₀、C₀。它们只提供一次测量信息；迭代的作用是改进非线性螺旋线关系的线性化。直线投影只用于 S 的初值，收敛结果必须满足完整螺旋线方程。'),
    src='[PV8] ChargedConstraint / FitChargedProductionVertex / Gain。')

page('八维误差如何随动量旋转一起输运',
    P('得到拟合结果后，v0.6.3 默认返回产生端。该端点的位置直接使用联合拟合中的 PV 估计，动量由衰变端旋转回来。输入事件里的 PV 对象不会被覆盖。'),
    P('设 Jx 描述原测量扰动如何改变受约束的十维状态，则 S 对原测量的导数为：'),
    eq(r'J_S=\frac{\mathbf v^T HJ_x}{\mathbf v^T\mathbf v}'),
    P('端点输运时，必须保留位置与动量对 S 的两组导数：'),
    eq(r'\frac{\partial\mathbf r_P}{\partial S}=-R(S)\mathbf p_D,\qquad \frac{\partial\mathbf p_P}{\partial S}=\frac{\partial R}{\partial S}\mathbf p_D',15),
    P('第一组说明飞行参数变动会改变落点；第二组说明飞行得更远也会多转过一点角度。只改位置、沿用旧动量，或者漏掉动量行对 S 的导数，都会让返回状态与协方差不再对应同一条track。'),
    eq(r'C=LL^T,\qquad C_{\rm new}=(JL)(JL)^T'),
    P('这里用协方差平方根的乘积传播，减少强相关项相减时的数值损失。特别是 PV 误差为零时，产生端的位置误差本来就应为零，不能靠随意添加小方差来掩盖数值问题。'),
    note('电荷和磁场必须跟随对象保存','ProductionState 除了八维参数，还保存 charge、magneticFieldByTesla 和 atProductionVertex。S 本身不能告诉程序该向哪边弯转。'),
    src='[PV8] TransportEndpoint / ValidCovariance；[Y] ProductionState。')

page('质量、弧长和寿命分别在什么时候改变',
    P('只进行无材料的磁场输运时，|p| 与 E 不变，所以质量也不变。S、L 与 cτ 表示同一个产生—衰变区间，换一个端点表示不会改变它们的数值。'),
    eq(r'm=\sqrt{E^2-p^2},\qquad L_{\rm arc}=Sp,\qquad c\tau=Sm'),
    P('给粒子增加一项 S，而不加入新条件，也不会改善质量分辨。质量可能发生变化的步骤是 PV 约束：位置残差通过输入协方差中的相关项，更新动量和能量。它没有把质量强行拉向一个指定值。'),
    P('带电时，下面三个量要分开：'),
    table([['量','实际含义'],['|SV−PV|','两端直线距离，不能直接代替弧长'],['pD·(SV−PV)/p²','以衰变端切线做的投影，可作短弧初值'],['S 与 L_arc=S|p|','由螺旋线几何拟合得到的飞行参数及有符号弧长']],[220,279]),
    Spacer(1,12),
    P('L 和 cτ 的误差仍按上一章的八维梯度传播，但此时必须使用带电拟合产生的协方差。端点输运前后，这些标量及其一阶误差应在数值精度内一致。'),
    note('比较质量谱时保持候选样本一致','先用完全相同的候选比较约束前后质量，再单独看 PV χ² 与 S/σS 选择。否则质量峰变窄可能只是部分候选被筛掉，不能全部归因于状态估计变好。'),
    src='[PV8] GetFlightObservables；test_charged_production.cxx 的标量／误差不变性检查。')

page('如何从超氚共同顶点接到新接口',
    P('先用物理动量和正确电荷构造两条子track：氦三 q=+2，π⁻ 的 q=−1。调用通用 VertexFitter 完成共同顶点和七维协方差，再将母粒子传给新接口：'),
    code('ProductionVertexInput input;\ninput.mother = fittedMother;\ninput.vertex = primaryVertex;\n\nProductionState topo;\nbool ok = FitProductionVertex(\n    input, ProductionVertexConfig(),\n    geometricChi2, geometricNdf,\n    topo, +1, +0.5);  // charge, By in tesla\n\nif (ok) {\n  FlightObservables flight;\n  GetFlightObservables(topo, flight);\n  TransportToDecayVertex(topo);\n}'),
    P('如果母粒子与 PV 有已知相关性，还要填写 input.motherVertexCovariance。示例的 tracking 测量与 PV 是独立生成的，因此该交叉块为零；实际数据不能仅凭代码默认值就认定独立。'),
    table([['新增字段','作用'],['charge / motherTopoCharge','以 e 为单位的电荷整数'],['magneticFieldByTesla / motherTopoByTesla','与状态匹配的带符号 By；ROOT 输出格式为 6']],[266,233]),
    Spacer(1,10),
    note('别把刚度当作动量','氦三测得的刚度 p/2 必须按输入定义转换成物理动量 p；相应协方差也要一起转换。不能只把 charge 设成 2，而让数组里继续装 p/2。'),
    src='examples/hypertriton_by/core_bridge.cxx；[PV8] 显式 charge、By 重载；[A] 输出字段。')

page('用同一套超氚样本核对正式核心',
    P('本页保留 v0.6.3 合入时的验证记录：使用 2026-09-29 独立测试的相同随机种子、输入和选择。v0.6.4 的逐候选对照见输运优化章节。共同顶点及母粒子误差仍来自 C++ 核心；新增的带电 PV 步骤也直接调用正式 C++。原 Python 算法仅用于数值对照，不作为失败后的替代结果。'),
    table([['验证','本次结果'],['电荷／场／S 正负号的独立对照','16 组：八参数 GLS 与完整协方差一致；洛伦兹力 ODE 一致'],['C++ 两轨核心输出的 GLS 对照','32 个超氚信号：状态差小于 3×10⁻⁸ 个标准差'],['全部通过共同选择的组合','97,974 个；两实现同为 97,973 个有效、1 个无效'],['全样本 C++ 与 Python 对照','最大状态差小于 8×10⁻⁸ 个标准差；协方差尺度化差小于 2×10⁻⁸'],['同一批 5,410 个信号的质量 RMS','约束前 2.692，约束后 2.650 MeV/c²；复现原独立测试']],[238,261]),
    Spacer(1,11),
    P('这里的输入包括 By=+0.5 T、母电荷 +1、子电荷 +2/−1。TOF 只用于时间一致性，动量由 tracking 提供；没有施加超氚质量或寿命约束。输入谱形和分辨率取自原玩具测试设定，不能当作 CEE 实测性能。'),
    note('原有长尾并没有被这次合入消除','S pull 的中央宽度接近 1，但全样本标准差约 1.19，仍有 9 个信号的 |pull|>5；疑似子track顶点分支问题尚待确认。本次验证的是带电核心的数值一致性，并未完成探测器级闭合。'),
    src='examples/hypertriton_by/output/{formal_vs_python,validation,summary}.json；MIGRATION_0_6_3.md。')

page('带电接口的边界与失败状态',
    P('本版完成的是均匀 By 下的带电母粒子状态、PV 约束和端点输运。默认 V0 Finder 与 ROOT 分析仍沿用原物种配置；完整超氚 PID、氦三输入适配、Ξ/Ω 组合及级联反向重拟合不在本次范围内。'),
    table([['项目','本版约定'],['旧接口只传电荷、不传 By','非零电荷返回状态 11 MissingMagneticField'],['实际零场','必须显式传 charge, 0.0；回到原中性求解分支'],['局部带电飞行范围','默认 |L_arc|≤1000 cm、|κS|<π'],['超出范围','状态 12 FlightRangeExceeded；不猜测绕行圈数'],['原状态 9 UnsupportedCharge','保留编号供旧结果识别；新带电入口不再用它拒绝电荷'],['材料、非均匀磁场、磁场误差','尚未纳入此端点与约束模型']],[234,265]),
    Spacer(1,12),
    P('这些范围限制不等于已经解决所有局部解问题。候选仍可能落在错误的子track交点分支；需要结合探测器的实际路程、顶点分布和模拟 pull 继续检查。无效结果会保留状态码与 NaN。'),
    P('原中性 V0 的直线 pointing、母粒子 DCA 等选择，不能直接套到带电超氚上。包内超氚示例使用共同顶点选择，再用带电 PV χ² 与 S/σS；没有经过旧的直线母粒子选择。'),
    note('更新后需要重新编译','ProductionState 新增了电荷和磁场字段，不能混用旧对象文件。按 MIGRATION_0_6_3.md 清理并重编两个模块。本环境尚未验证实际 CEE/ROOT/FairRoot 链接及 ROOT 文件读写。'),
    src='[PV8] ProductionVertexConfig / FitProductionVertex；[Y] ProductionFitStatus。')
