page('零场固定路程的完整 6×6 导数矩阵',
    P('上一页写出了 r′=r+su，其中 u=p/p、p=|p|。现在固定路程 s，考察六个输入分量变化时六个输出怎样改变。位置对动量的导数来自单位方向 u 对动量的变化。'),
    eq(r'\frac{\partial u_i}{\partial p_j}=\frac{\delta_{ij}-u_i u_j}{p}'),
    FixedPathMatrix(),
    P('行变量是 (x′,y′,z′,p′<sub>x</sub>,p′<sub>y</sub>,p′<sub>z</sub>)，列变量是未输运的六维状态。左上和右下单位阵表示平移前后位置与动量的直接对应；右上块表示方向误差经过距离 s 后造成的位置误差。'),
    P('这个矩阵适用于 CEE 的固定弧长参数。若 s 本身是由输入求出的最近路程，位置导数还包含 u·∂s/∂X 的贡献。v0.6.4 保留 TransportToClosestApproach 的原中心差分：每次扰动后重新求路径，将路径导数纳入传播。解析导数优化只用于固定路径 PropagateState。'),
    src='[T] PropagateParameters / PropagateState / TransportToClosestApproach；本页为零场解析对照。')

page('均匀 By 输运的位置公式与适用范围',
    P('沿用上一页的 ω 和 θ=ωs，积分动量方向可以得到三个位置分量。式中的 pₓ、pᵧ、p<sub>z</sub> 均指初始参考点的动量分量，p 是不变的动量大小。'),
    eq(r'x(s)=x_0+\frac{p_x\sin\theta-p_z(1-\cos\theta)}{p\omega}',18),
    eq(r'y(s)=y_0+s\frac{p_y}{p}',18),
    eq(r'z(s)=z_0+\frac{p_z\sin\theta+p_x(1-\cos\theta)}{p\omega}',18),
    P('当 ω 趋近零时，sin(ωs)/ω 趋近 s，(1−cos(ωs))/ω 趋近零，于是回到 r(s)=r₀+su。实际实现直接为极小旋转率使用直线分支，避免在这个极限做不稳定的除法。'),
    P('曲线公式描述一个给定的均匀 By 模型。寻找最近点时，搜索也受 maximumAbsPathCm 限制；不能把区间内的最小距离结果推广成任意场、任意材料和无限路径下的全局解。'),
    P('在 x-z 平面绕回原方位，并不等于三维位置也相同：若 p<sub>y</sub> 非零，每绕一圈 y 都会改变。因此，最近点搜索必须同时检查沿场的 y 距离。输入输出坐标及协方差顺序仍为 (x,y,z,px,py,pz)，没有交换坐标轴。'),
    src='[T] PropagateParameters；FindClosestApproachParameters。')

page('用一张图理解 Newton 的下一步',
    P('Newton 方法在当前位置用切线近似函数，并把切线与横轴的交点作为下一步。下面选 g(s)=s³−2 作为纯数学示例，让图中的横坐标与切线公式使用同一比例。'),
    NewtonGraph(),
    P('取 s<sub>n</sub>=1.6，g=2.096，g′=7.68，则下一步约为 1.3271，已经比起点更接近真根 1.2599。但这次成功不能证明任意函数或任意起点都能如此收敛。'),
    P('对实际track距离问题，g(s)=(r(s)−v)·u(s)，而 g′(s)=1+(r(s)−v)·u′(s)。导数太小会使一步跳得很远，所以代码要求 Newton 的提议保留在已括住根的区间内；否则取区间中点。'),
    P('螺旋线在 x-z 平面的周期为 T=2π/|ω|。相同平面方位的候选中，沿场 y 距离应尽量小。代码据此将搜索约化到 y 方向最优附近至多两个周期，再用 g′ 的根划分单调区间。保护 Newton 只负责区间内求根，最终仍要比较距离。'),
    src='[T] EvaluatePath / RefineStationaryPath / FindClosestApproachParameters。')

page('直线顶点初值与各层终止条件',
    P('顶点拟合开始前需要一个位置初值。取两条参考直线 r₁+s₁u 与 r₂+s₂v，令 w=r₁−r₂，并定义 b=u·v、d=u·w、e=v·w。令两条直线的间隔分别垂直 u 和 v，得到：'),
    eq(r's_1=\frac{be-d}{1-b^2},\qquad s_2=\frac{e-bd}{1-b^2}',18),
    eq(r'\mathbf v_0=\frac{(\mathbf r_1+s_1\mathbf u)+(\mathbf r_2+s_2\mathbf v)}{2}',16),
    P('当 1−b² 接近零时，两条直线近乎平行，初值会对测量波动非常敏感。程序使用最小开角进行保护；这个几何种子随后还要经过输运和带误差的顶点求解。'),
    table([['层次','核心构造默认值'],['开角保护','最小开角模 π 为 10⁻⁴ rad'],['单根细化／POCA','最多 12 次；驻点残差容差 10⁻⁵ cm；最大 |s|=1000 cm'],['共同顶点迭代','最多 12 轮；顶点变化容差 10⁻⁴ cm'],['测量方差下限 ε','10⁻⁶ cm²；共同法向残差使用 2ε'],['TOF 外层','最多 10 轮；顶点变化 10⁻³ cm、时刻变化 10⁻⁴ ns']],[149,350]),
    Spacer(1,8),
    P('这些是不同函数的终止与稳定性设置。initialVertexVarianceCm2 仅为旧接口兼容保留，当前对称批量顶点算法不用它建立一份新的顶点先验。'),
    src='[F] RawLineSeed；[V] ClosestLineSeed、VertexFitterConfig；[T] TransportConfig。')

page('局部 7×12 投影矩阵的每一行',
    P('把拟合后的输入按 (r₁,p₁,r₂,p₂) 排列，把母粒子输出按 (SV,P,E) 排列。上一页的局部投影矩阵记为 B，其全部元素如下；这里 pᵢ 和 Eᵢ 指共同约束后的子粒子值。'),
    LabeledGrid('projection'),
    P('前三行取两条已更新位置的平均，所以对应位置列为 1/2。中间三行对动量求和，所以对应动量列为 1。能量行中只有动量列非零，每项为相应分量 p<sub>ji</sub>/E<sub>i</sub>。其余格子为零。'),
    eq(r'(7\times12)\,(12\times12)\,(12\times7)\ \longrightarrow\ (7\times7)',15),
    P('B C<sub>f</sub>Bᵀ 会使用 C<sub>f</sub> 中共同约束产生的跨track交叉项。如果只取两条子粒子各自的动量协方差再相加，会漏掉这些项。最终输出采用后文重新运行顶点估计器得到的 J；两种矩阵相同的尺寸不意味着相同的输入映射。'),
    src='[V] ConstructMotherAtCommonVertex 中 projection[7][12]；B 是本书记号。')

page('TOF 的单track误差怎样组成六维输入',
    P('由 β=ℓ/(cΔt) 求动量时，代码先把路程误差 σ<sub>ℓ</sub> 和该条 TOF 测量的时间误差 σ<sub>t</sub> 视为独立。这里的 σ<sub>t</sub> 只表示单个停止探测器的测量误差，公共 T0 另行处理。'),
    eq(r'\sigma_\beta^2=\left(\frac{\sigma_\ell}{c\Delta t}\right)^2+\left(\frac{\ell\sigma_t}{c\Delta t^2}\right)^2',18),
    eq(r'\sigma_p=m\gamma^3\sigma_\beta,\qquad\gamma=(1-\beta^2)^{-1/2}',17),
    P('β 接近 1 时，γ³ 会放大一个很小的速度误差。这解释了为什么程序需要检查 β 的有效范围，也解释了高速度区单凭 TOF 得到的动量误差可能迅速变大。'),
    P('方向误差由一个小角度标准差 σ<sub>θ</sub> 近似，动量协方差分成平行方向和两个横向方向：'),
    eq(r'C_p=\sigma_p^2\mathbf u\mathbf u^T+p^2\sigma_\theta^2(I-\mathbf u\mathbf u^T)',17),
    P('沿track的变化主要来自动量大小误差，垂直track的变化主要来自方向误差。当前 TOF 状态以位置块 C<sub>r</sub> 和动量块 C<sub>p</sub> 构造六维协方差，初始位置–动量交叉块设为零。这是输入近似；后续输运和共同顶点仍会产生关联。完整跟踪输入则可以保留上游给出的六维交叉项。'),
    src='[Y] EstimateTofMomentum / BuildTofStateCovariance；[F] BuildTrackState。')

page('完整 7×12 Jacobian 的输入列与输出行',
    P('本矩阵每一格 J<sub>ij</sub> 是第 i 个母粒子输出对第 j 个拟合器输入分量的导数，行列编号从 0 开始。列顺序也就是 24 次扰动诊断掩码所使用的顺序。'),
    LabeledGrid('jacobian'),
    eq(r'J_{ij}=\frac{\partial F_i}{\partial X_j},\qquad i=0,\ldots,6,\quad j=0,\ldots,11',16),
    P('例如第 0 行第 9 列表示 SV 的 x 坐标对第二条输入track pₓ 分量的敏感程度。求这一项时，不是直接读取一个固定几何公式，而是改变该输入后重新运行顶点拟合，再比较母粒子位置。'),
    P('每格的单位是“输出变量单位／输入变量单位”。位置对动量求导可以具有 cm/(GeV/c) 的单位；把 J C<sub>12</sub>Jᵀ 完整相乘后，母粒子协方差各格重新具有对应两个输出变量的乘积单位。'),
    P('输入的分块协方差、外层 J、最终 C<sub>7</sub> 分别对应输入误差、算法响应和输出误差。明确这三个对象，才不会把导数矩阵误当协方差，或把一个输入分量误当成一条新的测量。'),
    src='[V] RecalculateMotherCovariance；列 0-5 为第一条track，列 6-11 为第二条。')
