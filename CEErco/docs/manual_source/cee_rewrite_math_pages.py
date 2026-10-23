# Core algorithm pages for the Chinese CEE v0.6.1 manual.
# The host builder supplies page, P, sub, eq, note, table and code.

page('先把track移到适合比较的位置',
 P('两条track的参考点通常位于探测器内，距离衰变顶点并不近。直接比较这两个参考点，没有明确的顶点意义。因此，程序首先沿各自的track模型移动参考点，把状态带到待考察顶点附近。这个操作称为<b>输运</b>。输运本身没有加入新测量，也不等于顶点约束。'),
 P('当前零磁场情形最容易理解：忽略材料作用时，粒子沿直线运动，动量保持不变。设初始位置为 r，单位动量方向为 u，沿track移动的有符号路程为 s，则'),
 eq(r'\mathbf{r}(s)=\mathbf{r}+s\mathbf{u},\qquad \mathbf{u}=\mathbf{p}/|\mathbf{p}|',size=18),
 P('正的 s 表示沿动量方向向前走，负的 s 表示反向回溯。给定一个待考察的空间点 v，最近位置由直线投影直接得到：'),
 eq(r's_{\mathrm{near}}=(\mathbf{v}-\mathbf{r})\cdot\mathbf{u}',size=18),
 P('这里不需要 Newton 迭代。程序还检查路程是否超过允许范围，以及返回点是否满足最近点条件。位置改变时，位置与动量之间的相关误差也要随之传播；不能只改坐标而原样保留协方差。'),
 note('本书的路程记号','这里的 s 对应 signedPathCm，单位为 cm。它不是 STAR 粒子状态中的第八维 S，也不是 STAR 接口所用的路径除以动量形式的 dS。'),
 src='[T] PropagateParameters；FindClosestApproachParameters；TransportToClosestApproach')

page('有磁场时，track为什么变成螺旋线',
 P('CEE 的束流沿 z 轴，磁场沿 y 轴。本修正版采用 B=(0,By,0)：带电粒子的动量在 x-z 平面内转向，沿磁场的 p<sub>y</sub> 保持不变。因此，x-z 投影是圆，三维track是绕 y 轴的螺旋线。By 的正负分别表示 +y 和 -y。'),
 P('程序用电荷 q、磁场 By 和动量大小 p 确定旋转率。采用 cm、GeV/c 和 T 为单位时，定义'),
 eq(r'\omega=0.00299792458\,qB_y/p,\qquad \theta=\omega s',size=18),
 P('在代码的洛伦兹力符号约定下，动量按下面的方向旋转。这个符号关系应与输入电荷、磁场方向一起检查，不能只凭track“向左弯还是向右弯”的直觉判断。'),
 eq(r'p_x(s)=p_x\cos\theta-p_z\sin\theta',size=18),
 eq(r'p_z(s)=p_z\cos\theta+p_x\sin\theta,\qquad p_y(s)=p_y',size=18),
 P('正电荷沿 +z 运动、By 为正时，应向 -x 弯转，这与 q(v×B) 一致。位置沿路程积分，例如 y(s)=y(0)+s p<sub>y</sub>/p。磁场不改变动量大小；本模型未加入平均能损或散射。零场或中性粒子回到直线形式。'),
 note('当前使用边界','有场运行默认要求已测得电荷，TOF 动量默认不允许用于有场情形。当前 TOF 路程仍按直线投影构造，不能因为螺旋线输运已经存在，就认为螺旋线路程的 TOF 动量求解也已经完成。'),
 src='[T] PropagateParameters；[F] ChargesAllow；BuildTrackState')

page('最近点不是顶点，而是拟合的准备',
 P('track 上离给定点最近的位置称为 POCA（Point of Closest Approach）。给定点可以是当前估计的 SV，也可以是已知 PV。求出 POCA 只说明“沿这条track走到哪里离目标最近”，并不说明粒子确实经过该目标，更没有强制它经过目标。'),
 P('把目标记为 v，track位置记为 r(s)，距离平方为 D²(s)。最近点处，距离对路程的一阶变化为零：'),
 eq(r'D^2(s)=|\mathbf{r}(s)-\mathbf{v}|^2',size=18),
 eq(r'g(s)=\frac{1}{2}\frac{dD^2}{ds}=(\mathbf{r}(s)-\mathbf{v})\cdot\mathbf{u}(s)=0',size=17),
 P('物理图像很直接：连接目标与最近点的线段，应当垂直于track当地的切线。直线情形可以解析求解；螺旋线情形一般需要数值求根。但驻点不一定是最小距离点，也可能是最大距离点，所以求出 g=0 后还要继续判断。'),
 P('为了传播协方差，程序会逐分量扰动输入位置或动量，并重新寻找最近点。这会把路径随输入的变化包含在数值导数中。v0.6.4 保留这一过程，只优化协方差矩阵乘法。不能直接使用固定路径导数代替。'),
 note('两个输运函数的区别','PropagateState 对指定的固定路程传播；TransportToClosestApproach 对重新求出的最近路程传播。后者的 6×6 导数矩阵包含最近路程对六维输入状态的变化。'),
 src='[T] EvaluatePath；PropagateState；TransportToClosestApproach')

page('Newton 为什么需要二分法保护',
 P('螺旋线可能绕过不止一圈，因此距离函数可能出现多个驻点。从一个起点直接做 Newton 迭代，即使收敛，也未必找到允许范围内最近的位置。当前实现先确定需要搜索的有限区间，再分段找根，最后比较候选点的距离。'),
 P('程序利用 x-z 平面运动的周期性，以及沿磁场 y 方向最接近目标的位置，缩小搜索区间。这里沿场方向不是束流 z 方向。随后用 g 的导数零点划分单调段，在有符号变化的区间中细化根。Newton 提出的下一步为'),
 eq(r's_{\mathrm{new}}=s-\frac{g(s)}{g^\prime(s)}',size=20),
 P('如果导数过小、结果不是有限数，或者新位置跑出已知根所在的区间，代码就改用区间中点。这是二分法发挥作用的地方：它保证下一步仍留在可控范围内。每走一步，程序重新缩小括住根的区间。'),
 P('找到驻点后，还要排除极大值；导数接近退化时，会比较左右邻域的距离。保留下来的局部最小值再与路径边界比较，选出距离最小者。因此，这里的稳健性来自分段搜索、受保护的求根和最终比较共同作用。'),
 note('不要把几个次数混在一起','每个根默认最多细化 12 次；顶点内部默认最多 12 轮；TOF 外层默认最多 10 轮；母粒子协方差的 24 次则是输入正负扰动的重拟合次数。它们解决四个不同的问题。'),
 src='[T] RefineStationaryPath；IsLocalMinimum；FindClosestApproachParameters；TransportConfig')

page('两条track怎样共同确定一个顶点',
 P('现在可以回到顶点问题。程序先给出一个 SV 估计，再分别把两条原始输入状态输运到离这个估计最近的位置。两条track使用同一个冻结的顶点估计，随后同时参与本轮更新。这样不会因为先处理质子还是先处理 pion 而改变更新顺序。'),
 P('一条track沿自身方向可以前后移动，真正约束顶点的是两个横向方向。为第 i 条track选两个垂直动量的单位向量，把它们排成 H_i，得到一个 2×3 投影矩阵。其二维测量协方差为'),
 eq(r'R_i=H_iC_{r_i}H_i^T+\epsilon I_2',size=18),
 P('程序将两条track的横向信息累积起来。误差较小的方向贡献较大的权重；误差较大的方向则允许顶点有更大变化。在线性化的这一轮中，'),
 eq(r'A=\sum_i H_i^TR_i^{-1}H_i,\qquad \mathbf{b}=\sum_iH_i^TR_i^{-1}H_i\mathbf{r}_i',size=16),
 eq(r'\mathbf{v}=A^{-1}\mathbf{b}',size=20),
 P('A 是 3×3 信息矩阵。得到新顶点后，程序再从原始输入重新输运和线性化，直到顶点变化足够小。这里是在求同一个问题的稳定解，并非把同一份测量重复加入、不断压小误差。'),
 src='[V] BatchUpdateVertex；VertexFitter::Fit')

page('为什么两条track的几何拟合只有一个自由度',
 P('看到三维空间中的两条track，人们很容易以为必须独立检验三个位置分量是否相等。但沿每条track的位置不是固定的：求共同顶点时，可以分别调整两条track上的路程。这两个可调路程会消去两个表面上的位置差。'),
 P('另一种计数更直观。每条track提供两个横向位置测量，两条track一共提供四个；待求顶点有三个坐标。因此，两体顶点拟合最终剩下一个有效自由度：'),
 eq(r'\mathrm{NDF}_{\mathrm{geom}}=2\times2-3=1',size=20),
 P('在局部直线近似下，两条track都不能通过沿自身方向移动来消除的差异，正是公垂线方向上的间隔。代码用两条单位方向的叉积构造法向量 n，再取这个方向上的位置差作为残差：'),
 eq(r'\mathbf{n}=\frac{\mathbf{u}_1\times\mathbf{u}_2}{|\mathbf{u}_1\times\mathbf{u}_2|},\qquad z=\mathbf{n}\cdot(\mathbf{r}_2-\mathbf{r}_1)',size=17),
 P('近乎平行或反平行时，叉积太小，共同顶点沿某些方向难以确定。程序在初始化时检查这一病态几何，而不是等矩阵求逆失败以后再补救。这里的一维残差和三维顶点并不矛盾：前者是剩余不相容程度，后者是拟合参数。'),
 src='[V] ClosestLineSeed；ConstructMotherAtCommonVertex；VertexFitter::Fit')

page('共同顶点约束怎样调整两条 daughter',
 P('顶点位置稳定以后，程序把两条已输运的六维状态合在一起，形成十二维向量 X=(r1,p1,r2,p2)。此时的输入协方差 C_X 是 12×12：两条track各占一个 6×6 对角块，尚未加入它们因共同顶点约束而产生的相关性。'),
 P('上一页的公垂线残差可以写成 z=hX。h 是 1×12 的行向量，位置部分分别为 −n 和 n，动量部分为零。残差的方差 V 是标量，Kalman 增益 K 则是 12×1：'),
 eq(r'h=(-n_x,-n_y,-n_z,0,0,0,n_x,n_y,n_z,0,0,0)',size=15),
 eq(r'V=hC_Xh^T+2\epsilon,\qquad K=C_Xh^T/V',size=18),
 eq(r'X_f=X-Kz,\qquad C_f=C_X-\frac{C_Xh^ThC_X}{V}',size=17),
 P('虽然残差直接比较的是位置，更新也可能改变动量。原因是输运后的协方差已经包含位置与动量的相关项：某个位置偏差，可能通过轻微改变方向来解释。与此同时，两条 daughter 之间会出现新的交叉协方差，因为它们现在共享同一个顶点条件。'),
 eq(r'\chi^2_{\mathrm{geom}}=z^2/V',size=20),
 P('这个 χ² 描述两条track形成共同顶点是否相容。程序在约束更新前记录两条 POCA 位置的间距，作为 daughterDcaCm；它不会先把track调到一起，再把得到的零间距当成筛选变量。'),
 src='[V] ConstructMotherAtCommonVertex；VertexFitter::Fit')

page('从两条 daughter 得到七维母粒子',
 P('共同顶点约束之后，母粒子的动量是两条 daughter 动量之和，能量是按当前质量假设计算的能量之和。母粒子不是再独立测量一次，而是由这两条相关的 daughter 状态构造出来：'),
 eq(r'\mathbf{P}=\mathbf{p}_{1f}+\mathbf{p}_{2f},\qquad E=E_{1f}+E_{2f}',size=18),
 eq(r'E_{if}=\sqrt{|\mathbf{p}_{if}|^2+m_i^2}',size=20),
 P('当前实现将两条已更新位置的平均值定义为 mother.positionCm，从而保持对 daughter 顺序的对称性。这个平均发生在输运、共同顶点求解和约束更新之后，不能据此把整个算法概括成“找两个最近点，再取中点”。'),
 P('最终中心状态的顺序固定为 Y=(x,y,z,Px,Py,Pz,E)。从十二维拟合状态到这七个量的局部导数矩阵记为 B，维度是 7×12：位置对相应 daughter 位置的导数为 1/2，总动量对相应动量的导数为 1，能量对各 daughter 动量的导数为 p_i/E_i。'),
 eq(r'C_{7,\mathrm{local}}=B C_f B^T',size=20),
 P('这先得到一张内部使用的七维协方差。进入最终候选输出前，Finder 还会通过输入扰动重新计算母粒子协方差，以包含最近点和顶点重新求解带来的变化。后面将单独解释这一步。'),
 src='[V] ConstructMotherAtCommonVertex；[F] BuildCandidate')

page('TOF 动量为什么离不开衰变顶点',
 P('零场直线跟踪可以较好地给出位置和方向，却不能靠曲率确定动量大小。当前 TOF 模式因此利用飞行路程与飞行时间估计速度，再结合粒子质量假设换算动量。这里的质量假设不可省略：同一个速度，作为质子和作为 pion 时，对应的动量不同。'),
 P('TOF 测量给出的是粒子到达读出面的时刻。daughter 并不是在碰撞发生时就从 PV 出发，而是在母粒子到达 SV、发生衰变以后才出发。因此应扣掉的起始时刻是 t_decay，路程也应从 SV 算起：'),
 eq(r'\ell_i=\mathbf{u}_i\cdot(\mathbf{r}_{\mathrm{TOF},i}-\mathbf{SV}),\qquad \Delta t_i=t_{\mathrm{TOF},i}-t_{\mathrm{decay}}',size=16),
 eq(r'\beta_i=\frac{\ell_i}{c\Delta t_i},\qquad p_i=\frac{m_i\beta_i}{\sqrt{1-\beta_i^2}}',size=19),
 P('方向 u_i 来自输入track，并在有 TOF 信息时定向为朝向 TOF 匹配位置。代码使用的是沿该方向的路程投影，要求路程和时间都为正，并要求 β 小于配置上限；默认上限为 0.9995。'),
 P('问题也由此出现：算 daughter 动量需要 SV 和衰变时刻；求 SV、母粒子速度和衰变时刻又需要 daughter 动量。接下来的外层迭代，就是让这些相互依赖的量逐步达到一致。'),
 src='[F] OrientedDirection；BuildTrackState；[Y] EstimateTofMomentum')

page('第一次迭代从哪里开始',
 P('程序尚不知道 daughter 动量大小时，仍然知道参考位置和track方向。因此，先把两条输入方向看作直线，求两条直线最近位置的中点，作为 SV 的初始估计。这个种子只负责启动计算，不直接决定最终重建顶点。'),
 P('设两条单位方向为 u、v，参考位置差为 w=r1−r2。两条直线上各自的最近路程由几何关系求出，再取对应位置的中点。若两条方向近乎平行或反平行，交会位置对微小测量误差极敏感，程序会拒绝这种初始几何。'),
 P('有了 SV 种子以后，还需要一个母粒子到达这里的时间估计。代码使用'),
 eq(r't_{\mathrm{seed}}=t_0+\frac{|\mathbf{SV}_{\mathrm{seed}}-\mathbf{PV}|}{\beta_{\mathrm{seed}}c}',size=19),
 P('默认 β_seed=0.7。它只是尚无母粒子动量时采用的启动值；进入迭代后，将被当前重建的母粒子速度所替代。因此，“程序假设所有 Lambda 的速度为 0.7”是错误理解。'),
 note('种子与物理约束的区别','初值是数值求解的出发点；物理约束会在最终解中持续起作用。当前 seedParentBeta 属于前者。实际使用时应检查收敛情况，而不能只因为初值合理，就认为所有候选都已得到可靠解。'),
 src='[F] RawLineSeed；SolvePair；FinderConfig')

page('TOF 外层怎样形成一个自洽解',
 P('给定本轮 SV 和衰变时刻后，程序先为两条 daughter 构造六维状态：方向沿用跟踪结果，使用 TOF 的track按本轮路程和时间求动量。随后调用前面已经解释过的共同顶点拟合，得到新的 SV、母粒子动量 P 和能量 E。'),
 P('由此可以更新母粒子在实验室系中的速度。母粒子从 PV 到 SV 的时间采用沿动量方向的有符号位移投影：'),
 eq(r'\beta_M=\frac{|\mathbf{P}|}{E},\qquad L_{\parallel}=\frac{\mathbf{P}}{|\mathbf{P}|}\cdot(\mathbf{SV}-\mathbf{PV})',size=17),
 eq(r't_{\mathrm{decay}}=t_0+\frac{L_{\parallel}}{\beta_Mc}',size=20),
 P('新的 SV 和 t_decay 又改变下一轮 daughter 的 TOF 路程与时间。程序反复进行“构造 daughter→共同顶点拟合→更新时间”，直到新旧 SV 和时间都足够接近。这属于不动点迭代，也称固定点迭代：将当前估计送入整个流程后，返回的估计几乎不再变化。'),
 P('这里没有把母粒子的标准寿命或标准质量拿来约束结果。β_M 来自本轮重建的 P/E。程序还要求 β_M 有效且 L_parallel 非负；如果 L_parallel 为负，SV 就落在以 PV 为起点、沿母粒子动量向前飞行的反方向，当前模型无法给出非负飞行时间。横向是否相容，随后另用拓扑量检验。'),
 note('一个容易遗漏的区别','初始化时间用的是 |SV_seed−PV|；后续时间更新用 L_parallel。输出的 decayLengthCm 也仍是三维距离，而不是这里的有符号投影。'),
 src='[F] SolvePair；UpdatedDecayTime')

page('怎样判断 TOF 求解真正收敛',
 P('收敛不能只看质量峰或动量是否变化很小。当前实现同时比较两件事：新顶点相对本轮输入顶点的位移，以及更新后的衰变时刻变化。两者都低于阈值，才宣布本轮达到收敛条件。'),
 table([['TOF 外层配置','构造函数默认值'],['最大迭代轮数','10'],['顶点变化阈值','10⁻³ cm'],['衰变时刻变化阈值','10⁻⁴ ns'],['要求收敛','开启']],[305,194]),
 P('达到条件以后，程序还会用最后一次更新的 SV 和时间重新构造 daughter，并额外拟合一次。这样可以避免输出的 mother 对应上一轮输入，而 decayTimeNs 已经对应下一轮。若原先声明收敛，额外拟合结果也必须继续满足相同阈值。'),
 P('输出 tofIterations 记录外层循环次数，不包含这次最后的一致性拟合。若两条track都使用 tracking 动量，则只做一次顶点拟合，直接将 tofConverged 置为真；这个标记此时表示没有需要求解的 TOF 闭环，不能理解成实际进行过多轮 TOF 检验。'),
 P('最大轮数与数值阈值由配置结构设置，当前二进制没有对应的直接命令行选项；诊断开关只能放宽是否要求 TOF 收敛。应以输出元数据中实际生效的配置为准。'),
 src='[F] SolvePair；FinderConfig')

page('为什么母粒子协方差还要重新计算',
 P('构造母粒子时已经得到一张局部七维协方差，为什么还要增加计算？关键在于，母粒子状态不仅取决于已输运的两条状态，还取决于“最近位置在哪里”和“共同顶点最后在哪里”。输入轻微改变时，这些求解结果也会一起改变。'),
 P('局部投影矩阵 B 描述：在本轮拟合状态附近，daughter 参数怎样组合成母粒子参数。最终的外层导数 J 描述的则是：从原来的两条输入六维状态出发，重新运行整个顶点估计器以后，母粒子怎样改变。两者维度同为 7×12，含义却不同。'),
 eq(r'Y=F(X_0),\qquad J=\frac{\partial F}{\partial X_0}',size=20),
 P('这里 X0=(r1,p1,r2,p2) 是十二维输入，F 包含最近点求解、顶点重新线性化、共同约束和母粒子构造。最终用 J 将输入误差传播到输出：'),
 eq(r'C_7=J C_{12}J^T',size=21),
 P('这张重新计算的 C7 会替换先前内部得到的母粒子协方差，不能再把两者相加。两份结果描述的是同一组输入误差，经不同近似层次传播后的结果，并不是两个独立误差来源。'),
 note('“完整”必须说明范围','这里重跑的是 VertexFitter。它从最终构造的两条六维状态出发，并不对所有原始 TOF 和 PV 测量重新执行整套外层闭环。'),
 src='[V] RecalculateMotherCovariance；[F] BuildCandidate')

page('24 次扰动怎样组成 7×12 导数矩阵',
 P('输入 X0 有十二个分量。程序依次只改变其中一个分量：先增加一个很小的步长，再减少同样步长，两次都重新调用顶点拟合。七个母粒子输出之差除以两倍步长，就得到 J 的一整列：'),
 eq(r'J_{:j}=\frac{F(X_0+h_j e_j)-F(X_0-h_j e_j)}{2h_j}',size=19),
 P('十二列各需要正负两个拟合，所以一共是 24 次额外重拟合。这不是随机抽样，也不是让同一候选重复参与统计；它是中心差分求导。每次重拟合内部仍会求最近点，并可能执行多轮顶点更新。'),
 table([['被扰动的输入','默认差分步长'],['位置分量','10⁻³ cm'],['动量分量','max(10⁻⁶ GeV/c, |该分量| × 10⁻⁴)']],[195,304]),
 P('输入 C12 由两条 C6 放在对角块上构成，并在每条track的位置对角上加入配置的测量方差下限。程序使用 J C12 Jᵀ 得到 C7，再检查数值有限性及半正定性。'),
 P('公共事件时间误差没有在这里当成两条独立输入误差重复加入。它会在单独步骤中通过改变 T0、重跑完整 TOF 外层求解来传播。'),
 src='[V] RecalculateMotherCovariance；VertexFitterConfig；[F] AddSharedTimeUncertainty')

page('中心值成功，不等于误差计算成功',
 P('一个候选可以在原始输入上成功拟合，却在某次微小正负扰动后失败。例如，扰动可能使几何接近退化边界，或者使最近点与顶点迭代无法收敛。这种情形必须与原始中心拟合本身失败分开记录。'),
 P('当前实现会继续尝试所有列的正负扰动，记录尝试数、成功数、失败数、正负侧失败掩码，以及第一次失败的列和拟合状态。只要有一次扰动拟合失败，就不再把不完整的导数矩阵当成有效结果。母粒子协方差会被置为无效值。'),
 P('默认设置要求母 Jacobian 完整，因此这类候选不进入最终有效样本。如果显式关闭 requireCompleteMotherJacobian，可保留这一种特定失败作为 QA：中心质量、动量和已成功的几何结果仍可观察，但依赖最终母粒子协方差的误差量与拓扑量不会被伪造出来。'),
 table([['此类 QA 候选','输出含义'],['motherJacobianComplete','false'],['covarianceDependentCutsEvaluated','false'],['质量误差、长度误差及显著性','NaN'],['PV 与总拓扑 χ²','NaN']],[285,214]),
 P('这项开关不是“忽略所有失败”。输入协方差非法、过程噪声不受支持、公共时间导数失败等其他错误仍然拒绝。分析时必须显式区分有效样本和这类诊断样本。'),
 src='[V] RecalculateMotherCovariance；[F] InvalidateCovarianceDependentObservables；PassCuts')

page('为什么质量误差需要能量与动量的相关项',
 P('母粒子质量中心值由总能量和总动量计算，而不是从粒子名称读取。只要组合的质量平方为正且数值有限，就有'),
 eq(r'M=\sqrt{E^2-|\mathbf{P}|^2}',size=21),
 P('但是质量误差不能只把两条 daughter 的动量误差相加。共同顶点会关联两条 daughter，E 和 P 也来自同一组动量，因此它们之间存在不可忽略的相关性。七维母粒子协方差把这些相关项保存在同一个状态中。'),
 P('对 Y=(x,y,z,Px,Py,Pz,E)，质量的一阶变化由下面的梯度给出。前三项为零，是因为在给定 E、P 时，质量不直接依赖当前参考位置：'),
 eq(r'g_M=\left(0,0,0,-\frac{P_x}{M},-\frac{P_y}{M},-\frac{P_z}{M},\frac{E}{M}\right)',size=17),
 eq(r'\sigma_M^2=g_M C_7 g_M^T',size=21),
 P('展开二次型以后，既有能量和动量的方差，也有它们的协方差。交叉项可能增加或减小最终质量方差；这正是“相关误差”实际影响物理量的地方。'),
 note('数值检查的含义','母粒子协方差无效时，不输出可信的质量误差。方差只有在机器舍入精度允许的极小负值范围内才可归零；真正的负方差意味着误差传播存在问题，不能靠截成零来掩盖。'),
 src='[Y] MotherMassAndError；[F] BuildCandidate')

page('先把几种距离和角度分清楚',
 P('候选形成以后，程序既关心两条 daughter 是否能形成一个共同顶点，也关心该顶点与 PV 的关系。几个常见几何量名字相似，却回答不同问题。将它们分开理解，后面的 χ² 才不会失去物理含义。'),
 table([['变量','当前实现的含义'],['daughterDcaCm','共同约束更新前，两条 POCA 位置的间距'],['daughterDcaToPvCm','每条 daughter track到 PV 的最近距离'],['decayLengthCm','三维点间距离 |SV−PV|'],['pointingCosine','SV−PV 与母粒子动量的夹角余弦']],[200,299]),
 P('令 d=SV−PV，L=|d|。指向关系和沿动量的投影可以写为'),
 eq(r'\cos\theta=\frac{\mathbf{d}\cdot\mathbf{P}}{L|\mathbf{P}|},\qquad L_{\parallel}=L\cos\theta',size=18),
 P('当母粒子明确指回 PV 时，cosθ 接近 1，三维距离与投影长度接近；偏离较大时，两者可以有明显差别。当前 decayLengthCm 存 L，时间更新使用 L_parallel，不能将这两个量互换。'),
 P('衰变长度误差目前将 PV 和 SV 的位置协方差相加，再沿 d 的单位方向投影：σL²=eLᵀ(CSV+CPV)eL。这个公式包含两端位置误差，但没有 PV 与母粒子联合状态的交叉协方差。'),
 src='[F] MeasureDaughterAtPv；BuildCandidate；UpdatedDecayTime')

page('daughter 偏离 PV 的 χ² 在检验什么',
 P('弱衰变的 daughter 应当从 SV 出发，它的track通常不会很好地穿过碰撞 PV。因此，要求 daughter 与 PV 有一定程度的不相容，是一种常用的次级粒子选择。但单看厘米单位的 DCA，无法区分测量精度不同的track。'),
 P('当前代码先将每条 daughter 输入状态输运到 PV 的 POCA。然后在局部动量的两个横向方向上取位置残差，并把track位置误差与 PV 误差相加、投影到同一个二维平面。得到的二次型为'),
 eq(r'\chi^2_{\mathrm{daughter,PV}}=\mathbf{r}_\perp^T V_\perp^{-1}\mathbf{r}_\perp',size=20),
 P('同样的 DCA，如果track和 PV 的误差很小，χ² 就会较大；误差较大时，则可能仍然相容。对 daughter 的这个量通常施加下限；对母粒子是否指回 PV，则通常要求相容，也就是施加上限。两种 cut 的物理方向不同。'),
 table([['核心成员：两元素数组','保存的量'],['daughterChi2FromPv','原始 χ²'],['daughterPvSignificance','sqrt(χ²)'],['daughterDeviationFromPv','旧版兼容字段，仍是 sqrt(χ²)']],[280,219]),
 P('ROOT 分支使用 daughterPvChi2{0,1}、daughterPvSignificance{0,1}，旧兼容名为 daughterDeviation{0,1}。此处使用最终的子粒子输入状态，未加入随后在母粒子协方差层传播的公共时间贡献。'),
 src='[F] MeasureDaughterAtPv；PassCuts；[Y] Candidate 字段定义')

page('母粒子指回 PV 的二维相容性',
 P('检验母粒子是否来自 PV 时，沿母track方向的路程可以自由选择。真正有辨别力的是垂直母粒子动量的偏离，因此 CEE 用两个横向残差来描述这项相容性。这里没有先建立一个新的 PV 约束母粒子对象。'),
 P('令 d=SV−PV，取垂直 P 的两个单位向量 e1、e2，残差为 r_a=e_a·d。输入变化时，SV 会改变，母粒子动量方向也会改变。代码以局部的一阶近似同时计入这两个影响：'),
 eq(r'J_a=\left(\mathbf{e}_a^T,-\frac{L_{\parallel}}{|\mathbf{P}|}\mathbf{e}_a^T,0\right)',size=18),
 P('J 由两行组成，维度为 2×7。若 H 的两行为 e1ᵀ、e2ᵀ，残差协方差和 PV 相容性统计量为'),
 eq(r'V=JC_7J^T+HC_{\mathrm{PV}}H^T+\epsilon I_2',size=18),
 eq(r'\chi^2_{\mathrm{PV}}=\mathbf{r}^TV^{-1}\mathbf{r}',size=21),
 P('式中的动量导数项很有意义：飞得越远，相同方向误差造成的横向偏移越大。因此，只拿 SV 和 PV 的位置方差去除一个几何 DCA，并不等同于当前这个指标。'),
 note('相容性检验不等于状态更新','这里的原统计量函数返回 χ²，不修改 mother，也不创建 S。v0.6.2 的八维副本由后续独立函数生成。它采用的仍是当前 CEE 模型下的近似误差描述。'),
 src='[F] PrimaryVertexConstraintChi2')

page('三种 χ² 应该怎样读',
 P('两条 daughter 能形成一个共同顶点，并不自动意味着母粒子来自事件 PV。因此，代码将共同顶点的几何 χ² 与指回 PV 的 χ² 分别保存，再构造总拓扑指标。分析和画图时，应先确认读取的是哪一个字段。'),
 eq(r'\chi^2_{\mathrm{topo}}=\chi^2_{\mathrm{geom}}+\chi^2_{\mathrm{PV}}',size=19),
 eq(r'\mathrm{NDF}_{\mathrm{topo}}=\mathrm{NDF}_{\mathrm{geom}}+2=3',size=19),
 table([['输出字段','含义'],['geometricChi2Ndf','共同顶点 χ² / 1'],['primaryVertexConstraintChi2Ndf','仅 PV 残差 χ² / 2'],['topologicalChi2Ndf','(共同顶点 χ² + PV 残差 χ²) / 3']],[290,209]),
 P('自由度中的 2 来自去除沿母track路程后的两个横向条件。当前这样命名和计数，是为了使选择变量的用途接近 STAR 的拓扑指标；它并不说明两套算法已经执行了相同的状态更新，或在同一输入上必然给出相同数值。'),
 P('STAR 的产生顶点约束还会更新母粒子状态和相关协方差。CEE 当前只是计算相容性，再将 χ² 相加，原始母质量和 SV 不因此改变。若要升级到真正的约束副本，应作为新的算法工作验证，不能通过修改手册措辞来当作已完成。'),
 src='[F] PrimaryVertexConstraintChi2；BuildCandidate；[Y] Candidate 字段定义')

page('衰变时刻、固有寿命和 S 是三件事',
 P('当前 decayTimeNs 是实验室系中的衰变时刻：从事件时间原点 t0 出发，加上母粒子飞到 SV 所需的时间。只有在 t0=0 时，它才恰好等于实验室系的飞行时间。它不是粒子静止系中的固有寿命 τ。'),
 eq(r't_{\mathrm{decay}}-t_0=\frac{L_{\parallel}}{\beta_M c}',size=20),
 P('利用现有的质量、动量和顶点中心值，可以另行计算长度形式的固有衰变时间 cτ：'),
 eq(r'c\tau=\frac{M L_{\parallel}}{|\mathbf{P}|}=\frac{M L\cos\theta}{|\mathbf{P}|}',size=20),
 P('这说明七维状态并不妨碍得到一个 cτ 中心值。v0.6.2 增加的 topoFlight.ctauCm 来自独立的八维约束副本，并带有一阶传播的误差。它与本页由原始七维量得到的简单投影估计，要按各自状态来源区分。'),
 P('STAR 的第八维 S 保存的是产生端与衰变端之间的飞行参数，并带有它同位置、动量和能量的协方差。CEE v0.6.2 的 motherTopo 也保存 S 及这些交叉项，提供中性粒子的端点输运接口。它的构造采用 CEE 联合约束拟合，具体过程见下一章。'),
 note('与前一本 STAR 手册的衔接','可以借助 Ξ→Λπ 的两级衰变理解每一级粒子的产生端和衰变端，但不能将那套完整端点管理，直接写成 CEE v0.6.2 已完成整条级联链。本章先讨论中性粒子；v0.6.3 的带电螺旋线接口另见后文。'),
 src='[F] UpdatedDecayTime；[Y] MotherState；Candidate')

page('公共 T0 误差为什么要单独传播',
 P('两条 daughter 的 TOF 测量各有自己的时间误差，但事件 T0 的误差是公共的。如果同一个 T0 偏晚，两条 daughter 的飞行时间会一起改变，动量也会产生相关变化。把 T0 方差分别加进两张独立track协方差，会丢失这种共同变化。'),
 P('当前代码先用每条track自身的 TOF 时间误差构造输入协方差。若开启公共时间传播，再改变事件 T0，重跑整个 TOF 外层求解。由母粒子七个输出随 T0 的变化，得到一个 7×1 导数向量 a：'),
 eq(r'\mathbf{a}=\frac{\partial Y}{\partial t_0},\qquad C_7\leftarrow C_7+\mathbf{a}\sigma_{\rm common}^2\mathbf{a}^T',size=18),
 P('程序分别使用 ±h 和 ±h/2，总共重跑四次 SolvePair，比较两种步长得到的导数是否稳定。通过检查后采用半步长结果。默认 h=10⁻³ ns，一致性容差为 5%。这四次是完整 TOF 闭环求解，与前面的 24 次顶点扰动是两个独立步骤。'),
 P('当前使用 σ_common²=σT0²+additionalDecayTimeSigmaNs²，将两项都按同一种公共时间平移来传播。额外时间项因此是一种模型化的公共误差，并不表示程序已经把衰变时间作为一个独立状态参数进行联合拟合。'),
 note('默认模拟配置','默认 t0=0、σT0=0，公共时间传播关闭。分析程序切换到 real 模式时会开启公共时间误差传播；使用前仍需核对时间单位、零点及不确定度，并按标定结果修改占位误差。'),
 src='[F] AddSharedTimeUncertainty；FinderConfig；[Y] AddCommonNuisanceCovariance')

page('怎样理解当前误差模型的可信范围',
 P('一张有限、对称、半正定的协方差矩阵，是正确误差传播的必要条件，但不能单靠这些数值检查证明误差模型已经完整。还要检查输入误差是否真实，哪些共同来源已加入，哪些反馈和相关性仍在近似处理。'),
 P('未提供 geometryCovariance 的 TOF 输入把位置块与动量块分开构造；Beno 零场接口则保留前四维转换得到的位置—方向交叉项。TOF 动量大小的误差仍来自路程和单 track 时间模型。最终 24 次扰动从这些六维状态出发，重跑 POCA 与共同顶点；公共 T0 则另行重跑 TOF 闭环。这两层各有明确覆盖范围。'),
 P('PV 会参与母粒子飞行时间的计算，因此在 TOF 模式下，母粒子状态并不天然与 PV 完全独立。原七维输出的长度和 PV 相容性仍沿用分块近似。v0.6.2 对新八维副本另建含 PV 响应和交叉块的联合输入，但仍没有补齐上游未提供的相关性，也没有重做所有原始 TOF 测量的联合误差模型。'),
 P('材料方面，默认位置和角度过程噪声均为零。输运接口虽有简化噪声参数，但最终母粒子协方差重算会拒绝非零过程噪声，因为单靠输入的 J C Jᵀ 无法包含沿路新产生的随机误差。能量损失、材料分布和完整多重散射模型也不属于当前实现。'),
 note('如何检验这些近似','应在有真值的模拟样本中检查残差、pull 宽度、效率和失败比例，并观察它们对动量、衰变位置及开角的依赖。数值收敛说明求解结束了；误差是否合适，还需要这些重建性能检查。'),
 src='[Y] BuildTofStateCovariance；[V] RecalculateMotherCovariance；[F] AddSharedTimeUncertainty；BuildCandidate')
