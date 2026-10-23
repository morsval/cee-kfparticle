CEE v0.6.5 重建 macro 示例
日期：2026-10-22

1. 这个 macro 调用什么

runCeeV0Reco.C 组装配置并调用已编译的 CeeV0Analysis.bin。
真正的事件循环、TChain 读取、TOF 匹配和轻量版 KF 拟合仍在现有分析程序中。
当前正式分析入口重建 Lambda -> p + pi-，只使用 TPC-TPC 组合。
它不会自动切换到超氚 toy 的 3He + pi- 重建，也不会从原始 TPC hit 重新做 tracking。
输入须是已经保存重建 track 的 ROOT 文件。

输入默认树：ceesim。
必要分支：CeeTpcTrack、CeeTpcVertex。
零场 TOF 模式还需 CeeiTOFHit / CeeEtofHit 中至少一个。
simulation：生成器时间原点为 0，不读取 CeeT0Hit，不传播公共 T0 误差。
MCTrack 是可选真值标签输入，缺失时继续重建，不通过真值做选择。
tracking：使用完整六维种子、协方差和拟合电荷，当前不绑定 TOF hit。
real+hit：需要 CeeT0Hit，并明确提供时间单位转换。
这些分支需位于同一输入树；分散在不同文件的 track/hit/PV 不是本示例的输入形式。
所有 list 中的文件应具有相容树结构与数据类字典。

2. 放置与编译

把 runCeeV0Reco.C 放入 CEErco/CeeV0Analysis/macros/。
在已初始化 ROOT / CEE / FairRoot 环境的 shell 中运行：

  cd /path/to/CEErco/CeeV0Analysis
  make -j4

编译需要 root-config、VMCWORKDIR、FAIRROOTPATH 以及相容的 CEE 数据类和 ROOT 字典。
如果只修改本 macro 的 cut 数值，无需重新编译 CeeV0Analysis.bin。
若修改重建核心或数据类接口，则需要重新编译相应程序/库/字典。

3. 单个 ROOT 文件

  root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.root","lambda.root",10000)'

第三个参数为最多处理的事件数。-1 表示全部事件，数目是整条 chain 的总数。
输出由程序以 RECREATE 打开，同名输出会覆盖。

4. ROOT 文件列表

创建 files.list，每行一个 ROOT 文件绝对路径，例如：

  /data/cee/reco_0001.root
  /data/cee/reco_0002.root
  /data/cee/reco_0003.root

然后运行：

  root -l -b -q 'macros/runCeeV0Reco.C("/data/files.list","lambda_all.root",-1)'

支持 .list 和 .txt。允许空行与以 # 开头的整行注释。
不要给行内路径加引号，不要在路径后接注释。当前解析器只接收以 .root 结尾的行。
相对路径按启动分析程序时的当前目录解释，不按 list 文件所在目录解释。
本示例因此建议填绝对路径。程序也支持一个目录中的 .root 文件，按名称排序。
路径里的 ~ 和环境变量请先展开成实际路径。

5. 有磁场模拟

  root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.list","lambda_By.root",-1,"tracking",0.5)'

第五个参数是带正负号的 By，单位 T；CEE 的磁场沿 y，束流沿 z。
tof 要求 |By| <= 1e-6 T；tracking 要求 |By| > 1e-6 T。
这只切换已有拟合 track 的输入方式，不会为文件补出缺失的协方差或电荷。

6. 怎么改 cut

打开 macro，修改“配置和 cut”区域即可。以下只示范写法，数值不是优化工作点：

  const int minimumTpcHits = 15;
  const double maximumTrackChi2Ndf = 5.0;
  const double maximumDaughterDcaCm = 1.0;
  const double maximumGeometricChi2Ndf = 10.0;

它们分别传给 --min-tpc-hits、--max-track-chi2ndf、
--max-daughter-dca、--max-geo-chi2ndf。
candidate.daughterDcaCm > cuts.maximumDaughterDcaCm 的拒绝条件，
就是通过 maximumDaughterDcaCm 这个配置设置阈值。
candidate.geometricChi2Ndf > cuts.maximumGeometricChi2Ndf 同理。
不要修改比较运算本身来设置 cut。

其他候选参数：
  minimumMassGeV / maximumMassGeV          母粒子质量上下限，GeV/c^2
  maximumTopologicalChi2Ndf               旧几何项+旧 PV 项的总 chi2/NDF
  minimumDecayLengthCm                   约束前几何长度 |SV-PV|，cm
  minimumDecayLengthSignificance         约束前几何长度的 L/sigmaL
  minimumPointingCosine                   约束前母粒子指向余弦
  minimumDaughterChi2FromPv               两条子 track 各自到 PV 的原始 chi2 下限

此 macro 先应用 legacy，然后显式覆盖全部上面这些候选 cut。
默认上限 1e30 表示数值上放开；长度显著性下限保留 0，不加 >3 等额外选择。
仍有输入有效性、中心解、TOF 收敛、完整母粒子 Jacobian 等基础条件。
legacy 不等于“完全无选择”。
零场默认还有 proton DCA_to_PV <= 1.2 cm、pion DCA_to_PV > 1.2 cm 的角色规则。
它和两子 track 之间的 daughterDcaCm 是不同的 DCA。
可用 disableDcaRole=true 关闭角色规则；tracking 默认自动关闭它。

7. 重建时的 cut 与输出后的 cut

重建 cut：控制哪些 track 和候选进入结果，可以减少无效组合的计算。
输出后 cut：用 LambdaCandidates 中的字段选样，便于扫描阈值。
输出树只保存通过重建配置的候选；已经被前面删掉的候选无法靠后处理放宽找回来。
因此准备做 cut 扫描时，先保存覆盖扫描范围的较宽样本。
分析程序将实际运行配置写入输出的 TNamed("CeeV0AnalysisConfig")。

特别区分：
  --min-ldl 切 decayLengthSignificance，不切 motherTopoDecayLengthSignificance。
  --max-topo-chi2ndf 切 topologicalChi2Ndf，不切 motherTopoPvChi2/motherTopoPvNdf。
  primaryVertexConstraintChi2Ndf 也是旧 PV 统计项，不是新增八维拟合的自动别名。
  massGeV 是约束前母粒子质量；motherTopoMassGeV 是 PV 更新副本的质量。
  motherTopoDecayLengthCm 是带符号的 p*S，并非原始无符号 |SV-PV|。
当前 CLI 没有针对新增 motherTopo 分支的独立 cut 开关。

例：在 ROOT 交互环境中，对输出候选统一要求八维 PV chi2/NDF < 5，
然后画约束前后的谱，不加 L/sigmaL cut：

  TFile *f = TFile::Open("lambda.root");
  TTree *t = (TTree*)f->Get("LambdaCandidates");
  TCut common = "motherTopoValid==1 && motherTopoFlightValid==1 && "
                "motherTopoPvNdf>0 && motherTopoPvChi2<5*motherTopoPvNdf";
  t->Draw("massGeV>>hBefore(200,1.08,1.18)", common);
  t->Draw("motherTopoMassGeV>>hAfter(200,1.08,1.18)", common);

这里两次使用同一个 common，候选名单相同。对 Lambda 其 PV 有效约束 NDF 为 2。
这演示的是用本次 Lambda 输出进行比较，不是重新生成 PPT 的超氚 toy。

若今后要给八维长度显著性加 >3，可在同一个 common 中追加：

  common = common && "motherTopoSigmaDecayLengthCm>0 && "
                     "motherTopoDecayLengthSignificance>3";

本示例默认没有追加这一项。motherTopoFlightValid 用于确认派生量有效。
若只想加普通几何 cut，则可用：

  TCut geometric = "daughterDcaCm<=1.0 && geometricChi2Ndf<=10.0";
  t->Draw("massGeV>>hMass(200,1.08,1.18)", geometric);

8. 真实数据

先根据标定修改 realT0ScaleToNs：数值已为 ns 时填 1，是 ps 时填 0.001。
realT0SigmaNs 和 realAdditionalTimeSigmaNs 也需核对，默认 0.10 ns 为占位。
然后运行：

  root -l -b -q 'macros/runCeeV0Reco.C("/data/real.list","lambda_real.root",10000,"tof",0,"real")'

有场 real 入口当前仍走通用事件时间选择，默认也需要 CeeT0Hit。
若实际数据采用明确的固定时间约定，可在 macro 中将 realT0Source 改为 fixed，
并正确设置 realFixedT0Ns 及相关误差；不能仅为绕过缺失分支而假定时刻。

9. 调试与运行位置

第七个参数指定 CeeV0Analysis.bin 的路径，默认 ./CeeV0Analysis.bin。
从其他目录运行时，传入程序绝对路径，或先 cd 到 CeeV0Analysis。
第八个参数 true 只打印将执行的命令，不处理数据：

  root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.root","lambda.root",100,"tof",0,"simulation","./CeeV0Analysis.bin",true)'

macro 返回实际子程序退出码。若批任务依赖 shell 退出码，在 ROOT 会话中可调用：

  .L macros/runCeeV0Reco.C
  gSystem->Exit(runCeeV0Reco("/data/reco.root","lambda.root",10000));

运行后可用原项目的 QA macro：

  root -l -b -q 'macros/plotCeeV0.C("lambda.root","lambda_QA.pdf")'

10. 手册对应位置

已核对 CEE_V0Lite_STAR_KFParticle_Principles_Teaching_Manual_v0.6.5.pdf。
以下是 PDF 页码，章节标题数字比页码小 4：
  82 页：编译、最小运行、ROOT/list/目录输入、宏与二进制共用实现。
  83 页：真实数据 T0 单位、偏移和误差。
  84 页：有磁场模式和六维输入要求。
  85 页：legacy/loose/star 预设与逐项覆盖。
  86 页：输入、树名、事件数、文件数等选项。
  87 页：TOF/T0 时间及匹配。
  88 页：track 质量、角色规则与误差回退。
  89 页：候选 cut，含 daughter DCA、几何 chi2、长度及显著性。
  90 页：LambdaCandidates 输出结构。

手册说明了这些接口，但没有逐项展开本示例这样一份可编辑 cut 的完整 .C。
原 runCeeV0Analysis.C 只有较少入口参数。本示例是它的独立扩展示例，
没有修改项目原宏、重建核心、PPT 或手册。

验证范围：C++11 编译与实际 main.cxx 参数解析器验证，检查参数传递、
cut 覆盖、路径引号、模拟 T0 合约、有符号 By 与子进程返回码。
当前环境没有 ROOT/CEE 运行库，未执行真实 ROOT 文件 I/O 或重建闭合。
输入分支和运行说明依据当前 v0.6.5 CeeV0Analysis 源码与上述手册。

ROOT 官方宏运行说明：
https://root.cern.ch/manual/root_macros_and_shared_libraries/

11. v0.6.5 有磁参考位置默认行为

magneticTrackReference="pv-dca"：在正式 Finder 中，每个事件/每条 track
在配对前转换一次，参考位置为测量螺旋线到所选重建 PV 的最近点。
位置、动量、方向、完整 6x6 协方差一起更新，不修改输入 ROOT 文件。
SetVertexSeed 只保存传入值，不能只把位置改成 DCA 点而保留旧动量/协方差。
referenceMaximumAbsPathCm=1000 是从输入参考点起的有符号弧长搜索半范围。
DCA 求解必须收敛；无效输入/失败 track 被跳过并计入 hPvDcaReferenceStatus。
另存 hPvDcaReferencePath、hPvDcaReferenceDca。原 hTrackReferenceSource 仍表示适配层来源。
需要与旧参考约定比较时设 magneticTrackReference="stored"，其他 cut 保持相同。
零场模式忽略这个参考转换开关。母粒子 PV 约束仍由独立开关控制。
本次环境无 ROOT/FairRoot：核心/CLI/语法与宏的 dry-run 已检查，真实 ROOT I/O 待实验环境验证。
