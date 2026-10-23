# CEE KFParticle v0.7.0：用户分析层

原来已有 `CeeV0Analysis::Analyzer`，它在一个类里完成读 ROOT、选 track、TOF 匹配、Finder 重建和固定树输出。它能跑分析，但不适合频繁替换 PID、输入格式或母粒子树结构。原 `AnaMaker` 也是旧分析实现，并非这里所需的格式无关 KFP 用户接口。

v0.7.0 新增 `CeeKFParticleAnalysisMaker.h/.cxx`，借鉴 STAR 的 Maker → Interface → 重建结果这一层次；不引入 StMaker、singleton、StChain、TMVA 或 STAR 数据格式依赖。物理拟合仍由现有 CeeV0Lite 完成。

## 组件的职责

|组件|用户在哪里修改|职责|
|---|---|---|
|`CeeKFParticleRootSource`|原始输入|读 ROOT/list、PV、track，零场匹配 TOF|
|`CeeKFParticleBenoSource`|BenoDst 输入|读 ROOT/list，复用保存的匹配、对齐位置和已扣 T0 的时间|
|`ConfigureInputModes` / `OpenCeeInput`|macro|独立选择格式与动量来源，打开 Source|
|`CeeKFParticleAnalysisMaker`|派生类或 callback|事件和 track 选择、PID、组对前的选中列表、重建后的母粒子选择|
|`CeeKFParticleInterface`|`maker.GetInterface().SetConfig(...)`|连接 Finder，接收选中 track，提供原始候选及参考点诊断|
|`CeeV0Lite::Finder`|`FinderConfig` 的基础 cut|子粒子组合、共同顶点、Kalman 更新、原 7D 母粒子和附加 8D PV 约束副本|
|`CeeKFParticleRootWriter` / 自定义 Sink|`SetVariables`、`AddColumn` 或自己写树|选择分支、保存自定义标量/协方差，或完全接管 ROOT 输出|

PV = Primary Vertex；SV = Secondary Vertex；POCA = Point of Closest Approach；DCA = Distance of Closest Approach。带磁场 track 的默认参考点仍为 track–PV DCA 点；这里只改变参考位置、动量及完整协方差的表示，不对子 track 加 PV 拟合约束。

## 从 macro 直接控制

在配置好 ROOT、FairRoot 和 CEE 的环境中：

```bash
cd CEErco/CeeV0Analysis
make shared
root -l -b -q 'macros/runCeeKFParticle.C("input.root","mothers.root","raw","tracking",0.5)'
root -l -b -q 'macros/runCeeKFParticle.C("files.list","tof.root","benodst","tof",0.,1000,true)'
```

参数依次是 `input, output, inputFormat, momentumMode, ByTesla, maximumEvents, realData, overwrite`。format 为 raw/benodst，mode 为 tof/tracking/auto。旧数值第三参数 By 的重载仍可用。Beno 实验数据已减 T0，但 `config.t0.sigmaNs` 仍是共同误差，TOF hit 误差中不能重复包含它。`maximumEvents=-1` 为全部；列表一行一个 `.root` 路径，允许空行和以 `#` 开头的注释，路径相对于运行目录；支持 ROOT 可访问的远程文件路径。磁场采用带符号的均匀 **By**，束流为 z。宏直接构造并调用类，不再通过 shell 启动分析可执行文件。ROOT 6 的 Cling 加载共享库；Makefile 使用本机 root-config 指定的 C++ 标准。既有 CLI 与 `runCeeV0Reco.C` 仍保留。

```cpp
AnalysisConfig config;
// 在这里确定输入模式、By、TOF/T0 标定及 Finder 基础 cut。
CeeKFParticleRootSource source(config);
CeeKFParticleAnalysisMaker maker(config.finder);
maker.SetPIDSelector(myPidFunction);
maker.SetMotherSelector(myMotherCut);
CeeKFParticleRootWriter writer("mothers.root");
writer.SetVariables({"massGeV", "DcaProton", "DcaPion"});
std::string error;
if (!maker.Run(source, writer, error)) throw std::runtime_error(error);
```

完整可运行宏在 `CeeV0Analysis/macros/runCeeKFParticle.C`；本节短片段只示意接口，不能省略输入和模式配置。

## PID 和 track 选择由用户负责

`AnalysisTrack` 同时保存格式无关的 `TrackInput` 与 `nHits`、`chi2Ndf`、`fitConverged`、原始 `dEdx`，以及可选的 `nSigmaProton`、`nSigmaPion`、`tofMass2GeV2`。**当前 CEE reader 不虚构 PID 标定**：后三者默认 NaN。有磁模式也不会用直线投影随意关联 TOF hit。需要用户在 reader 或自己的 PID 函数中接入实际标定、匹配结果。

```cpp
maker.SetPIDSelector([](const AnalysisTrack& t,const AnalysisEvent&) {
  unsigned mask=0;
  if (std::isfinite(t.nSigmaProton) && std::fabs(t.nSigmaProton)<3.)
    mask |= CeeV0Lite::kProtonHypothesis;
  if (std::isfinite(t.nSigmaPion) && std::fabs(t.nSigmaPion)<3.)
    mask |= CeeV0Lite::kPionHypothesis;
  return mask;
});
```

这段示例要求输入已经提供标定后的 nSigma；直接接到尚未填 nSigma 的 reader 会拒绝所有 track。返回 0 拒绝 track，返回一个 bit 只接受对应假设，两个 bit 同时保留歧义。实测电荷由 Finder 检查；没有通过真值决定 PID。默认 IdentifyTrack 仅保留输入 hypothesisMask，**不是已经做过 PID 选择**。

可以通过 `MakerTrackCuts` 设置 hits、chi2/NDF、fit convergence，也可重写 `AcceptTrack` / `AcceptEvent`。派生类覆盖这些方法时，若仍要保留基类配置，应显式调用基类。示例宏覆盖 IdentifyTrack 和 AcceptMother；其余走基类。

## 只重建你指定的 track

三个入口对应不同阶段：

1. `source.SetTrackIndexSelector(...)`：返回原 ROOT TClonesArray 的索引，先筛选再适配和 TOF 匹配；适合已有 good-track ID 列表。callback 可读取原始 CeeTpcTrack。原数组索引保存在 `TrackInput.id` 中，不重新编号。
2. `maker.SetTrackListSelector(...)`：返回当前 `AnalysisEvent.tracks` 的索引，只对这些对象执行 track/PID 选择，再送入 Finder。注意它不是原 ROOT 数组索引；需要原始 ID 时查看 `tracks[i].input.id`。
3. 未来 picoDst reader 可以直接只填需要的 `event.tracks`，随后调用 `maker.Make(event,error)`。

空列表明确表示不选任何 track；重复或越界索引、重复或负的 track ID 是错误。每个事件都会清空上一事件的选中 track、候选和诊断。读取原始全事件并进行首次筛选通常仍需一次线性遍历；减少的是后续 pair 数量及昂贵拟合。没有声称本版已有实测加速倍数。

零场模式的 track–TOF 唯一匹配在 reader 中进行。若某些 track 不应参与匹配，应使用 source 的索引选择器提前排除；Maker 的后续 PID 不会重新分配已匹配的 hit。`maximumTracksPerEvent` 超限在新 reader 中报错，不再按输入顺序悄悄截断。原 CLI 的历史行为保留。

## 基础 cut 与母粒子 cut 的顺序

Finder 的 `lambdaCuts` / `k0Cuts` 是重建内的基础筛选，作用于原始候选，提前拒绝的对象不会进入 Maker。之后才构造附加母粒子 PV 约束副本，再由 Maker 的 `AcceptMother` / `SetMotherSelector` 检查完整结果。用户需要某候选时，应先把 Finder 的对应基础 cut 放宽，不能靠 Maker 恢复已被前一级拒绝的候选。

```cpp
maker.SetMotherSelector([](const Candidate& c,const AnalysisEvent&) {
  if (!(c.massGeV>1.05 && c.massGeV<1.25)) return false;
  if (!(DcaProton(c)>0.1 && DcaPion(c)>0.1)) return false;
  if (!c.motherTopo.valid) return false;
  if (!(c.motherTopo.pvConstraintNdf>0 &&
        c.motherTopo.pvConstraintChi2/c.motherTopo.pvConstraintNdf<5.)) return false;
  return true;
});
```

这些阈值只是展示写法，不是优化好的 CEE 物理工作点。示例主宏不默认加 L/sigmaL 或 S/sigmaS cut。显式 `std::isfinite(x)` 或 `!(x > threshold)` 可防止 NaN 漏过筛选；不要只写 `if (x < threshold) reject` 来检查可能无效的误差量。

`maker.GetCandidates()` 返回通过用户母粒子选择的候选；`maker.GetInterface().GetCandidates()` 返回通过 Finder 基础 cut、尚未经过用户母粒子选择的候选。两者都是当前事件的只读引用，下次 Make 会更新；需跨事件保留时复制。接口还提供 `GetReferenceDiagnostics()`，可统计参考点失败原因。

## DcaProton / DcaPion 原本在哪里

原量早已存在于 `Candidate::daughterDcaToPvCm[2]`，旧 ROOT 树名为 `daughterDcaToPv0Cm/1Cm`。Lambda 中 daughter 0 是质子、daughter 1 是 pion；新增 `DcaProton(c)` / `DcaPion(c)` 依据绝对 PDG 查角色，反 Lambda 同样适用，避免写死数组下标。K0S 的两个 pion 角色有歧义，`DcaPion` 返回 NaN；请使用显式的 0/1 分支。

具体实现：`MeasureDaughterAtPv(solution.inputStates[i],...)` 把该质量假设下、共同顶点更新前的子 track 状态输运到 PV 的 POCA，再算 `|r_POCA - r_PV|`。有磁场使用螺旋线，零场使用直线；单位 cm，无符号、三维，**不是**两条 daughter 的相互 DCA，也不是 PV chi2。零场动量仍依赖已收敛的 TOF 解。

`Chi2PrimaryProton/Pion` 对应 `daughterChi2FromPv`，用横向残差及 track/PV 协方差计算原始 chi2；`daughterPvSignificance` 才是它的平方根。与 STAR 的命名习惯相似，不保证其输入误差、几何定义或数值阈值可直接照搬。

## 取得结果、选择树分支或自己写树

```cpp
writer.SetVariables({"massGeV","DcaProton","DcaPion","motherTopoValid",
                     "topoMassGeV","topoLCm","topoSigmaLCm"});
writer.AddColumn("massMeV",[](const Candidate& c,const AnalysisEvent&){
  return 1000.*c.massGeV;
},"MeV/c^2");
writer.SaveCovariances(true);
```

`SetVariables` 替换整个可选列表（也会清掉先前 AddColumn 的列）；因此先 SetVariables、再 AddColumn。未知/重复名字立即报错；Begin 后配置锁定。固定分支 `eventId` 为 64 位整数、`candidateIndex` 为整数；89 个可选内建标量和自定义列为 double。元数据记录版本、列名/单位/阶段/意义和用户提供的配置文本。

完整变量表见 `CandidateVariables_v0.7.0.md/.csv`。`SaveCovariances(true)` 增加 row-major 的母粒子 7×7、PV 副本 8×8、两条 daughter 6×6 矩阵；无效母粒子矩阵写 NaN。通过 Candidate 本身还可读取 SV 的 3×3 协方差、完整 daughter 状态和 Jacobian 故障诊断。

对于 event 数组、自定义对象树、分渠道树等需求，实现 `CeeKFParticleCandidateSink` 的 Begin/Write/Finish；`examples/UserTreeSink.h` 是用户完全拥有分支和 Fill 逻辑的例子。Sink 在每个接受的事件上调用一次 Write，即使候选数为 0。默认平面母粒子树只为候选写行，空事件不会产生虚假的母粒子。

ROOT 默认 CREATE 保护已有输出，只有明确 `SetOverwrite(true)` 或宏最后一个参数 true 才 RECREATE。读写失败会停止并返回错误，Run 仍调用 Finish 关闭输出；发生错误时留下的文件可能是部分结果，不能当作完整处理成功。Run 返回 false 时检查 error。类不可复制的 Source/Writer 独立拥有 I/O 对象，不要让多个线程共享同一 Maker/Source/Writer 实例。

## 替换为 picoDst 类输入

`examples/PicoDstAdapterExample.h` 给出可编译的 callback 适配器。映射 `AnalysisEvent` 的 PV/误差/T0，以及每条选中 track 的参考位置、动量、实测电荷、完整 6×6 协方差、PID 和原始 ID 即可。磁场输入的 r/p/C 必须位于同一参考点，动量为物理 GeV/c 而非刚度；不要另行对子 track 加 PV 约束。未来输入的读取方式只出现在新 EventSource，Maker、Finder 和 Writer 无须知道 picoDst 的类名。该例子没有虚构一个目前不存在的实验 picoDst schema。

## 迁移范围与验证

新接口支持 Finder 已有的 Lambda、反 Lambda、K0S 开关；默认宏只开 Lambda。超氚仍为已有独立 toy/拟合示例，本次没有把 helium PID 和任意衰变通道加入通用 Finder。旧 Analyzer/CLI 固定 Lambda tree、QA 及 truth 匹配仍保留；新 Reader 默认不读 MCTrack，新 Writer 不自动复制旧 QA 直方图。需要时在自定义 Sink 中增加 QA。

本版以先前交付的 v0.6.5 为基线，并首次把此前本地交付的轻量模块接入此 GitHub 分支。没有调整 Finder、transport、vertex fitter 或 production fitter 的计算公式。测试包含：9 组核心回归、旧分析层检查、Maker 对直接 Finder 的零场/有磁数值一致性、PID/子集/重复 ID/事件清空/母粒子 cut、真实 ROOT 的写出再读回，以及真实 ROOT + 合成 CEE 形状字典的多文件读入测试。真实实验的 FairRoot/CEE 字典链接及真实数据文件仍需在用户环境验证；合成字典不等于实验 EDM 验证。
