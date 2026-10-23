# CEE 超氚：L/sigmaL 选择与 7/8 维一致性检查

2026-10-14，基于 CEErco v0.6.4。核心拟合源码保持不变；本目录增加受控 toy 比较和独立参考实现。

## 1. 这次采用什么量

令 p=|p|，S 为飞行参数，L=pS 为从 PV 到 SV 的有符号三维螺旋线弧长，单位 cm。sigmaL 是 L 的标准差，L/sigmaL 是衰变长度显著性。正的较大数值表示拟合出的衰变位移相对于其误差较明确；它不是信号概率，也不是不变质量谱的 S/sqrt(S+B)。

使用完整 8×8 协方差的一阶误差传播：

    sigmaL² = p² Var(S) + S² Var(p) + 2 S p Cov(S,p)
    Var(p)  = uᵀ Cpp u,  u=p_vector/|p_vector|
    Cov(S,p)= CS,p_vector u

等价梯度为 (0,0,0,S ux,S uy,S uz,0,p)，sigmaL²=gᵀ C8 g。不能只用 p*sigmaS，也不能忽略相关项。L/sigmaL 与 S/sigmaS 一般不相同。

v0.6.4 的 GetFlightObservables 已经按这个完整梯度计算，本次直接调用正式 C++ 接口，并以标量展开式交叉核对，最大相对方差差 4.16e-14。因此没有必要替换核心误差算法。

STAR KFParticleBaseSIMD::GetDecayLength 同样包含这三个贡献，不能据此宣称本实现比 STAR 更准确。用户提供的 StKFParticleAnalysisMaker 使用的是 PV 重拟合前的 GetDistanceToVertexLine，那个量与本次 PV 重拟合后的有符号弧长不同。

## 2. 受控质量谱比较

两侧完全共用输入、共同 SV 拟合、带电母粒子 PV 重拟合、质量计算及 chi2PV/2<5。右侧只增加 L/sigmaL>3。没有加入质量约束或寿命约束，阈值 3 是本次演示的工作点，未声称最优。

| 选择 | 信号总数 | 本底总数 | 质量窗内信号 | 质量窗内本底 |
|---|---:|---:|---:|---:|
| PV 重拟合有效，chi2PV/2<5 | 5381 | 90373 | 5200 | 10172 |
| 同上，再加 L/sigmaL>3 | 4805 | 2551 | 4647 | 281 |

质量窗为 |M−2.99131|<0.006 GeV/c²。相对左侧，信号保留 89.2957%，本底剩余 2.82275%。总数包含绘图质量范围外的候选。两图使用相同的 2.95–3.06 GeV/c² 范围、0.5 MeV/c² 分箱和 0–1000 纵轴。

相同 PV 条件下，旧 S/sigmaS>3 选择为 4803 个信号、2544 个本底。改为 L/sigmaL 后有 9 个候选改变（2 个信号、7 个本底），不能只把旧图的轴标/标签改名。

toy 沿用随机种子 20260929、6000 个信号、120000 个组合本底、By=+0.5 T，母粒子电荷 +1、子粒子电荷 +2 和 −1，PV 各坐标误差 0.05 cm。输入接受后 109706 个组合；共同 SV 选择要求正式拟合成功、几何 chi2<10、daughter DCA<1 cm，余 97974 个。该 toy 的组合本底含 80% prompt 设定，抑制率不等于真实 CEE 性能预测。输入细节见 results/toy_config.json。

## 3. 正确的 7/8 维比较

profiled_7d.py 以同一七维母粒子测量和三维 PV 测量为输入，包含完整母粒子/PV 交叉协方差。S 是自由飞行参数；每次线性化后解析消去它的变化量，而不是令 S=0。

若三维几何残差对测量的雅可比为 H，对 S 的导数（差一个符号）为 v，W=H C Hᵀ，则消去 S 后的权重为：

    Q = W⁻¹ − W⁻¹ v (vᵀ W⁻¹ v)⁻¹ vᵀ W⁻¹

用 Q 更新并返回 q7 及其边缘协方差。正式 C++ 通过等价的二维横向约束拟合并输出完整八维状态。比较在相同衰变端、相同先验和相同 PV cut 下进行。显式存储第八维没有加入新的测量。

97974 个候选：97973 个双方均有效，1 个双方均失败，无单侧成功，PV cut 判定全部一致。共同有效候选的最大差异：

| 量 | 最大差异 |
|---|---:|
| 七维状态，按各分量 sigma 归一 | 2.60e-12 |
| 协方差，按 sigma_i*sigma_j 归一 | 1.02e-11 |
| 质量 | 4.89e-15 GeV/c² |
| chi2PV | 2.81e-12 |

verify_explicit_8d.py 另用显式八参数（产生端位置、动量、能量、S）的独立非线性最小二乘验证，完整白化十维测量协方差。112 个 toy 抽样（含 cut 边界附近）与 16 个含电荷/磁场符号和非零交叉协方差的测试案例全部通过。最大状态差 3.92e-8 sigma，尺度化协方差差 5.58e-10，质量差 6.01e-11 GeV/c²，chi2 差 5.33e-14。详见 JSON；这是数值交叉验证，不是对所有输入的精确数学恒等性保证。

S 的不确定性和相关性在联合拟合或正确消元中已经处理。chi2PV 衡量母粒子轨迹与 PV 的相容性，L/sigmaL 衡量沿飞行方向的位移显著性；二者提供不同选择信息，但不要求统计独立。不能把拟合后的 sigmaS 再当作独立测量误差额外加进同一次 PV 残差，否则可能重复计入不确定性。

## 4. 复现

需要 g++、make、Python 3、NumPy、SciPy、Matplotlib。在解压后的 CEErco 目录执行：

```bash
cd CeeV0Lite/examples/hypertriton_by
bash build.sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 python3 simulate.py --out output_LsigmaL_input
python3 flight_study/run_flight_study.py --cee-dir ../../.. --input output_LsigmaL_input/candidates.npz --out flight_study/results
python3 flight_study/verify_explicit_8d.py --cee-dir ../../.. --input output_LsigmaL_input/candidates.npz --out flight_study/results
```

simulate.py 保留原历史例子的 selected（S/sigmaS）标记。本次 run_flight_study.py 不读取那个标记，而是重新构建 left_PV_chi2 和 right_LsigmaL，因此不会混用两个工作点。输出目录单独命名，避免覆盖旧质量谱。本次生成阶段复用了先前已验证的同输入 v0.6.4 共同顶点结果缓存，并重算 PV 拟合；前 40 个输入另与正式核心逐位核对。上述完整复现命令不依赖缓存。

results/candidate_comparison.csv 包含共同 SV 选择后的逐候选 ID、真值类别、有效性、左右选择、旧 S 选择、质量、L、sigmaL、S、sigmaS 与两种实现的 chi2。summary.json、histograms.json 是 PPT 新页的直接数据来源；PNG/PDF 是独立质量谱导出。

## 5. 实际分析中使用哪个分支

已有 CeeV0Analysis 输出提供 motherTopoDecayLengthCm、motherTopoSigmaDecayLengthCm、motherTopoDecayLengthSignificance。新的量是 PV 重拟合副本的值；旧 decayLengthSignificance 与 --min-ldl 仍沿用旧定义，不能混用。

对应本次工作点，可以在已完成共同候选选择的输出上使用：

```cpp
motherTopoValid && motherTopoFlightValid &&
motherTopoPvNdf == 2 && motherTopoPvChi2 / motherTopoPvNdf < 5 &&
motherTopoSigmaDecayLengthCm > 0 &&
motherTopoDecayLengthSignificance > 3
```

画对应的母粒子质量应取 motherTopoMassGeV。这里的母粒子 PV 增量 chi2 不是 daughter chi2primary，也不是包含先前 SV chi2 的总 motherTopoChi2Ndf。

本补充没有改变全局默认选择；它提供上述明确的受控工作点与验证程序。
