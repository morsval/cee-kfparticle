# 七维与八维：相同 L/sigmaL 选择

2026-10-15。本目录补充 v0.6.4 的验证，正式 C++ 重建核心未修改。

## 结论与比较口径

不把 S 存进粒子的状态向量，并不等于不能计算 sigmaL。只要保留输入母粒子与 PV 的误差信息、相关性和拟合响应，就能将 L 作为派生标量，单独计算其误差。只保留拟合后的七维边缘协方差、丢弃 PV 及相关信息时，一般无法完整恢复同一个 sigmaL。

之前的 flight_study 比较的是相同 PV 条件下加不加 L/sigmaL cut，反映长度显著性选择的效果。本目录改为左右两侧均有：

1. 同一输入、共同 SV 选择和带电母粒子螺旋线 PV 重拟合。
2. chi2PV/2 < 5。
3. L/sigmaL > 3。

左侧为独立 NumPy 七维状态实现，S 作为局部飞行参数消去，不构造八维状态或 8×8 协方差，不读取右侧的 L 或 sigmaL。右侧为正式 C++ 八维结果及其完整协方差。

## 七维实现怎样计算 sigmaL

输入 y=(母粒子 q7, PV xyz)，联合协方差为 C0。PV 是外部测量，粒子状态仍为 q7。消去飞行参数是求解方式，与最终是否保存第八维不同。

螺旋线几何约束写成 g=rD-rP-A(s)pD=0。在线性化位置，令 H=dg/dx、v=d(Ap)/ds、W=H C0 H^T。消去 s 的局部变化后：

    Q  = W^-1 - W^-1 v (v^T W^-1 v)^-1 v^T W^-1
    Jx = I - C0 H^T Q H
    Js = (v^T W^-1 H) / (v^T W^-1 v)

Jx 表示拟合状态对原始输入的响应，Js 表示临时飞行参数对原始输入的响应。取 u=p_vector/|p_vector|：

    Jp = u^T Jx[动量三行,:]
    JL = p Js + s Jp
    L = p s
    sigmaL² = JL C0 JL^T

这里直接计算 L 的梯度，不把 sigmaS 或 S 与七维状态的交叉协方差拼成 C8。所有计算都是收敛点的一阶协方差传播，与八维实现的线性化约定一致。七维输出另有 L、sigmaL 两个派生量，正如质量与质量误差也不必各占一个状态维度。

八维实现把同一响应整理成 q8=(q7,S) 和 C8，再用 L=pS 的梯度传播。只要输入信息和近似相同，两种表达应一致。

## 实际结果

同一 toy：By=+0.5 T，种子 20260929，生成 6000 信号和 120000 组合本底。共同 SV 选择后 97974 个候选，97973 个双方均有效，1 个双方均失败。

| 条件/结果 | 七维＋派生 L/sigmaL | 八维＋完整 C8 |
|---|---:|---:|
| 母粒子 PV 重拟合 | 有 | 有 |
| chi2PV/2 | <5 | <5 |
| L/sigmaL | >3 | >3 |
| 选择后信号总数 | 4805 | 4805 |
| 选择后本底总数 | 2551 | 2551 |
| 质量窗内信号 | 4647 | 4647 |
| 质量窗内本底 | 281 | 281 |

质量窗为 |M-2.99131|<0.006 GeV/c²，半宽 6 MeV/c²，全宽 12 MeV/c²。窗内计数逐候选计算，不把窗口边界近似到 bin 中心。

两侧有效性、PV cut 和最终 cut 判定全部一致，独立填充的质量谱逐 bin 一致。最大 sigmaL 相对差 2.81e-13，L/sigmaL 绝对差 1.84e-11，L 差 1.81e-12 cm，质量差 4.89e-15 GeV/c²。另有 16 个包含电荷/磁场符号和母粒子/PV 非零交叉协方差的案例，最大 sigmaL 相对差 1.42e-14，全部通过。

因此，本次比较不支持“增加第八维本身提高重建能力”的结论。第八维的实用价值是集中保存飞行参数及其相关协方差，便于后续输运、组合和误差传播，避免下游丢失信息。七维实现若保留等价信息，也能得到相同的长度显著性选择。

这里比较的是正确处理误差的等价实现。把七维一侧的相关项删掉、将 S 固定为零，或改用不同几何量，会同时改变误差模型或约束，不能把这些差异归因于状态维数。

## 复现

依赖与原 toy 相同：g++、make、Python、NumPy、SciPy、Matplotlib。在 CEErco/CeeV0Lite/examples/hypertriton_by 下执行：

```bash
bash build.sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 python3 simulate.py --out output_LsigmaL_input
OPENBLAS_NUM_THREADS=1 python3 flight_study/equivalent_selection/run_equivalent_selection.py --input output_LsigmaL_input/candidates.npz --out flight_study/equivalent_selection/results
OPENBLAS_NUM_THREADS=1 python3 flight_study/equivalent_selection/check_correlated_fixtures.py --cee-dir ../../.. --out flight_study/equivalent_selection/results
python3 flight_study/equivalent_selection/plot_equivalent.py --results flight_study/equivalent_selection/results
```

若已有同一 toy 的 candidates.npz，可直接复用，省去生成和共同 SV 拟合。本目录分析不使用旧的 selected 标记。equivalent_candidates.csv 记录逐候选的两侧有效性、选择、质量、L、sigmaL、L/sigmaL 与 chi2PV。

配套 PPT 为 CEE_KFP_Intro_Revised_v10.pptx，第 15 页展示这个等价选择比较。图形沿用灰色填充本底、红色阶梯信号、蓝色总计数与 sqrt(N) 误差棒，两边均标出质量统计窗口。toy 本底含 80% prompt 组合假设，不能直接作为真实 CEE 的本底抑制率预测。
