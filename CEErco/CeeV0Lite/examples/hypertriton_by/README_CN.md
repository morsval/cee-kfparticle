# v0.6.5 默认参考更新

`bash build.sh` 后运行 `python3 simulate.py`，默认通过正式 `PrepareTrackReference` 将测量 track 转到重建 PV 的 DCA 参考点。保持原始固定种子、20 cm 处的测量噪声生成和 PID；改变的是传给拟合器的参考状态。`--track-reference stored` 可重现旧参考约定。参考失败候选 raw status=-300；诊断写入 NPZ/JSON。不使用 truth 来确定参考点、seed 或选择。

两条 daughter 电荷为 +2、−1，母粒子 +1；正式 CeeV0Analysis 默认仍是 Λ 通道。本例不是实验 ROOT I/O 验证。旧版本预生成图/数字保留作历史记录；最新同信号回放见 `../../../v0.6.5_review`（从本示例目录向上三级到 CEErco）。

# 用正式 v0.6.3 核心重建带电超氚

这个例子从 2026-09-29 的 `Hypertriton_By_Toy_20260929.zip` 接续而来。
输入、随机种子、tracking/TOF 设定和选择保持一致，PV 步骤改为调用正式 C++ 核心。
它展示 `³ΛH → ³He + π⁻`，电荷为 `+1 → +2 + (-1)`，磁场为 `By=+0.5 T`。

## 代码怎样连接

`core_bridge.cxx` 先调用 `VertexFitter::Fit()` 和完整的母粒子协方差重拟合，
再用 `FitProductionVertex(..., charge, By)` 生成带电八维状态。
示例返回衰变端表示以方便与原测试比较；正式接口默认返回产生端。
两者通过 `TransportToDecayVertex()` 转换，动量和完整协方差一起变化。

`charged_topology.py` 现在只是 C++ 函数的 Python 包装。
`reference_charged_topology.py` 保留原 Python 螺线算法，用于独立数值对照；
正式拟合失败时不会用它补一个结果。

本示例直接使用通用两轨核心，不经过旧中性 V0 的直线 pointing/DCA 选择。
氦三输入采用物理动量 p 及相应协方差，不能使用未转换的刚度 p/2。
TOF 时间差只作 PID 一致性检查；公共起始时刻相消，不要求用户提供 T0。
输入谱形、误差和简化本底沿用原玩具测试，不是 CEE 实际探测器响应。

## 本次结果

- 生成 6,000 个信号和 120,000 个无关组合；97,974 个通过共同重建选择。
- 正式 C++ 与原 Python 对 97,974 个候选逐一对照：同为 97,973 个有效、1 个无效。
- 有效候选中，最大状态差为 7.4×10⁻⁸ 个标准差，协方差尺度化差为 1.94×10⁻⁸。
- 同一批 5,410 个信号：质量 RMS 为 2.692 → 2.650 MeV/c²，与原独立测试一致。
- 额外要求 χ²PV/2<5、S/σS>3 后，保留 4,803 个信号、2,544 个本底。
  这种本底抑制取决于玩具本底设定，不能当作实验性能预测。
- 没有施加超氚质量或寿命约束。只保存 S 本身不会改变不变质量。

原有长尾仍存在：S pull 全样本标准差约 1.19，中央 68% 半宽约 1.00，
有 9 个信号的 |S pull|>5。原分析怀疑存在子轨迹顶点分支问题，但尚未证明根因。
本次没有通过真值删除这些候选，也没有把母粒子的带电升级当作该问题已解决。

## 重现

在当前目录运行，需要 g++、OpenMP、Python、NumPy/SciPy/Matplotlib；无需 ROOT：

```bash
bash build.sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 python3 simulate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 python3 validate.py
OPENBLAS_NUM_THREADS=1 python3 compare_formal_core.py
OPENBLAS_NUM_THREADS=1 python3 plot_results.py
```

本包包含图、汇总、逐候选 CSV、分箱 CSV 和全部源代码。
约 82 MB 的 `output/candidates.npz` 是可重建中间文件，不随包交付；
先运行 `simulate.py`，其余检查便可读取它。

`validation.json` 的 GLS 使用独立 SciPy 八参数拟合，样本为 32 个信号。
其中明确标成 Python reference 的端点检查用于核对参考公式；
正式 C++ 端点、S 导数、零方差 PV、正负电荷与磁场等检查位于
`../../tests/test_charged_production.cxx` 与 `../../tests/verify_charged_reference.py`。
实际 CEE/FairRoot 链接、ROOT 文件读写、材料效应、真实效率和误差闭合仍未完成。

## 输入来源

本次复用原独立测试的 `data/spectrum_fit.json`、谱形 CSV 与原始设置。
详细物理依据已附在 `ORIGINAL_TOY_NOTES_20260929.md`。
谱形来自原测试恢复的 STAR 3 GeV、0–10% 曲线；分辨率系数是玩具参数，
不是经过验证的 CEE 或 STAR 氦三 embedding 响应。
