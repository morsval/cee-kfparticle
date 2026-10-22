# v0.6.4：保守的输运计算优化

基线：2026-09-30 正式发布的 `CEErco_v0.6.3.zip`。
本次要求是优先保留能确认收益、且对现有拟合影响很小的修改。

## 保留的修改

- `PropagateState`：固定路径的直线/均匀 By 输运使用解析 Jacobian，
  省去构造数值导数的 12 次参数输运。状态中心值仍使用原公式。
  极小动量、直线/螺旋线分支阈值附近或非有限导数回退到原差分。
- `PropagateCovariance`：将四重累加改为两次矩阵乘法 `A=JC; C'=AJ^T`，
  不丢弃位置与动量的交叉项。噪声加入与对称化顺序沿用原版。
- POCA 的三个临时数组预留容量，避免反复扩容，不改变候选顺序或根的选择。
- 输出运行元数据标记 `ceeV0LiteVersion=0.6.4` 与 `analyticFixedPathJacobian`。

## 保留的拟合流程

`TransportToClosestApproach` 仍执行中心求解及 12 次扰动求解，包含路径对输入的导数。
共同 SV 的内部迭代、零场 TOF 外层迭代、Kalman 残差与更新、最终母粒子
7×12 Jacobian 的 24 次重拟合、PV 约束和候选 cuts 均沿用 v0.6.3。
几何求根的阈值、近平行判断、迭代上限和失败政策没有修改。

没有采用 STAR 成对 POCA 替换，也没有把最近点中点直接当作一个新的拟合定义。
试验过的隐式 POCA 导数在全样本对照中改变了 6 个拟合状态与 2 个最终选择，
少数质量误差变化较大，因此撤回。交付代码不含该实验分支。

## 配置与重新编译

新增 C++ 配置字段默认开启：

```cpp
CeeV0Lite::TransportConfig cfg;
cfg.useAnalyticFixedPathJacobian = true;
// false restores the old fixed-path finite-difference derivative.
// It does not change the POCA derivative; the matrix optimization stays on.
```

没有新增命令行开关。该字段也可从 `FinderConfig::transportConfig` 设置。
ROOT 输出树格式仍为 6。因为配置结构新增字段，必须清理并重新编译相关模块：

```bash
cd CEErco/CeeV0Lite
make clean
make -j4 test
make reference-test
make all
cd ../CeeV0Analysis
make clean
make check
# In the configured CEE / ROOT / FairRoot environment:
make -j4
```

当前环境验证的是独立核心、分析层检查及合成输入，未完成实际 CEE/FairRoot 链接和 ROOT I/O。

## 验证证据

见 `../v0.6.4_review/comparison.json`、测试日志及新版手册的“输运导数优化”章节。
完整 v0.6.3→v0.6.4 文本补丁及逐文件 SHA-256 也在该目录。
八项核心测试包括新增 1,200 组导数对照；独立中性/带电 GLS 与洛伦兹力积分参考测试保留。

本包 `examples/hypertriton_by/output` 中的原图和表是 v0.6.3 的历史结果，未冒充本次重跑。
v0.6.4 的同输入逐候选验证单独存放在 `v0.6.4_review`。
保存 S、施加 PV 约束及重建后的 S/σS 选择仍为不同步骤。
S=L/|p|，带电母粒子的 L 为有符号三维弧长，S 不是“显著性”的缩写。
