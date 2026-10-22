# CEE v0.6.1：磁场方向修正（2026-09-20）

CEE 束流沿 z 轴，磁场沿 y 轴。本补丁将原先错误的均匀 Bz 实现
改为 B=(0,By,0)。By 的单位为 T；正值表示 +y，负值表示 -y。
磁场数值及正负号应使用实际运行条件，本补丁不指定实验场强。

## 必须同步修改的接口

- C++：`config.magneticFieldBzTesla` 改为 `config.magneticFieldByTesla`。
- 命令行：`--bz-tesla VALUE` 改为 `--by-tesla VALUE`。
- `PropagateState`、`TransportToClosestApproach`、`VertexFitter::Fit` 的
  位置参数数量不变，但磁场数值现在明确表示 By。
- 旧命令 `--bz-tesla 0` 仅为零场兼容保留，会打印提示；非零值拒绝运行。
  不把非零 Bz 悄悄当作 By，也不保留含义错误的 C++ 字段别名。
- 替换后重新编译两个模块。v0.6.1 版本号及 outputFormat=4 保持不变。

```bash
# 零场模拟，原有物理配置保持不变
./CeeV0Analysis.bin --input /data/cee_sim.root --output lambda_legacy.root \
  --mode simulation --momentum tof --by-tesla 0 \
  --cut-preset legacy --events 10000 --debug 1

# 有场示例；0.5 仅演示参数写法，应换成实际带符号 By
./CeeV0Analysis.bin --input /data/cee_tracking.root --output lambda_by.root \
  --mode simulation --momentum tracking --by-tesla 0.5 --cut-preset legacy
```

## 物理和数据约定

- 正电荷沿 +z 运动、磁场为 +y 时，洛伦兹力指向 -x。
- 动量在 x-z 平面旋转，py 和 |p| 保持不变；沿 y 匀速运动。
- 最近点算法中的“纵向”现在指沿磁场的 y，不能理解成束流方向 z。
- 所有输入输出仍使用实验室 (x,y,z,px,py,pz)，协方差仍按同一顺序；
  不交换 ROOT 分支、TPC 参数或数据类的坐标含义。
- 数值雅各比在修正后的输运/最近点函数上重新求导，包含交叉协方差。
- 有场仍要求可靠的跟踪动量、非零电荷和协方差。TOF 动量路径依旧只
  用于零场，不因此获得带曲线路程的 TOF 重建能力。
- 输出配置摘要记录 ByTesla、magneticFieldAxis=y、beamAxis=z、
  fieldConvention=CEE-By-20260920，旧的非零 Bz 结果需要重新运行。
- 本补丁仅支持均匀 By，没有新增磁场图或材料模型；STAR 参考代码不改。

## 检查

```bash
make -C CeeV0Lite test
make -C CeeV0Analysis check
```

新增测试以独立 RK4 洛伦兹力积分校验正负电荷、正负磁场、正反向路径、
沿场直线、全六维协方差和多圈最近点；现有顶点及 Lambda 重建测试保留。
CLI 测试使用真实参数解析器，只替换 ROOT 事件读取部分。
这些检查不能替代实际 CEE/ROOT 环境的链接和真实数据物理验证。
