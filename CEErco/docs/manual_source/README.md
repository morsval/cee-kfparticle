# v0.7.0 中文手册源稿

`rewrite_cee_manual_v070.py` 生成 129 页 PDF。新增用户分析层为第 5-24 页，其中第 14-23 页是 89 个可选标量表，第 24 页补充完整状态和诊断。原物理原理、TOF 迭代、Kalman 更新、带磁场输运与 PV 状态章节保留。

```bash
CEE_MANUAL_OUTPUT=/absolute/path/manual.pdf python3 rewrite_cee_manual_v070.py
```

目录应保持在 `CEErco/docs/manual_source`，以便计算源码指纹并读取相邻的 CandidateVariables JSON。依赖 Python reportlab、matplotlib、fontTools；字体子集和 OFL 随目录提供。`cee_rewrite_page_map.json` 为构建产物，记录实际页码。
