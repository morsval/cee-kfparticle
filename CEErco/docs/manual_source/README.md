# v0.7.0 用户与算法手册

当前入口为 `build_user_book_v070.py`，按书籍工作流组织：软件结构、编译、输入、用户分析、衰变配置、macro/ROOT、数据与 cut 表、共同基础、零场 TOF、有磁 tracking、PV/S、验证、附录。

包含自动页码目录、章节编号和 PDF 书签；保留经核对的原推导与示意图，连续排版，不再按版本追加章节。`book_user_sections.py` 是新的使用部分。旧 `rewrite_cee_manual_v070.py` 是兼容图文素材来源。

```bash
CEE_MANUAL_OUTPUT=/absolute/path/manual.pdf python3 build_user_book_v070.py
```

依赖 reportlab、matplotlib、fontTools；附字体子集和 OFL。新增中文文字时需要重新生成覆盖这些字符的字体子集。`book_page_map.json` 记录实际章节页码。
