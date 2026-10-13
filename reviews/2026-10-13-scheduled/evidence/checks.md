# 检查记录

- `python check/physics_regression.py`：PASS 3380；汇总 FAIL 0 / WARN 0 / 3381 项（含 1 INFO）。本轮包括新附录 AK 的回归。
- `python check/figure_layout.py`：12 PASS、0 WARN、0 FAIL；插图文字最小字号 6.5 pt，未见重叠/裁切。
- `python check/pdf_outline.py`：0 FAIL；目录 482 条匹配 482 个书签、顺序一致、目录链接一致；WARN：19/482 书签页差一页；附录 F 的 35 条链接中 29 条与 `.aux` 页记录不一致。
- `python check/full_check.py`：37 项检查 0 FAIL / 0 WARN；编译 PDF 487 页，5 个轻微 overfull（最大 4.8 pt，可忽略），真实缺字 0；退出码 0。
- Git 基线：`3926bfccf76627572893a4367ed72bb0ab0aae66`；评审前工作区干净。

物理回归是自动化核对，不构成 37 个附录的全量物理推导审查。
