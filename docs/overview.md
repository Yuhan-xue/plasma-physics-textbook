# 项目一览（一页看完这本书和这个仓库）

![项目一览](overview.png)

上图可直接点开看大图；下面是把图里的数字用文字再说一遍，方便检索与复制。

## 这本书现在是什么

| | |
|---|---|
| 规模 | **657 页**：12 章正文 176 页 + 52 个专题附录 449 页，另有封面目录 18 页、前言 5 页、结语与索引 9 页 |
| 配图 | **95 张图**，全部是实拍照片或公开示意图，逐张登记作者、许可与来源页（[`figures/CREDITS.md`](../figures/CREDITS.md)） |
| 交付 | 三个 PDF：[完整本 657 页](../pdf/plasma_physics_textbook_v2.pdf)、[卷一 208 页](../pdf/plasma_physics_textbook_v2_vol1_chapters.pdf)、[卷二 449 页](../pdf/plasma_physics_textbook_v2_vol2_appendices.pdf) |
| 版本 | v1.27.3；每版都发到 [Gitea releases](http://zsyq.hxlab.tech:3000/yuhanxue/plasma-physics-textbook/releases) |
| 动手 | 正文不印程序，需要动手的地方只给路径：[`code/`](../code/)、[`tools/calc/`](../tools/calc/) |

## 版面是怎么变的

- 页数：**676 → 674 → 664 → 668 → 657**（v1.26.31 → v1.27.3）。减页靠的是把"交代来源与版权"
  的文字改成细排、删掉每节重复的复算交代、收紧间距，**没有删掉任何物理内容**。
- 图数：**70 → 87 → 95**。图越换越实、越换越新（NET 装置、ITER 主机厂房、NSTX-U 加热束等）。
- 唯一的例外是 v1.27.2：补了 8 张图，页数从 664 回到 668；v1.27.3 又用版面手段收了回来。

## 质量门：六套检查，失败 0 项

| 检查 | 项数 | 管什么 |
|---|---:|---|
| `check/physics_regression.py` | 7195 | 把书里印出的数字**独立重算**一遍，并要求数字留在它自己那一节 |
| `check/full_check.py` | 49 | 结构、内容、编译、PDF、一致性、分册 |
| `check/consistency_audit.py` | 14 | 符号、图表编号、README 与正文的对应 |
| `check/ui_readability.py` | 59 | 句长、段长、颜色对比度、图文距离 |
| `check/figure_layout.py` | 16 | 插图是否越界、图注与来源是否齐全 |
| `check/pdf_outline.py` | 11 | 书签层级与印刷目录页码是否对得上 |

跑法（Windows 下先设 UTF-8，否则中文输出会崩）：

```bash
python -X utf8 check/full_check.py --no-compile
python -X utf8 check/physics_regression.py
```

## 附录体检：还想再薄，能动哪里

见 [`appendix-audit.md`](appendix-audit.md)——52 个附录逐个体检，并列出四处可动的地方
（最大的是附录 I 学习策略手册 22 页）与"动内容必须同步改断言"的代价说明。
