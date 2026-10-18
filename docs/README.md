# docs/ —— 项目一览、图册、体检表与可复跑的检查报告

这里放两类东西：**给人看的一页**（一览图、图册、体检表），和**脚本生成的检查快照**（随时可以重跑核对）。

| 文件 | 内容 |
| --- | --- |
| [`overview.md`](overview.md) ＋ [`overview.png`](overview.png) | 项目一页看完：规模、配图、交付物、六套检查、版面演进 |
| [`figures.md`](figures.md) | 图册：118 张插图按章/附录排开，含缩略图、中文说明、作者、许可与来源页 |
| [`appendix-audit.md`](appendix-audit.md) | 48 个附录逐个体检（页数、图表数）＋"还想再薄能动哪里"的四处候选 |
| [`self-study-audit.md`](self-study-audit.md) | 自学标准全量审计：12 章、48 个附录（另含配套手册 I／J／K／L 源文件）逐项对照自学要素 |

## 怎么重生成

- 图册：`python -X utf8 tools/make_figures_gallery.py --write`（读 `figures/CREDITS.md`）
- 附录体检：`python -X utf8 tools/appendix_audit.py --json` 再 `python -X utf8 tools/make_appendix_audit_doc.py --write`
- 自学标准审计：`python -X utf8 tools/selfstudy_audit.py --write`

## 判据在哪

判据写在 [`check/selfstudy_criteria.py`](../check/selfstudy_criteria.py)，三处共用同一份：
生成这份报告、全量检查里的自学标准一类、以及缺口打印。
换言之：报告里的每一项都能在检查里被重新核对，不是"看起来差不多"。

## 三分对象、三套标准

不拿一把尺子量所有东西：

- **章节**：正文 + 本章习题精选合并看（引导段、核心目标、前置依赖、例题与答案、先猜再看、
  解题策略、概念图、量纲自检、要点回顾、随堂自测、配图）；
- **工具型附录 A–I**（公式、常数、书单、术语、学习策略等；J–L 在配套手册里）：看结构、可检索性与互链；
- **专题附录 M–AZ**：看估算练习、答案、出处、适用边界、数据表与语义盒。

答案块的写法本书有三种——计算型「解答」、案例型「分析＋结论」、应用型「数值估算」——
判据同时接受，不为了统一而改写正文。
