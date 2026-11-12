# QA.md —— 构建、门禁与可复现性

本文件说明六套门禁各管什么、在什么前提下跑、通过线是什么，以及本机没有 TeX 时哪些结论不能下。
命令与通过线以 `check/` 下的脚本为准；本文件不复制脚本实现。

## 1. 构建前提

- TeX Live 2024+ 或 MiKTeX，必须 **XeLaTeX**；构建入口 `src/plasma_physics_textbook_v2.tex`。
- 两个后处理脚本依赖 PDF 与 `.aux`：`tools/fix_outline_pages.py`、`tools/fix_section_links.py`。
- 分册切分依赖完整本 PDF 与附录 A 的起始页；切完卷一 1–207、卷二 208–610。
- 没有 XeLaTeX 的机器只能跑静态检查；不要用 `pdf/` 里的旧 PDF 冒充本轮构建产物。

## 2. 六套门禁

| 门禁 | 命令 | 管什么 | 通过线 |
|---|---|---|---|
| 物理回归 | `python3 check/physics_regression.py` | 逐项重算书中数字，并核对数字留在本节 | FAIL 0；完整构建 7195 项 |
| 全量结构 | `python3 check/full_check.py --no-compile` | 环境配对、盒、图、标签、分册、README 数字 | FAIL 0（无 TeX 时允许 C/D 两条环境 FAIL） |
| 一致性 | `python3 check/consistency_audit.py` | 术语变体、符号卡、标签唯一、附录字母、README 对齐 | FAIL 0 / WARN 0 |
| 可读性 | `python3 check/ui_readability.py` | 句长、段长、标题层级、颜色对比、盒嵌套 | FAIL 0；WARN 需逐条看 |
| 插图版式 | `python3 check/figure_layout.py` | 图宽高、越界、图注、来源、图号 | FAIL 0 |
| 书签 | `python3 check/pdf_outline.py` | 书签层级、印刷目录页码、目录链接 | FAIL 0 / WARN 0 |

## 3. 项数口径

项数会随「有没有构建产物」变化，读报告时先看它读到了什么：

- `full_check.py`：发布构建约 45 项；无 `.log` 与 `src/*.pdf` 时 32 项，因为 C（编译日志）与 D（PDF 产物）两段提前返回。
- `physics_regression.py`：有 `.aux` 时 7195 项；无 `.aux` 时附录 F 的章号归属 2 项不执行，脚本以 1 条 WARN 代替，总数 7193。
- `ui_readability.py`：发布构建约 45 项；无 PDF 时部分版面项跳过，本机静态运行 32 项。
- `figure_layout.py`：36 张图；无 `.aux` 时个别图的页码无法核对，出 WARN。
- `pdf_outline.py`：发布构建约 11 项；本机对旧 PDF 运行 10 项。

因此，README 与 `docs/overview.md` 里的「45 / 45 / 36 / 11」是**发布构建**口径；`STATUS.md` 里的数字是**本轮静态运行**口径。两者不是矛盾，前提不同。

## 4. 本机限制（2026-11）

- 没有 `xelatex` / `latexmk`，没有 `src/plasma_physics_textbook_v2.log`、`.aux`、`src/*.pdf`。
- `full_check.py` 的 C、D 两条 FAIL 是环境缺失，不是内容错误；其余 21 项 PASS。
- `pdf/` 下四册是编辑前的 v1.33.9：610 / 207 / 403 / 47 页；`figure_layout.py` 与 `pdf_outline.py` 跑的是这批旧产物。
- 本轮源码改动后未重建 PDF；页数、书签、图表页码必须等下一次构建后复核。

## 5. 失败怎么定位

1. 先看报告分类字母：A 环境配对、B 结构、C 编译日志、D PDF、E 文档一致性、F 可读性、G 盒、H 概念图、I 量纲自检、J 先猜再看、K 学习路径、L 附录导航、M 附录 B、N 解题策略、O 交叉一致性、P 代码清单、Q 分册、R 自学标准。
2. 物理回归 FAIL 先看 `[O]` 与数字所在节：多数是数字被搬到别的节，或改写时动了数值。
3. 一致性 FAIL 先看术语变体与 README 的「12 章正文 + N 个附录」「附录（N 个）」两处。
4. 书签/目录问题先确认是否跑过两个后处理；源码直接编译出来的页码会早一页。
5. 不要改 `physics_regression.py` 的期望值来消 FAIL；先判断是正文错、期望过期，还是环境缺失。

## 6. 可复现性

- 检查脚本只读 `src/`、`pdf/`、`figures/`、`docs/` 与 `check/fixtures/`，不写书稿。
- `check/fixtures/` 存可证伪样例；改检查逻辑时要同步夹具，不要只改通过线。
- 生成物（`docs/overview.png`、`docs/appendix-audit.md`、`docs/figures.md`）由 `tools/` 重跑；重跑前先构建 PDF。
- 每次发布把六套门禁的项数、FAIL/WARN 数写进 `CHANGELOG.md` 与 `STATUS.md`；只写「全部通过」不算记录。
