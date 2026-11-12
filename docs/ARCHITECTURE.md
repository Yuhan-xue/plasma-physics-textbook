# ARCHITECTURE.md —— 目录、构建与接口

本文件说明书稿的物理布局与构建拓扑。目标是让协作者知道「东西放哪、谁读谁、改哪里会牵动哪里」，
不重排 `src/`，也不改变任何 `\input` 路径或构建入口。

## 1. 根目录

| 路径 | 角色 | 稳定性 |
|---|---|---|
| `README.md` | 读者入口：成品、结构、构建、质量、许可 | 面向读者，可改文案，不改数字口径 |
| `AGENTS.md` | 协作者硬规则与目录索引 | 规则入口，改动需同步本文件 |
| `CHANGELOG.md` | 版本记录，只追加 | 历史条目不改写 |
| `ERRATA.md` | 勘误台账，八列 | 每条带位置、证据、状态 |
| `STATUS.md` | 当前进度与门禁 | 每轮更新，不替代 CHANGELOG |
| `LICENSE` | 教材 CC BY-NC-SA 4.0；`code/` MIT | 不动 |

## 2. `src/` 的布局

| 文件/组 | 内容 | 挂载方式 |
|---|---|---|
| `plasma_physics_textbook_v2.tex` | 前言、符号说明、第 1–12 章、工具型附录 A–H、结语 | 主构建入口 |
| `ch1_supplement.tex` … `ch12_supplement.tex` | 每章速查卡、检查清单、补充例题 | 各章末 `\input{chN_supplement}` |
| `appendix_m_*.tex` … `appendix_az_*.tex` | 44 个专题/工具附录中的 M–AZ | 主文件末尾按字母顺序 `\input` |
| `appendix_i_*.tex` … `appendix_l_*.tex` | 学习策略、教材对照、题库 | 由 `companion_handbook.tex` 挂载 |
| `companion_handbook.tex` | 47 页配套手册 | 独立构建入口 |

附录 A–H 直接写在主文件里（约 5871–6848 行），不单独成文件；这是历史布局，不要为了「整齐」拆出去。
`appendix_*` 的字母顺序由主文件 `\input` 顺序决定，检查脚本也按这个顺序核对，不要调换。

## 3. 构建拓扑

```
plasma_physics_textbook_v2.tex
├── 前置：封面、版权、目录、符号说明、前言
├── 正文：ch1 … ch12（每章末 \input{chN_supplement}）
├── 工具附录：A–H（内联）
└── 专题附录：M … AZ（逐个 \input）
        ↓ XeLaTeX × latexmk
src/plasma_physics_textbook_v2.pdf（610 页）
        ↓ fix_outline_pages.py / fix_section_links.py
pdf/plasma_physics_textbook_v2.pdf
        ├── 按附录 A 起始页切分 → 卷一 1–207、卷二 208–610
        └── 书签与印刷目录页码随 PDF 一起保留

companion_handbook.tex → pdf/plasma_physics_handbook.pdf（47 页）
```

三个派生事实：

- 卷一、卷二的印刷页码与完整本一致，不是重新编号。
- `pdf_outline.py` 检查的是 PDF 书签与印刷目录，不是 `.tex` 源码；必须在构建与两个后处理之后跑。
- `full_check.py` 的页数/版本一致性读取 `src/*.pdf` 与 `.log`；无构建时这两条会 FAIL，属环境限制。

## 4. 接口与不变量

- **标签**：`\label{...}` 全书唯一（当前 1718 个），`\ref`/`\eqref`/`\cref` 只指向存在的标签；不要重命名已有 label。
- **输入**：`\input{...}` 使用相对 `src/` 的文件名；不要改路径、不要加子目录。
- **图**：源码引用 `../figures/<name>`；图片版权与来源登记在 `figures/CREDITS.md`，新增图必须同时登记。
- **数字**：正文数字与 `check/physics_regression.py` 的期望值绑定；改数字先看回归项，不要改期望值来迁就正文。
- **盒子**：`theorem`、`insight`、`warning`、`tip`、`predict`、`strategy` 等语义盒都是两参数调用（标题 + 唯一标签）；不要裸写 `\begin{tip}`。
- **章末装置**：12 章各一条 `\item \textbf{量纲自检}`、各一张 `fig:cm-chN`、各一只「解题策略」盒；新增内容不要破坏这些计数。
- **自测盒**：全书 18 只 `miniquiz`（16 只带答案 + 2 只导读盒）；新增一节不要顺手加盒，除非同步检查脚本。

## 5. 生成物与脚本

| 生成物 | 生成命令 | 读什么 |
|---|---|---|
| `docs/overview.png` | `python3 tools/overview_fig.py` | `pdf/`、`figures/CREDITS.md`、`check/` 项数 |
| `docs/appendix-audit.md` | `python3 tools/appendix_audit.py --json` + `make_appendix_audit_doc.py --write` | 48 个附录的页数、图表数 |
| `docs/figures.md` | `python3 tools/make_figures_gallery.py --write` | `figures/CREDITS.md` |
| `docs/self-study-audit.md` | `python3 tools/selfstudy_audit.py --write` | `check/selfstudy_criteria.py` |

`check/` 只读源码与 PDF，不写书稿；`tools/` 可以生成文档，但不应改 `src/`。两处的职责不要混。

## 6. 改内容时的最短路径

1. 先定位章节：正文在主文件，章末补充在 `chN_supplement.tex`，专题在 `appendix_*.tex`。
2. 公式、符号、数字、label、`\input` 路径不动；只改叙述时，保持 `\ref` 目标不变。
3. 新增一节：插在既有 `\section` 之间，用现有语义盒；不要新增 `\usepackage`，不要新增章末量纲自检条目。
4. 改完跑 `consistency_audit.py` 与 `full_check.py --no-compile`；改数字再跑 `physics_regression.py`。
5. 记录到 `ERRATA.md`（若属勘误）或 `CHANGELOG.md`（若属版本改动），并在 `STATUS.md` 更新门禁状态。
