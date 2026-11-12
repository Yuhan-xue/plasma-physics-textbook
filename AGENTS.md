# AGENTS.md

本文件是协作者（含 AI 写代理）进入本仓库的入口：先读这里，再按索引跳到对应规范。
它不替代 [`README.md`](README.md)（读者入口）与 [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)（目录与接口）。

## 0. 当前状态

- 版本 **v1.33.9**；交付 PDF 在 [`pdf/`](pdf/)，源码在 [`src/`](src/)。
- 本轮为源码编辑轮：源码已改、PDF 待构建。逐条改动见 [`ERRATA.md`](ERRATA.md)，进度见 [`STATUS.md`](STATUS.md)。
- 本机无 XeLaTeX 时不要声称 PDF 验收通过；静态门禁与不可复现项以 `STATUS.md` 为准。

## 1. 硬约束

- 只改本仓库；不要改其他仓库、不要动 `reports/`（自己的 changes 文件除外）。
- 不要 `git commit` / `push` / 切分支，除非任务明确要求。
- 不要升版本号，除非任务是发布；版本号同时出现在主文件封面、`pdftitle` 与 `README.md`。
- 不要改 `\input` 路径、label 名、构建脚本行为；不要在主文件正文新加 `\usepackage`。
- 不要改 `check/physics_regression.py` 的 7195 项期望值来「让门禁变绿」；若修改会让某项读数变化，先停下并记录。
- 物理或事实无法核实时，不改正文，写「待作者」并给出两个选项与影响。
- 预计超过约 30 秒的命令，按第 5 节用进度 CLI 包装。

## 2. 目录索引

| 路径 | 内容 | 先读 |
|---|---|---|
| [`src/plasma_physics_textbook_v2.tex`](src/plasma_physics_textbook_v2.tex) | 主文件：前言、12 章正文、工具型附录 A–H | 文件头的宏与章节顺序 |
| [`src/ch1_supplement.tex`](src/ch1_supplement.tex) … `ch12_supplement.tex` | 每章速查卡、检查清单、补充例题 | 对应章末的 `\input` |
| [`src/appendix_*.tex`](src/) | 44 个独立附录文件：专题 M–AZ 40 个 + 手册 I–L 4 个（A–H 在主文件内） | [`docs/appendix-audit.md`](docs/appendix-audit.md) |
| [`src/companion_handbook.tex`](src/companion_handbook.tex) | 配套手册（47 页） | 手册目录 |
| [`figures/`](figures/) | 插图与 [`CREDITS.md`](figures/CREDITS.md) | 图片许可 |
| [`pdf/`](pdf/) | 四个交付 PDF：完整本 610、卷一 207、卷二 403、手册 47 | [`README.md`](README.md) |
| [`docs/`](docs/) | 读者文档、架构、文风、验收、图册 | [`docs/README.md`](docs/README.md) |
| [`check/`](check/) | 六套门禁与夹具 | [`docs/QA.md`](docs/QA.md) |
| [`code/`](code/) | 书里引用的可运行脚本 | 各脚本文件头 |
| [`tools/`](tools/) | 生成、采集与后处理脚本 | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |

## 3. 构建与验收

构建（需要 TeX Live 2024+ 或 MiKTeX，必须 XeLaTeX）：

```bash
cd src && latexmk -xelatex -interaction=nonstopmode -f plasma_physics_textbook_v2.tex
```

交付前还要跑两个后处理与分册切分，命令见 [`README.md`](README.md) 的「快速开始」。
本机没有 XeLaTeX 时只能跑静态检查；不要用旧 PDF 冒充新构建。

六套门禁（最近一次发布构建全部 FAIL 0，项数见 [`README.md`](README.md)）：

```bash
python3 check/full_check.py --no-compile
python3 check/physics_regression.py
python3 check/consistency_audit.py
python3 check/ui_readability.py
python3 check/figure_layout.py
python3 check/pdf_outline.py
```

`full_check.py` 在无 `.log` / 无 `src/*.pdf` 时会报两条 FAIL（编译日志与 PDF 产物缺失）；
这是环境限制，不是内容错误，但必须在 `STATUS.md` 写明。

## 4. 文风与内容

- 文风硬指标与禁用词表：[`docs/STYLE.md`](docs/STYLE.md)；版面色板与版式：[`docs/palette.md`](docs/palette.md)、[`docs/layout.md`](docs/layout.md)。
- 不写自我削弱、套话、过程叙述与结论复述；保留适用范围、不确定性与教学装置（先猜再看、量纲自检、公式卡、概念图）。
- 公式、符号、单位、数字、label、`\ref`、`\cite`、`\input` 路径不要为了「读起来顺」而改。
- 图注只改评论性句子，不改物理内容；图片来源必须能追溯到 [`figures/CREDITS.md`](figures/CREDITS.md)。
- 勘误按 [`ERRATA.md`](ERRATA.md) 的八列记录，证据分级 A/B/C；改完给出 grep 或 diff 证据。

## 5. 长任务进度上报

预计运行超过约 30 秒的命令，一律用进度 CLI 包起来执行：

    node "$env:DSH_PROGRESS_CLI" run --task <任务id> -- <你的命令>

- CLI 会自动宣布任务开始、从命令输出里读取进度（行内出现 `12%` 或 `12/88` 都能识别）、
  退出时按退出码写终态（0 = done，非 0 = failed），不需要手写上报脚本。
- `--task` 的任务 id 必须作为完整单词出现在命令行里（面板靠它把这行进度和 job 对上）；
  id 只用字母、数字、`.`、`_`、`-`，1–40 个字符。
- 包装器表达不了的情况（脚本自己算进度），把进度事件追加到
  `$DSH_PROGRESS_DIR/<任务id>.jsonl`（一行一条 JSON，例如
  `{"v":1,"task":"<id>","state":"running","pct":42,"msg":"linking"}`），
  或用 `node "$env:DSH_PROGRESS_CLI" emit --task <id> --pct N --msg "..."` 上报、
  `done` / `fail` 收尾。
- 一个任务 id 只留一个写入者；结束必须写终态，不要把进度文件读回来——它是给人看的。

## 6. 发布

1. 六套门禁全部 FAIL 0，且 `full_check.py` 无「编译日志/PDF 缺失」外的 FAIL。
2. 版本号在同一轮改三处：主文件封面、`pdftitle`、`README.md` 的 `当前 **vX.Y.Z**`。
3. `CHANGELOG.md` 追加本轮条目（每条带文件路径或 ERRATA ID），不要改写历史条目。
4. 重出四个 PDF，更新 `pdf/` 与 `docs/overview.md` 的页数、图数、书签数。
5. 打 tag 前确认 `pdf_outline.py` FAIL 0 / WARN 0，卷一、卷二页码与完整本一致。
