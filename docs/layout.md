# 目录规范：每个文件该放哪里

**一句话：根目录只放四份读者最需要的东西，其余按"源 / 图 / 交付 / 文档 / 工具"五类归位。**

## 1. 顶层结构

```
plasma-physics-textbook/
├── README.md            一页看懂：是什么、4 个 PDF、怎么读、怎么建
├── CHANGELOG.md         版本记录：每一轮改了什么、为什么改
├── LICENSE              许可
├── AGENTS.md            协作者约定（含长任务进度上报）
├── .gitignore           只让交付物与源文件进仓库
├── src/                 书稿源：LaTeX 源文件，编译入口在根
├── figures/             实拍与公开许可配图 + CREDITS.md
├── pdf/                 交付 PDF（4 个，带版本号）
├── docs/                给人看的文档
├── check/               自动检查（PDF/几何/文风/一致性）
├── code/                书里引用的可运行代码
└── tools/               生成与采集脚本（calc/ 是数值脚本）
```

## 2. 各目录放什么

| 目录 | 放什么 | 不放什么 |
|---|---|---|
| `src/` | `plasma_physics_textbook_v2.tex`（主文件：导言区 + 12 章正文）、`ch*_supplement.tex`、`appendix_*.tex`、`companion_handbook.tex` | 编译产物、图片、临时脚本 |
| `figures/` | 实拍照片与公开机构示意图（jpg/png），一张图一个文件，小写下划线命名 | 自画图（那些在 `src/` 的 TikZ 里）、生成物 |
| `pdf/` | 只放**发布版** PDF：完整教材、分册卷一/卷二、配套手册 | 中间构建产物、草稿 PDF |
| `docs/` | 读者文档：`README` `overview` `figures` `appendix-audit` `self-study-audit` `palette` `style` `layout` `figures-sourcing` `review-log`；`figures-thumbs/` 是图册缩略图（生成物） | 源码、PDF |
| `check/` | 自动检查脚本与 `fixtures/` 自测数据 | 一次性分析脚本 |
| `code/` | 书里让读者复制运行的脚本（附录 G） | 内部工具 |
| `tools/` | 生成/采集/校验脚本；`calc/` 是附录数值的唯一来源（`*.py` 算、`*.txt` 存输出） | 读者的示例代码 |

## 3. 命名

| 对象 | 规则 | 例 |
|---|---|---|
| 章节补充 | `ch{N}_supplement.tex` | `ch3_supplement.tex` |
| 附录 | `appendix_{字母}_{主题}.tex`，字母重复时写全（`aa`/`az`） | `appendix_ag_equilibrium.tex` |
| 配图 | 小写英文 + 下划线，能看出内容 | `tokamak_interior.jpg` |
| 文档 | 小写英文 + 连字符 | `self-study-audit.md` |
| 缩略图 | 与配图同名 | `figures-thumbs/tokamak_interior.jpg` |

## 4. 不许进 git

编译中间物（`*.aux *.log *.out *.toc *.xdv *.fls *.fdb_latexmk *.bbl *.blg`）、
`__pycache__/`、`*.pyc`、`tmp/`、`src/*.pdf`（交付 PDF 在 `pdf/`）、
agent 运行时目录（`.cindy-worktrees/`、`.graphflow*/`、`graphflow-out/`、`.dsh-progress/`）。

## 5. 编译与产物

唯一构建入口（在 `src/` 下跑，图片用 `../figures/` 相对路径）：

```powershell
cd src
latexmk -xelatex -interaction=nonstopmode -f plasma_physics_textbook_v2.tex
```

产物 `src/plasma_physics_textbook_v2.pdf` 不入 git；要发布时复制到 `pdf/` 并打 tag。
配套手册同理（`companion_handbook.tex`）。

## 6. 新增一个附录 / 一张配图

1. 附录：新建 `src/appendix_{字母}_{主题}.tex` → 在主文件 `\input` 处按字母顺序插入 → 在 `check/` 的附录连续性与页数期望里登记 → 跑 `python -X utf8 check/full_check.py`。
2. 配图：图片放 `figures/` → 在 `figures/CREDITS.md` 追加一行（文件名 | 中文名 | 作者 | 许可 | 来源页 | 用在）→ 正文 `\includegraphics` + `\figsource{...}` → 跑 `python -X utf8 check/figure_layout.py`。
