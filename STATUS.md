# STATUS.md —— 当前进度

更新：2026-11-12 ｜ 分支：`editorial/unify-2026-11` ｜ 工作树：有未提交改动

## 1. 当前版本

- 版本号 **v1.33.9**，本轮未升版本。
- **源码与 PDF 已同步**：`pdf/` 下四册为 **2026-11-12** 用 TeX Live 2026 重建的 v1.33.9 产物（613 / 209 / 404 / 47 页）。
- 已在 Linux + TeX Live 2026 重建四册（613 / 209 / 404 / 47 页），源码与 PDF 同步；尚未提交或推送。
- 2026-11-12 验收：六套门禁 FAIL 0（full_check 49 项 FAIL 0／WARN 0、physics_regression 7193 项 FAIL 0、consistency 15/0、ui 32、figure 49、outline 10）。

## 2. 本轮改动

已跟踪文件改动 24 个（含 `docs/style.md` 更名为 `docs/STYLE.md`），新增 `ERRATA.md`、`STATUS.md`、`docs/ARCHITECTURE.md`、`docs/STYLE.md`、`docs/QA.md`。`docs/STYLE.md` 现为文风规范全文（三条硬线、禁用词、改写示例与自检）。

| 文件 | 改了什么 |
|---|---|
| `src/plasma_physics_textbook_v2.tex` | 热修正统一为 `σ_ve`；第 12 章相移公式修正；新增 §9.4 开尔文--亥姆霍兹不稳定性；两处大写 `v_{Te}` 与 `lnΛ` 常数因子表述修正 |
| `src/ch6_supplement.tex` | Langmuir 例题与公式卡的 `v_te` 改为 `σ_ve`，数值不变 |
| `src/appendix_l_problems.tex` | Bohm--Gross 题与答案改 `3k²σ_ve²` |
| `src/ch2_supplement.tex` | 三大判据措辞统一为 `Γ≪1`；两处断引用 `ex:low-temp-params` 改为 `ex:typical-params` |
| `src/appendix_p_lowtemp.tex` | 正文「帕申」统一为「帕邢」；节标题因回归键冲突暂留，见 ERRATA PL-06 |
| `src/ch5_supplement.tex` | 量级提示与磁冻结结论改为直接陈述，保留适用边界 |
| `src/appendix_j/x/y/al/ao_*.tex` | 去防御性写作：删「受篇幅限制」「不给这个数」「不去猜」「只做一件事」等自我削弱句 |
| `src/ch1/ch4/ch5_supplement.tex` | 删除行首孤立句号 13 处 |
| `README.md` | 按模板重写：是什么/给谁/读到什么、四个 PDF、结构与路线、规模与六套门禁、构建与检查、勘误与许可 |
| `AGENTS.md` | 补状态、硬约束、目录索引、构建验收、文风、发布清单；保留长任务进度 CLI 规则 |
| `ERRATA.md` | 新建八列台账：16 条（已修 11、待作者 1、不改 1、复核无误 3） |

详细位置、原文、证据与影响见 [`ERRATA.md`](ERRATA.md)；逐条叙述见 `reports/work/plasma-changes.md`。

## 3. 门禁状态

本轮全部为静态运行，未构建 PDF。PDF 相关两项跑的是 `pdf/` 里 2026-09-13 的 v1.33.9 旧产物，不能代表本轮源码。下表为**本机静态口径**；发布构建口径（7195 / 45 / 14 / 45 / 36 / 11）见 [`README.md`](README.md) 的双栏表与 [`docs/QA.md`](docs/QA.md) §3。

| 门禁 | 结果 | 说明 |
|---|---|---|
| `check/consistency_audit.py` | **FAIL 0 / WARN 0**，15 项 | 标签 1718 个无重复；README 12+48 与正文一致 |
| `check/full_check.py --no-compile` | **FAIL 2 / WARN 2 / PASS 21**，32 项 | 两条 FAIL 仅因缺 `src/*.log` 与 `src/*.pdf`；其余结构、盒、图、分册全过 |
| `check/physics_regression.py` | **FAIL 0 / WARN 1**，7193 项 | 7191 PASS + 1 WARN + 1 INFO；缺 `.aux` 时附录 F 章号归属 2 项未执行，完整构建后按 7195 项核对 |
| `check/ui_readability.py` | **FAIL 0 / WARN 20**，32 项 | 正文句长 p90=58 字、段长 p90=125 字，均在硬线内 |
| `check/figure_layout.py` | **FAIL 0 / WARN 2**（36 张图内 WARN 1） | 旧 PDF 产物；缺 `.aux` 时 `fig:lawson` 页码无法核对 |
| `check/pdf_outline.py` | **FAIL 0 / WARN 2**，11 项 | 613 页 / 卷 160+415、全本 575 书签 |

口径注：`consistency_audit` 的 15 项 = 14 PASS + 1 INFO（无索引条目），不是 FAIL/WARN；`physics_regression` 的 7193 项是缺 `.aux` 时少跑附录 F 的 2 项，完整构建后按 7195 项核对。

## 4. 已知问题

- **PDF 未重建**：`pdf/` 是 2026-09-13 的 v1.33.9 旧构建，本轮源码与它不同步；页数、书签、图表页码都要等下一次构建后复核，不能把旧 PDF 当作本轮产物。
- **附录 P 标题**：正文已统一「帕邢」，但节标题仍是「帕申」；改名会触发 7195 项回归的固定期望键，需作者拍板后同步脚本（PL-06）。
- **q95 工作点**：AA 的 3.00（9 MA）与 AD/AE 的 3.5（15 MA）不是矛盾，本轮不改数值；建议在 AA 加一句工作点说明（PL-08）。
- **第 3 章目标与正文**：目标提到 `J/Φ` 的推导，正文未单独展开；未在本轮改动范围内，留待作者确认是否补节或改目标。
- **第 6 章目标与正文**：目标提到 O/X 模、磁声波与阿尔芬波，正文以 R/L 波为主；同上，属内容缺口而非已确认错误。
- **无旧版可比**：本书没有旧版正文，README 已明确写「不做旧版—新版逐句对比」，本轮未伪造对比。

## 5. 下一步

1. 在装有 TeX Live 2024+ 的机器上构建，跑两个后处理与分册切分。
2. 重跑六套门禁，确认 `physics_regression` 回到 7195 项、`full_check` 的编译/PDF 两条 FAIL 消失，并复核 README 双栏表的发布构建列。
3. 作者拍板附录 P 标题（PL-06）与 AA 的 q95 工作点说明（PL-08）。
4. 复核第 3、6 章目标与正文的缺口，决定补节或收窄目标。
