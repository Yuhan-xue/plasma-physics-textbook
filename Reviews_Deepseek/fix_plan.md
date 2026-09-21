# 修复清单（可直接执行）

配套 `Reviews_Deepseek/README.md`。每条给出**文件:行号 → 当前内容 → 替换内容 → 验证方式**。

行号对应 commit `c0f24c4` 的工作树（已完成排版修复后的状态）。若行号漂移，用给出的原文定位。

---

## A. 已完成的修复（仅记录，无需再动）

| # | 文件 | 说明 | 状态 |
|---|---|---|---|
| F-1 | `src/ch10_supplement.tex:29/43/81, 475/477/479/481` | ①②③ 改用 `\ding{172..174}`；裸 α → `$\alpha$` | ✅ 已修 |
| F-2 | `src/plasma_physics_textbook_v2.tex` 导言区 + 7 个附录 `\appchapter` | 附录页眉 | ✅ 已修 |
| F-3 | 同上 | `\thechapter=\Alph{chapter}` | ✅ 已修 |
| F-4 | `:5433/5478/5508/5541/5563` | `fvextra` 断行 | ✅ 已修 |
| F-5 | `\difficulty` 宏 | 边注越界 | ✅ 已修 |
| F-6 | 导言区 | CJK 伪斜体、章首页页码、底部留白 | ✅ 已修 |

验证命令：
```bash
cd src
latexmk -xelatex -interaction=nonstopmode -f plasma_physics_textbook_v2.tex
grep -c "Missing character" plasma_physics_textbook_v2.log   # 期望 1（仅 U+000A）
grep -c "Overfull \\\\hbox" plasma_physics_textbook_v2.log   # 期望 2
pdfinfo plasma_physics_textbook_v2.pdf | grep Pages           # 期望 208
```

---

## B. P1 — 必修（导航层，直接影响自学者）

### B-1　前置知识清单：常微分方程指向错误的章节

**文件**：`src/plasma_physics_textbook_v2.tex:732`

**当前**：
```latex
常微分方程（可分离变量、线性ODE） & 第 1 章 \S 1.1 & 单粒子轨道、波动方程无法求解 \\
```

**替换为**（本书确无常微分方程复习内容，改为诚实指引）：
```latex
常微分方程（可分离变量、线性ODE） & 本书未系统复习，建议先修 & 单粒子轨道、波动方程无法求解 \\
```

**同时修改** `:853`：

**当前**：
```latex
第 3 章 & 第 1 章 \S 1.1（微分方程） & 第 2 章（准中性概念） \\
```

**替换为**：
```latex
第 3 章 & 第 1 章 \S 1.1（矢量分析） & 第 2 章（准中性概念） \\
```

**验证**：`第 1 章 §1.1` 的实际标题是「矢量分析复习」（`:1002`）；第 1 章只有 1.1 矢量分析 / 1.2 分布函数 / 1.3 傅里叶 / 1.4 δ函数，无微分方程内容。

> **备选方案**：若作者希望真的提供 ODE 复习，应在第 1 章新增一小节（如 §1.5），并同步更新 `:732` 与 `:853` 两处引用。

---

### B-2　路线 A：章节范围与文字说明矛盾

**文件**：`src/plasma_physics_textbook_v2.tex:773`

**当前**：
```latex
必读 & 第 4 章 & 只读双流体方程（\S 4.1--\S 4.3），
跳过 MHD 近似与磁流体力学 \\
```

**替换为**：
```latex
必读 & 第 4 章 & 只读 \S 4.1--\S 4.3（玻尔兹曼方程、矩方程、双流体方程组）；%
\S 4.4 广义欧姆定律可先跳过。第 5 章（磁流体力学）整章跳过 \\
```

**同时修改** `:876`（tip 盒中重复的同一表述），当前为：
```latex
第 1 章 $\to$ 第 2 章 $\to$ 第 4 章（仅 \S 4.1--\S 4.3）$\to$ 第 11 章。
```
**替换为**：
```latex
第 1 章 $\to$ 第 2 章 $\to$ 第 4 章（\S 4.1--\S 4.3：玻尔兹曼、矩方程、双流体）$\to$ 第 11 章。
```

**验证**：第 4 章实际小节为 4.1 玻尔兹曼方程、4.2 矩方程的推导、4.3 双流体方程组、4.4 广义欧姆定律。MHD 是**第 5 章**标题，不在第 4 章内。

---

### B-3　配套资源：二维码/网址不存在

**文件**：`src/plasma_physics_textbook_v2.tex:929` 附近的表格前文字

**当前**：
```latex
本书配有完整的线上资源，通过随书二维码或指定网址访问：
```

**替换为**（与 `:947` 已有的诚实说明保持一致）：
```latex
本书的配套资源情况如下。请注意：\textbf{下表所列的动图与论文导读目前尚未提供线下副本，
可离线复现的数值代码只有附录 G.1--G.3}。
```

**验证**：全仓库搜索无二维码图片、无 URL 字段；`:947` 已说明"本仓库不含其独立脚本"。

> **备选**：若资源确实存在，直接补上真实 URL。

---

## C. P2 — 应修

### C-1　Spitzer 电阻率系数与自身解析式不自洽（重要）

**文件**：`src/plasma_physics_textbook_v2.tex:2642`

**当前**：
```latex
\eta_{\parallel} \approx 5.2\times 10^{-5} \frac{Z \ln\Lambda}{T_e^{3/2}\,\mathrm{(eV)}}\,\Omega\cdot\mathrm{m}
```

**替换为**（补上 γ_E 说明，闭合矛盾）：
```latex
\eta_{\parallel} \approx 5.2\times 10^{-5} \frac{Z \ln\Lambda}{T_e^{3/2}\,\mathrm{(eV)}}\,\Omega\cdot\mathrm{m}
```

并在该式**之后**插入一段说明：
```latex
\textbf{关于这个系数的两种口径（务必分清）：} 上面写出的解析式不含动力学修正，
直接代入其常数得 $7.3\times10^{-5}$；而数值式中的 $5.2\times10^{-5}$ 是\emph{含}
Spitzer 动力学修正因子（$\gamma_E\approx0.58$ 一类，$1/1.4\approx0.71$）的标准值。
\textbf{两者相差 1.4 倍，属模型修正而非计算误差。}引用时必须连同修正一并说明，
不要把一个口径的公式配上另一个口径的数值。本节下方的算例统一采用 $5.2\times10^{-5}$。
```

**复算证据**（`Reviews_Deepseek/recompute.py`）：
```
coefficient implied by the book's DISPLAYED formula : 7.2928e-05
coefficient the book PRINTS                         : 5.2e-05
ratio                                               : 1.402
```
书中 $\eta$ 与 $\nu_{ei}$ 两式本身互相自洽（$\eta=m_e\nu_{ei}/ne^2$，比值精确 1.000000），故只需补说明，不必改解析式。

**相关**：`:2600` 已提醒过 $\gamma_E\approx0.58$，本条只是把该提醒与 §8.4 挂上钩。

---

### C-2　附录 B Spitzer 行量级与温度条件

**文件**：`src/plasma_physics_textbook_v2.tex:4870`

**当前**：
```latex
Spitzer 电阻率 & $\eta_\parallel$ & $\displaystyle\eta_\parallel = \frac{Z e^2 m_e^{1/2}}{12\pi^{3/2}\eps_0^2(k_\mathrm{B}T_e)^{3/2}}\ln\Lambda$ & $\sim 10^{-8}$ $\Omega\cdot$m & \S 8.4 \\
```

**替换为**（补下标定条件）：
```latex
Spitzer 电阻率 & $\eta_\parallel$ & $\displaystyle\eta_\parallel = \frac{Z e^2 m_e^{1/2}}{12\pi^{3/2}\eps_0^2(k_\mathrm{B}T_e)^{3/2}}\ln\Lambda$（数值式含 $\gamma_E$ 修正，见 \S 8.4） & $10^{-9}$（10\,keV）$\sim10^{-4}$（低温放电）$\Omega\cdot$m & \S 8.4 \\
```

**验证**：§8.4 算例（$T_e=10$ keV，$\ln\Lambda=17$）得 $8.8\times10^{-10}\,\Omega\cdot$m。

---

### C-3　附录 B 与第 8 章的 $\ln\Lambda$ 参数符号统一

**文件**：`src/plasma_physics_textbook_v2.tex:4867`–`:4868`

**当前**：
```latex
朗道长度 & $b_\pi$ & $\displaystyle b_\pi = \frac{e^2}{4\pi\eps_0 k_\mathrm{B}T}$ & $\sim 10^{-14}$ m & \S 8.1 \\
库仑对数 & $\ln\Lambda$ & $\displaystyle\ln\Lambda = \ln\!\left(\frac{\lambda_D}{b_\pi}\right)$ & $10\sim 25$ & \S 8.1 \\
```

**替换为**（加注与第 8 章的关系）：
```latex
朗道长度 & $b_\pi$ & $\displaystyle b_\pi = \frac{e^2}{4\pi\eps_0 k_\mathrm{B}T}$ & $\sim 10^{-14}$ m & \S 8.1 \\
库仑对数 & $\ln\Lambda$ & $\displaystyle\ln\Lambda = \ln\!\left(\frac{\lambda_D}{b_\pi}\right)$（$b_\pi$ 即 \S 8.1 的 $b_0$ 取热平均 $m_rv^2\to2k_\mathrm{B}T$ 后的形式） & $10\sim 25$ & \S 8.1 \\
```

**验证**：`b_pi/b_0 = 2.000000`；`(lambda_D/b_pi)/(3*N_D) = 1.00000000`（数值一致，仅命名不同）。

---

### C-4　23 个 TikZ 图补 caption 与 label

**现状**：全书仅 1 处 `\label{fig:`（`fig:langmuir-iv`、`fig:lawson` 除外时为 0），23 个 `tikzpicture` 中绝大多数裸嵌在正文里，无 `figure` 环境、无 caption、无 label。

**做法**（逐个图）：
```latex
\begin{figure}[htbp]
\centering
\begin{tikzpicture}
  ... 原图内容不变 ...
\end{tikzpicture}
\caption{<图的物理内容、坐标轴单位、关键标注的说明>}
\label{fig:<语义化名字>}
\end{figure}
```
并在正文首次引用处写 `如图~\ref{fig:xxx} 所示`。

**验证**：`grep -c '\\label{fig:' src/*.tex` 应从 2 增至约 25。

> 这是本清单中**工作量最大**的一项（3-4 小时），但也是提升可检索性最明显的一项。

---

### C-5　formulabox / checklistbox 标题重复

**文件**：`src/plasma_physics_textbook_v2.tex:176`、`:186`

**当前**：
```latex
\NewTColorBox{formulabox}{m m}{%
  ...
  title={核心公式速查：#1},
```
当调用处 `#1` 已是"第 3 章核心公式速查"时，标题成为"核心公式速查：第 3 章核心公式速查"。

**修复**：二选一
1. 改标题模板为 `title={#1}`，把"核心公式速查"字样移到调用处；
2. 或批量把调用处的 `#1` 改为纯主题词（如"Lawson 判据"）。

**验证**：`grep -n 'begin{formulabox}' src/*.tex` 检查每个 `#1` 是否已含"核心公式速查"。

---

## D. P3 — 打磨（可选）

| # | 内容 | 建议 |
|---|---|---|
| D-1 | 配色黑白打印不可辨（`defgreen↔historybrown` 灰度差 1） | 采用 `src/style_v2_preview.tex` 的等距灰度阶梯色板（间隔 12.0，改善 8 倍） |
| D-2 | `warnyellow` 2.04:1、`cororange` 3.79:1 低于 WCAG AA | 同上；重设计色板最低 3.13:1 |
| D-3 | 难度边注判定标准不可操作 | 每条边注补"卡住时怎么办"，替代再加抽象标准 |
| D-4 | 例题缺自检提示 | 每章选 1-2 个代表例题，补"量纲/数量级/极限/一致性"四项 |
| D-5 | 公式卡缺符号说明 | 每卡底部加"符号：…（与 §X.Y 一致）" |
| D-6 | "费米打赌"缺出处（`:2505`） | 改为有据可查的表述或补文献 |
| D-7 | TOC 点引线挤压、书签锚点警告 | 调 `\cftdotsep` / 修 hyperref 锚点 |
| D-8 | 软件清单表逐字换行（p209） | 加宽列或改 `tabularx` |

---

## E. 重设计原型落地（独立评估，非必须）

**产物**：`src/style_v2_preview.tex`（独立可编译，**未修改原稿**）

**若决定采纳**，建议按此顺序（低风险优先）：

1. **先换色板**：把 `:127-135` 的 9 个 `\definecolor` 替换为样张的 `sem1..sem8` + `inkBrand/inkBody/inkMuted`。
   → 黑白打印语义区分度立即改善 8 倍。
2. **再加冗余编码**：`tcbset` 加 `boxrule=0pt` + 各类 `borderline west={4pt}{0pt}{semN}`；`\newtcbtheorem` 的 title 前置 `\ding` 图标。
   → 这一步才是色觉无障碍的真正落点。
3. **再改层级**：章首页色带 + 节标题强调条。
4. **最后动正文流量参数**（行距/缩进/字号）——**会改变全书分页**，需整本重排并重核目录页码。

**风险提示**：
- 样张是 **11pt**，原稿 **12pt**。直接套用会改变分页。
- 封面用的是 `coverMid` 等独立色系，**只改语义色，不动封面色板**（封面已评为全书视觉质量最高的部分）。
- 需与 F-2 已修的 `\appchapter`/`\frontchapter` 机制协调，避免回退附录页眉修复。
- 建议先在 `src/` 之外试编译，逐页对比后再决定合并。

**作者需拍板的三点**（见 `ui_redesign.md` §7.3）：
1. 正文字号 11pt 还是 12pt；
2. 章首页是否采用全宽色带；
3. `sem8 提示` 的对比度取舍（维持阶梯等距 3.13:1，还是压暗至 AA 破坏等距）。

---

## F. 修复后验证清单

```bash
cd src
latexmk -xelatex -interaction=nonstopmode -f plasma_physics_textbook_v2.tex

# 1. 编译健康
grep -c "Missing character" plasma_physics_textbook_v2.log    # 期望 1
grep -c "Overfull \\\\hbox" plasma_physics_textbook_v2.log    # 期望 <= 2
grep -c "undefined" plasma_physics_textbook_v2.log            # 期望 0

# 2. 页数
pdfinfo plasma_physics_textbook_v2.pdf | grep Pages

# 3. 附录页眉（F-2 回归检查）
for p in 186 195 205; do pdftotext -layout -f $p -l $p plasma_physics_textbook_v2.pdf - | head -1; done
# 期望依次出现 "附录 A" / "附录 D" / "附录 G"

# 4. 正文章号未偏移（回归检查）
pdftotext -layout -f 20 -l 20 plasma_physics_textbook_v2.pdf - | head -1   # 期望 "Chapter 1"

# 5. 数值复算
cd .. && python Reviews_Deepseek/recompute.py    # 期望 "32 checks run, 0 outside tolerance"
```

---

## G. 流程性建议（比单条修复更重要）

### G-1　把 `Missing character` 计数纳入交付前检查

**这不是一条代码修复，而是一条流程修复。**

原稿 ①②③ 缺字形（F-1）的**根因不是疏忽，而是缺少验证步骤**：空白字形在 PDF 预览里几乎看不出来——图例只是少了圈号，文字还在，肉眼极易滑过。这正是它能逃过多轮评审的原因。

重设计样张的作者（ui-designer）在独立工作中**踩到了完全相同的坑**：初版直接写 `■◆●▲★✦▣◈`，编译立刻报 **19 处 `Missing character`**，改用 `pifont` 后归零。两次独立的同一类错误说明这是系统性问题。

**建议**：在编译流程中固定加一条断言，非零即失败：

```bash
cd src
latexmk -xelatex -interaction=nonstopmode -f plasma_physics_textbook_v2.tex
N=$(grep -c "Missing character" plasma_physics_textbook_v2.log)
if [ "$N" -gt 1 ]; then    # 1 = U+000A 换行符，属正常
  echo "FAIL: $N missing glyphs — 检查非 ASCII 装饰字符"
  grep "Missing character" plasma_physics_textbook_v2.log | sort -u
  exit 1
fi
```

**推广规则**：任何非 ASCII 装饰字符（①②③、■◆●、✓、箭头、单位符号等）在交付前，**必须以 log 的 `Missing character` 计数验证，不能只看 PDF 预览**。缺字形时优先走 `pifont`（`\ding{...}`）或中文字体，二者已验证可靠。

---

## H. 未在本清单中的事项（需单独排期）

- **第 11 章（低温工业等离子体）正文约 800 行 + 5 个补充例题的完整审查**——本轮未覆盖。这是路线 A（微电子方向）的核心章节。
- 全部 78 个例题的逐个复算（本轮抽样覆盖约 40 个）。
- 难度边注 59 处的分布合理性完整统计。
- 原子数据、经验定标与外部书目的逐条文献核查。
