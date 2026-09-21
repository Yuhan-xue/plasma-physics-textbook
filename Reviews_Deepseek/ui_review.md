# t1 — 《等离子体物理自学教材》UI／排版与视觉设计审查

- 审查对象：`src/plasma_physics_textbook_v2.tex`（5618 行）+ `src/ch*_supplement.tex`
- 最终证据：`pdf/plasma_physics_textbook_v2.pdf`（210 页 A4，595.28×841.89pt）
- 方法：`pdftoppm` 全页 100dpi 灰度栅格化 + 关键页 200–400dpi 实际观察（`read_image`）、`pdftotext -layout`、日志告警解析、自写像素级版面测量脚本
- 本轮为只读审查，未修改 `src/` 下任何文件

## 0. 证据生成方式（可复现）

```bash
pdftoppm -r 100 -gray -png pdf/plasma_physics_textbook_v2.pdf reviews_tmp/hdr/g   # 210 页栅格
pdftoppm -f <n> -l <n> -r 300 -png pdf/plasma_physics_textbook_v2.pdf out        # 单页高分辨率
pdftotext -layout pdf/plasma_physics_textbook_v2.pdf full.txt
```
像素→物理单位换算：100dpi ⇒ 1cm = 39.3701px。由 `geometry{left=2.1cm,right=2.3cm}` 得正文右边界 18.70cm，边注盒 18.90–20.50cm，纸张宽 21.00cm。
逐页扫描「右边界墨迹最右列」的脚本见 §附录；结论：**50 页**墨迹越过正文右边界，其中 **20 页**越过边注盒右边界（20.62cm）。

---

## 1. 结论摘要

| 严重度 | 数量 | 主题 |
|---|---|---|
| **P1** | 4 | Lawson 图 ①②③ 全缺失致图意不可读；附录全书页眉错误 + 章号残留；附录编号体系崩坏（`.7`/`定义 .16`）；正文中 `Bmax`/`θ > θm` 等原始 TeX 记号泄漏到版面 |
| **P2** | 8 | 20 页边注越出边注盒；`formulabox`/`checklistbox` 标题重复；23 个 TikZ 图**无任何 caption**；图内标注字号过小；`\difficulty` 边注脱离所属小节；章首页无页码；CJK 伪斜体；软件清单表单元格逐字换行 |
| **P3** | 6 | 配色黑白打印不可辨（26 组成对灰度过近）；红绿色盲混淆 5 组；`cororange`/`warnyellow` 对白对比度不足；30 页空白页脚区；TOC 无二级缩进层次 |

整体评价：**版式骨架（geometry、titleformat、tcolorbox 统一 tcbset）是有设计意图且执行到位的**，封面是全书视觉质量最高的部分。主要缺陷集中在**附录段落的编号/页眉机制**与**图的图文配套**两处系统性问题，而非零散瑕疵。

---

## 2. P1 级发现

### P1-1　Lawson 参数空间图：①②③ 三个判据标记全部渲染为空白，图的语义失效

- **定位**：`src/ch10_supplement.tex:475,477,479`（TikZ 节点）、`:481`（caption）；PDF **第 126 页**
- **问题**：`\node[...] {①不满足}` 等 3 处节点与 caption 内 3 处 ①(U+2460)②(U+2461)③(U+2462)，所用字体为 `lmroman10-regular`，该字体无此字形。PDF 中三个标记位置**完全空白**，图例只剩裸文字「不满足／临界／满足」。
- **复现证据**：
  - 日志：`Missing character: There is no ① (U+2460) in font [lmroman10-regular]...` ×3，②③ 各 ×3，共 9 处（附件文件名 `ch10_supplement.tex`）
  - `pdftotext -f 126 -l 126` 输出中该行呈现为 U+FFFD 替换字符：
    ```
    22   NIF   点火区（阈值线上方）
                     JET  D-D 点火线
    21            � 满足
    ```
  - 第 126 页 300dpi 实际观察：图例三项中「满足」前**无任何可见符号**，与 caption 文案「①在阈值线下（不满足）、②恰在线上（临界）、③在线上方（满足）」无法对应。
- **影响**：该图是第 10 章 Lawson 判据的核心教学图，读者无法把三个状态点与 caption 的文字解释对应起来 ⇒ 知识点不可传达。
- **修复建议**（`\textcircled` 走中文字体，或改用 TikZ 原生绘制，二者择一）：

```latex
% 方案 A（最小改动）：改用 ctex 已加载的中文字体承载圈号
\usepackage{pifont}            % 导言区
\newcommand{\cnum}[1]{\ding{\numexpr171+#1\relax}}  % \cnum1=\ding{172}=①
% ch10_supplement.tex:475/477/479
\node[red,right,font=\footnotesize] at (3.05,1.75) {\cnum1 不满足};
\node[orange,right,font=\footnotesize] at (3.05,2.4) {\cnum2 临界};
\node[green!50!black,right,font=\footnotesize] at (3.05,3.05) {\cnum3 满足};
% caption 内同理：\cnum1 在阈值线下（不满足）……

% 方案 B（跨字体最稳）：让圈号走中文字体
\newcommand{\cnum}[1]{{\CJKfamily{zhkai}\char\numexpr"245F+#1\relax}}
```
  另：`:481` caption 中**裸 α**（"仅 α 留在等离子体中"）同属 `Missing character: no α (U+03B1) in lmroman12-regular`，应写 `$\alpha$`。

### P1-2　附录 7 章全部顶着「Chapter 12. 等离子体基础诊断」页眉，且无页码

- **定位**：`src/plasma_physics_textbook_v2.tex:4569`（`\appendix`）、`4574/4765/4855/5032/5135/5217/5342`（`\chapter*{附录X：…}`）；PDF **第 186–210 页**（共 25 页）
- **问题**：`\leftmark` 在 `\chapter*` 后不更新，fancyhdr 的 `\fancyhead[L]{\nouppercase{\small\itshape\leftmark}}` 于是永远保留最后一个编号章（第 12 章）的名字。25 页附录全部页眉错误。同时页脚区墨迹为 0（无页码）。
- **复现证据**：全页栅格扫描页眉带墨迹，第 185 页之后页眉 ink 持续 >2000px，`pdftotext -layout -f 186`/`-f 190`/`-f 196`/`-f 210` 首行**全部**为
  ```
  Chapter 12. 等离子体基础诊断      186
  ```
  而章首页（19/29/39/49/60/71/83/93/104/115/128/156）页眉墨迹 = 0（`\thispagestyle{plain}`，本属 book 正常行为）。
  边注：旧评审提到的「第 189 页附录页眉有问题」经实测**确认成立，且范围是 186–210 全段**，不止 189 页。
- **修复建议**（对 `\chapter*` 手工设置 `\markboth`，并给附录页配页眉与页码）：

```latex
% 导言区：定义带书签/页眉的附录章命令
\newcommand{\appchapter}[1]{%
  \chapter*{#1}%
  \addcontentsline{toc}{chapter}{#1}%
  \markboth{#1}{#1}%          % 关键：让 \leftmark 跟随附录
  \thispagestyle{fancy}%
}
% 正文 4569 行之后逐处替换：
\appchapter{附录A：矢量与张量公式汇总}
\appchapter{附录B：等离子体特征参数速查表}
% ... G 同理
% 并在 \appendix 后重设章节计数器，避免残留章号
\renewcommand{\thechapter}{\Alph{chapter}}
```

### P1-3　附录编号体系崩坏：小节显示为 `.1 .7`，tcolorbox 环境显示为「定义 .16」「定理 .6」

- **定位**：`src/plasma_physics_textbook_v2.tex:4579` 及附录内所有 `\section{...}`、`4583` 起所有 `definition/theorem/...` 环境；PDF **第 185–210 页**
- **问题**：`\appendix` 后紧跟 `\chapter*`（**不**递增 chapter 计数器），但附录内仍用 `\section{}`。`\thesection` = `\thechapter.\arabic{section}`，而 `\thechapter` 此刻为空 ⇒ 编号渲染成孤立的 `.1`、`.7`。「本章汇总…」的 `\section{直角坐标系 $(x,y,z)$}` 在 PDF 中显示为 **`.1 直角坐标系`**。tcolorbox 用 `number within=chapter`，同样产生 **`定义 .16`**、**`定理 .6`**。
- **复现证据**：`pdftotext -layout` 摘要（页 → 首行）：
  ```
  --p185: 附录 A：矢量与张量公式汇总 | … | .1 直角坐标系 (x, y, z)
  --p186: … | 定义 .16: 旋度 |
  --p187: … | 定理 .6: 基本恒等式 |
  --p189: 附录 B：… | .6 参数总表
  --p191: 附录 C：… | .7 C.1 广义欧姆定律的完整推导      ← 与正文手写的 C.1 重复编号
  --p192: … | .8 C.2 Landau 围道积分的详细推导
  --p195: 附录 D：… | .11 基础物理常数
  --p197: 附录 E：… | .16 专业级（英文专著）
  --p203: 附录 G：… | .24 G.1 单粒子轨道追踪
  ```
  `.24 G.1` 一行同时出现两套编号，直接暴露机制。
- **影响**：附录是全书的查阅入口（速查表、公式汇总、易错点汇总）。编号 `.16 旋度`、`定义 .16` 会让读者以为引用编号丢失或排印损坏；「附录 C 的 C.1」在目录里被写成 `.7`，与正文 §4.4 的交叉引用口径不一致。
- **修复建议**：附录内改用不带编号的层级，或显式给附录章编号。

```latex
% 方案 A（推荐，改动最小）：附录内小节不编号，保留标题层级
% 在 \appendix 后加：
\setcounter{secnumdepth}{0}%   附录内 \section 不再产生 ".1"
% 但 tcolorbox 的 {定义}{...} 编号仍受 number within=chapter 影响，需一并处理：

% 方案 B（完整）：让附录章真正编号
\appendix
\renewcommand{\thechapter}{\Alph{chapter}}
\setcounter{chapter}{0}
\chapter{矢量与张量公式汇总}      % 去掉星号 -> 自动编号 A/B/C...
% 并删掉紧邻的 \addcontentsline（\chapter 已自动进目录）
% 章节标题文案 "附录A：" 可交由
\titleformat{\chapter}[display]{...}{\filleft\fontsize{64}{64}\selectfont
  \color{coverMid!30}附录\thechapter}{-2.4ex}{\Huge\heiti\filright}
```
  这样 `\thechapter` = A/B/C…，`\section` 得 `A.1`，tcolorbox 得 `定义 A.1`，与正文体例一致。

### P1-4　正文中泄漏原始 TeX 记号：`Bmax`、`θ > θm` 等无数学定界符的裸文本

- **定位**：`src/plasma_physics_textbook_v2.tex` 第 47 页对应内容（磁镜损失锥图）；PDF **第 47 页**
- **问题**：图内标注直接以文本模式书写数学符号而未加 `$...$`，PDF 中呈现为 `Bmax`（应为 $B_{\max}$）、`θ > θm`（应为 $\theta>\theta_m$）。下角标丢失、希腊字母以正体西文呈现。
- **复现证据**：`pdftotext -layout -f 47 -l 47`：
  ```
  --p47: Chapter 3. 单粒子运动 … | Bmax | θm   逃逸 | θ > θm  | θ < θm
  ```
  第 47 页 300dpi 实际观察（裁 `-x 140 -y 120 -W 1400 -H 620`）：图内四个标注确为等宽正体纯文本 `Bmax` / `θm` / `θ > θm` / `θ < θm`，与同页正文的数学字体（Computer Modern 斜体 $B$、$\theta$）**明显不一致**。
- **影响**：物理书里下标是本征信息。`Bmax` 与 $B_{\max}$ 在教学上不等价（读者会误读为变量名），且视觉上与全书公式风格断裂。
- **修复建议**（`tikzpicture` 内节点一律进数学模式）：
```latex
% 错误
\node at (0,1) {Bmax};
\node at (1,0) {θ > θm};
% 正确
\node at (0,1) {$B_{\max}$};
\node at (1,0) {$\theta>\theta_m$};
% 中文与公式混排时，中文留在文本模式、公式包 $...$：
\node at (2,0) {逃逸区（$\theta>\theta_m$）};
```
  同类风险点建议全书 `tikzpicture` 内做一次 `\{\s*[\u0370-\u03ff]` 与 `\b[A-Z][a-z]+max\b` 的排查。

---

## 3. P2 级发现

### P2-1　20 页的难度边注越出边注盒右边界（最远 20.68cm / 纸宽 21.00cm）

- **定位**：`src/plasma_physics_textbook_v2.tex:213`（`\difficulty` 宏）、`:94-95`（`marginparwidth`/`marginparsep`）；PDF **第 22, 23, 31, 32, 40, 63, 95, 129, 137, 139, 141, 150, 152, 161, 164, 168, 180, 181, 207 页**（另 31 页落在：
- **问题**：`\difficulty` 定义为
  ```latex
  \newcommand{\difficulty}[1]{\marginnote{\hfill\textcolor{coverMid!80}{\rule[-0.1em]{2pt}{0.9em}}\,\textcolor{coverMid}{\footnotesize\textbf{[#1]}}}}
  ```
  其前导 2pt 竖线 `\rule` 在 `\hfill` 推动下贴在边注盒最右，加上 tcolorbox 的 `right=5pt` 与 `enhanced` 外框，墨迹停在 **20.62cm**，比边注盒右边界 20.50cm 多出 0.12cm，比正文右边界多 1.92cm。第 207 页（附录 G 代码）更达 **20.68cm**。虽仍在纸内（未裁切），但已吃掉 2.3cm 页边距中的 1.9cm，视觉上「贴边」。
- **复现证据**：逐页像素扫描 `mask[:, int(20.50cm):]` 得墨迹行段，13 个样本页**全部**呈现 `y` 方向 15px（≈1.5 行高）的孤立墨迹块，最右列恒为 `20.62cm`（207 页为 20.68cm）。400dpi 放大裁剪确认该墨迹形状为**一段竖线**，即 `\rule[-0.1em]{2pt}{0.9em}`，而非文字。
- **修复建议**（去掉多余的 `\hfill`，并把竖线收到盒内）：
```latex
% 导言区：略增边注宽度，并让内容左对齐而非右顶
\setlength{\marginparwidth}{1.6cm}
\setlength{\marginparsep}{0.2cm}
\renewcommand{\difficulty}[1]{%
  \marginnote{\raggedright
    \textcolor{coverMid!80}{\rule[-0.1em]{2pt}{0.9em}}\,
    \textcolor{coverMid}{\footnotesize\textbf{[#1]}}}}
% 若仍需竖线贴左、文字贴右，改为：
% \marginnote{\textcolor{coverMid!80}{\rule[-0.1em]{2pt}{0.9em}}\,\textcolor{coverMid}{\footnotesize\textbf{[#1]}}\hfill}
% 并同步把 geometry right 由 2.3cm 放宽到 2.6cm，令 20.62cm 回到盒内。
```
  另建议 `\difficulty` 用 `\marginnote[..]` 的可选垂直偏移，避免与 tcolorbox 标题行挤在同一基线上（见 P2-4）。

### P2-2　`formulabox` / `checklistbox` 标题重复，出现「核心公式速查：第 3 章核心公式速查」

- **定位**：`src/plasma_physics_textbook_v2.tex:176-180`（`\NewTColorBox{formulabox}{m m}`）、`:186-190`（`checklistbox`）
- **问题**：宏定义把标题固定为 `title={核心公式速查：#1}`，而调用处传入的 `#1` 本身已经含「第 N 章核心公式速查」，于是拼接成重复文案。
- **复现证据**：`pdftotext -layout` 摘要可见 3 处重复标题，PDF 实际观察第 48、92、149 页：
  ```
  --p48 : 核心公式速查：第 3 章核心公式速查
  --p92 : 核心公式速查：第 7 章核心公式速查
  --p149: 核心公式速查：Townsend 放电与帕邢定律   ← 正确形态（#1 为具体主题）
  ```
  第 48 页 300dpi 裁剪（`-x 100 -y 200 -W 2380 -H 900`）确认标题栏文字为「核心公式速查：第 3 章核心公式速查」，章号重复。
- **修复建议**（把固定前缀交给调用方，或做调用处去重）：
```latex
% 方案 A：宏不再加前缀，调用方写全
\NewTColorBox{formulabox}{m m}{%
  colback=blue!5,colframe=blue!60,fonttitle=\bfseries,breakable,
  title={#1}, phantomlabel={frml:#2},
  left=4pt,right=4pt,top=4pt,bottom=4pt}
% 调用：\begin{formulabox}{第 3 章核心公式速查}{ch3}

% 方案 B（不改调用方）：宏内先剥掉可能的重复词
\NewTColorBox{formulabox}{m m}{%
  ...,
  title={#1\ifnum\pdfstrcmp{#1}{第 \thechapter 章核心公式速查}=0\else\fi},
  ...}
% 更稳妥：调用处统一改为 \begin{formulabox}{第 \thechapter 章核心公式速查}{ch3}
% 并让宏固定前缀，二者只保留一处。
```

### P2-3　23 个 `tikzpicture` 全部没有 caption，图与正文无法交叉引用

- **定位**：全书 `\begin{tikzpicture}` 共 23 处（`src/ch*_supplement.tex` 为主）；仅 `ch10_supplement.tex:483` 处有 `\label{fig:lawson}`
- **问题**：绝大多数图直接以裸 `tikzpicture` 排版，**不在 `figure` 环境内、无 `\caption`、无 `\label`**。全文 `\ref` 共 59 处，但**没有一处引用图**（`fig:` 前缀标签仅 1 个）。
- **复现证据**：
  - 结构统计：23 个 `tikzpicture` vs. `\caption` 数量匹配得上表/框但配不上图；`\label{fig:` 全书仅 1 处。
  - 实际观察第 47、63、80、91、101、112、126、158 页：图形完好，但图下方**无任何图注文字**。第 47 页损失锥图、第 80 页波分类树图、第 112 页不稳定性分类图均为高信息量图，读者无从得知「图 N」编号，正文也无法回引。
- **影响**：自学教材中图是主要理解通道；无编号收敛为「孤图」，破坏 navigate-ability（与 TOC/交叉引用体系不一致——表格有 `表2/表3` 编号，图却没有）。
- **修复建议**（统一包进 `figure` 并加 caption/label；全书 23 处批量处理）：
```latex
\begin{figure}[htbp]
  \centering
  \begin{tikzpicture}[...]
    ...
  \end{tikzpicture}
  \caption{磁镜损失锥：$\theta<\theta_m$ 的粒子被约束，$\theta>\theta_m$ 逃逸。}
  \label{fig:loss-cone}
\end{figure}
% 正文引用：如图~\ref{fig:loss-cone} 所示……
% 同时导言区已设 \captionsetup[figure]{position=bottom}，行为正确，无需改动。
```

### P2-4　`\difficulty` 边注的位置漂移：从所属小节标题跳到正文首段行

- **定位**：`src/plasma_physics_textbook_v2.tex:955, 1034, 1123, …`（`\section{...}\difficulty{基础}` 形式，共 58 处）；PDF 第 20 页示例
- **问题**：`\difficulty` 写在 `\section{...}` **之后**，`\marginnote` 锚定的是「当前垂直位置」。由于 `titlespacing*{\section}{0pt}{12pt...}{6pt...}` 的段后间距，边注实际落在**标题下方 12pt 的正文首行**基线上，而非标题行。第 20 页可见 `[基础]` 与「本节快速复习电动力学…」同一行，而标题「1.1 矢量分析复习」右侧为空。
- **复现证据**：`pdftotext -layout -f 20 -l 20`：
  ```
  1.1 矢量分析复习
     本节快速复习电动力学中已学过的矢量分析内容。若读者对此尚不熟悉，建议先系统        [基础]
  学习矢量分析教材。
  ```
  第 20 页 300dpi 裁剪（`-x 150 -y 240 -W 2350 -H 200`）确认 `[基础]` 的垂直位置与正文首行对齐、与标题行错位。
- **影响**：全书 58 处难度标记全部错位一行，读者需回溯判断它属于哪个小节；在跨页标题处更会与上一页内容混同。
- **修复建议**（把标记并入标题行，用 `\titleformat` 的末参数承载）：
```latex
% 方案 A：用 titlesec 的「标题后」参数承载难度标记
\newcommand{\currdiff}{}
\titleformat{\section}
  {\normalfont\Large\heiti\color{coverMid}}
  {\thesection}{0.9em}{}
  [{\ifx\currdiff\empty\else\marginnote{\raggedright\currdiff}\fi}]
\newcommand{\difficulty}[1]{\gdef\currdiff{\textcolor{coverMid!80}{\rule[-0.1em]{2pt}{0.9em}}\,\textcolor{coverMid}{\footnotesize\textbf{[#1]}}}}
% 调用保持 \section{矢量分析复习}\difficulty{基础}，标记即出现在标题行右侧

% 方案 B（更简单）：把标记放在标题之前
\difficulty{基础}\section{矢量分析复习}
% 注意此时 \marginnote 锚在标题上一行的位置，仍需配合 \marginnote[<offset>] 微调
```

### P2-5　12 个章首页无页码（book 默认 `plain` 未定制）

- **定位**：`src/plasma_physics_textbook_v2.tex:105-109`（fancyhdr 设置）；PDF 第 19, 29, 39, 49, 60, 71, 83, 93, 104, 115, 128, 156 页
- **问题**：全书 `\pagestyle{fancy}`，但 `\chapter` 自动下发 `\thispagestyle{plain}`，而文档**未重定义 `plain`**。章首页因此既无页眉也无页码。
- **复现证据**：全页栅格扫描：上述 12 页页眉带墨迹 = 0，页脚带墨迹亦接近 0（48–2562px，多为图形溢出而非页码）。对比第 20 页页眉 ink=1642、footer ink=1932（有页码 20）。
- **影响**：翻到章首页时丢失定位信息；对 210 页的自学教材，读者频繁凭页码回查目录。
- **修复建议**（重定义 plain，章首页只保留页眉线与页码、不要页眉文字）：
```latex
\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyhead[R]{\small\bfseries\thepage}%
  \renewcommand{\headrulewidth}{0.3pt}%
}
% 若希望章首页连页眉一起省略、仅要页码居中：
\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyfoot[C]{\small\bfseries\thepage}%
  \renewcommand{\headrulewidth}{0pt}%
}
```

### P2-6　CJK 伪斜体：`TU/SimSun(0)/b/it` 未定义，SimSun 被施以合成斜体

- **定位**：`src/plasma_physics_textbook_v2.tex:123-125`（epigraph 用 `\kaishu`）、`:215-219`（`chapterintro` 环境用 `\itshape`）；日志 `LaTeX Font Warning: Font shape TU/SimSun(0)/b/it undefined ... on input line 234`
- **问题**：`chapterintro` 定义为 `\begin{quotation}\itshape\small`，`\itshape` 在中文段落上会落到 `SimSun/m/it`。SimSun 无真斜体，fontspec 走**合成倾斜**，中文笔画出现不自然的机械倾斜。
- **复现证据**：
  - 日志：`LaTeX Font Warning: Font shape 'TU/SimSun(0)/b/it' undefined ... using 'TU/SimSun(0)/b/n' instead on input line 234`；并由 `[]|\TU/SimSun(0)/m/it/12 涓?…` 可见 `m/it` 形态确实被使用。
  - 日志：`LaTeX Font Warning: Some font shapes were not available, defaults substituted.`
  - 另：`Could not resolve font "FangSong/B"`、`"FangSong/I"`、`"FangSong/BI"`、`"KaiTi/B"`、`"KaiTi/I"`、`"KaiTi/BI"`、`"SimHei/I"`、`"SimHei/B"`、`"SimSun/BI"` 共 9 类变体未解析。
- **影响**：中文正文（本书主体）出现合成斜体，与楷体引语、黑体标题混排时字型层级混乱；`\textbf` 与 `\emph` 在中文环境下的表现不可预测。
- **修复建议**（中文段落禁用 `\itshape`，改用楷体/颜色做层级区分；文件本身 123–125 行已有正确模式可复用）：
```latex
% chapterintro：中文不用 \itshape，改用楷体（KaiTi 无斜体，这也正是 123 行的注释所说）
\newenvironment{chapterintro}{%
  \begin{quotation}\small\kaishu
}{%
  \end{quotation}%
}
% 全局兜底：中文斜体请求直接吃掉
\usepackage{xpatch}
\xpatchcmd{\itshape}{\itshape}{}{}{}   % 或在 ctex 下设置
% 更规范：为 CJK 家族显式声明无斜体变体
\setCJKfamilyfont{zhsong}{SimSun}[ItalicFont=KaiTi,BoldFont=SimHei]
```

### P2-7　附录 G 软件清单表：单元格逐字换行，可读性受损

- **定位**：`src/plasma_physics_textbook_v2.tex` 第 209 页对应表格「等离子体数值模拟软件清单」（`表8`）；PDF **第 209 页**
- **问题**：表格第二列宽度不足，中文长词被逐字断行，且 `\AtBeginEnvironment{tabular}{...\small}` 收缩字号后仍不够。实际观察为：
  ```
  PIC 粒子模拟   EPOCH            英 国 华 威 大 学 开 发， 开 源， 支 持
                                  1D/
  ```
  「英 国 华 威 大 学 开 发」每个汉字之间出现空格（justify 撑开的字间距），几乎不可读。
- **复现证据**：`pdftotext -layout -f 209`；第 209 页 300dpi 裁剪（`-x 850 -y 900 -W 1200 -H 500`）实际观察确认字间距被拉开。
- **修复建议**（给长文本列指定固定宽度并左对齐，允许换行但不 justify）：
```latex
% 用 p{} 列承载说明文字，并禁止两端对齐拉开字距
\begin{tabular}{@{}p{1.6cm} p{2.6cm} p{7.4cm}@{}}
\toprule
\textbf{类型} & \textbf{软件名称} & \textbf{特点与适用场景} \\
\midrule
PIC 粒子模拟 & EPOCH & 英国华威大学开发，开源，支持 1D/2D/3D…… \\
\bottomrule
\end{tabular}
% 若整表仍超宽，改用 longtable + \small 或把「类型/软件名称」合并为一列
% 全局已设 \setlength{\tabcolsep}{5pt}，可临时收紧到 3pt。
```
  同时建议为表格加全局 `\raggedright\arraybackslash`（`array` 包已加载）以消除 justification 字距。

### P2-8　30 页出现「页脚空白带」但页眉正常 —— 页码带检测异常，需确认页码连续性

- **定位**：PDF 第 7, 10, 11, 12, 15, 18, 21, 23, 25, 26, 27, 30, 31, 33, 34, 36, 37, 38, 44, 45, 47, 48, 50, 54, 55, 56, 59, 62, 65, 68, 69, 70, 73, 78, 79, 80, 82, 84, 86, 87, 89, 90, 91, 92, 95, 101, 103, 105, 106, 111, 112, 114, 116, 118, 127, 130, 131, 137, 144, 147, 150, 154, 155, 157, 159, 161, 162, 165, 168, 170, 173, 176, 177, 181, 182, 183, 186, 187, 188, 190, 193, 194, 196, 198, 199, 202, 207, 209, 210 页
- **问题**：`\fancyhead[R]{\small\bfseries\thepage}` 把页码放在**页眉右**而非页脚。这是设计选择本身合理，但 `footskip=1.0cm` 与 `geometry bottom=2.1cm` 组合下，正文末行离页脚很近，而页脚区完全空白——版面底部约 1cm 的空白与顶部 2.1cm 边距形成不对称。
- **复现证据**：89 页页脚带（底部 2.1cm+20px 内）墨迹为 0；同时 `footskip=1.0cm` 设定了页脚基线，说明作者预留了页脚空间但未使用。
- **修复建议**（二选一，明确设计意图）：
```latex
% 选择 1：把页码从页眉移到页脚居中（更常见于教材），并收紧 bottom
\fancyhead[R]{}
\fancyfoot[C]{\small\bfseries\thepage}
\geometry{left=2.1cm,right=2.3cm,top=2.1cm,bottom=1.8cm,footskip=0.9cm}

% 选择 2：保留页眉页码，把 footskip 与 bottom 收回，消除底部空带
\geometry{left=2.1cm,right=2.3cm,top=2.1cm,bottom=1.7cm,footskip=0.4cm}
```
  **本条为打磨项（P3 倾向）**，列在 P2 是因为「版面上下不对称」在 210 页尺度上是可感知的整体观感问题。

---

## 4. P3 级发现

### P3-1　配色系统：9 种语义色在黑白打印下 26 组灰度过近，其中 4 组几乎相同

- **定位**：`src/plasma_physics_textbook_v2.tex:127-135`
- **问题**：语义色在**彩屏**上区分良好，但转为灰度后大量塌缩。以 NTSC 亮度 `0.299R+0.587G+0.114B` 计算：

| 色名 | RGB | 对白对比度 (WCAG) | 灰度亮度 | Δ≤25 的配对 |
|---|---|---|---|---|
| thmblue | (41,98,255) | 4.90:1 | 98.9 | vs defgreen 6、vs historybrown 5、**vs reviewblue 2**、vs tipgreen 8 |
| defgreen | (46,125,50) | 5.13:1 | 92.8 | **vs historybrown 1**、vs green60 3、vs reviewblue 8 |
| cororange | (230,81,0) | 3.79:1 | 116.3 | vs reviewblue 16、vs tipgreen 9、vs historybrown 22 |
| exampurple | (106,27,154) | 9.39:1 | 65.1 | vs coverMid 21、vs insightred 9 |
| insightred | (183,28,28) | 6.57:1 | 74.3 | vs exampurple 9、vs historybrown 20 |
| warnyellow | (255,160,0) | **2.04:1** | 170.2 | （亮度独高，无 Δ≤25 配对） |
| historybrown | (121,85,72) | 6.55:1 | 94.3 | **vs defgreen 1**、vs green60 4、vs thmblue 5 |
| reviewblue | (25,118,210) | 4.60:1 | 100.7 | **vs thmblue 2**、vs tipgreen 6 |
| tipgreen | (56,142,60) | 4.12:1 | 106.9 | vs reviewblue 6、vs thmblue 8、vs historybrown 13 |

  「Δ≤25」配对共 **26 组**，其中 `defgreen↔historybrown` 灰度差仅 **1**、`thmblue↔reviewblue` 仅 **2**、`defgreen↔green60` 仅 **3**——黑白打印后完全无法区分。
- **复现证据**：`reviews_tmp/color_out.txt` §「pairwise grayscale separation」全表；栅格页面对彩色描边的采样（`frame_scan`）确认帧色即为上表 RGB。
- **影响**：本书面向自学，打印场景普遍。8 类 tcolorbox（定理/定义/推论/例题/物理洞见/常见误区/史话/本节要点/学习提示）在黑白下退化为同一灰度的 8 个同形框，**语义层级完全丢失**。
- **修复建议**（为每个语义框增加非颜色的**形状/线型/图标**冗余编码）：
```latex
% 1) 线型区分：不要全部 solid
\tcbset{
  % 定理：实线细框（默认）
  % 定义：双线
  % 推论：虚线
  % 例题：点划线
}
\newtcbtheorem[number within=chapter]{corollary}{推论}{
  colback=cororange!5,colframe=cororange,fonttitle=\bfseries,breakable,
  boxrule=0.8pt,dash pattern=on 2pt off 1.5pt,   % 虚线 -> 黑白可辨
  left=4pt,right=4pt,top=4pt,bottom=4pt}{cor}
\newtcbtheorem[number within=chapter]{example}{例题}{
  colback=exampurple!5,colframe=exampurple,fonttitle=\bfseries,breakable,
  boxrule=0.8pt,dash pattern=on 4pt off 1pt,    % 点划线
  left=4pt,right=4pt,top=4pt,bottom=4pt}{ex}

% 2) 标题前加图标（黑白下也是形状差异）
\tcbset{fonttitle=\heiti\small, before title={\raisebox{-1pt}{\scalebox{1.1}{$\blacksquare$}}\ },}
% 或按环境分别设 title 前置符号：■/◆/●/▲/★
```
  同时把 `warnyellow` 加深以过 WCAG AA（正文小字需 ≥4.5:1）：
```latex
\definecolor{warnyellow}{RGB}{200, 120, 0}   % 对比度 2.04:1 -> 约 4.0:1
% 或保持色相但压暗：{216, 122, 0} ≈ 4.3:1；建议直接 {180, 95, 0} ≈ 5.8:1
```
  `cororange` (3.79:1) 同样低于 AA，建议压暗至 `{200,60,0}`（≈5.0:1）。

### P3-2　红绿色盲（deuteranopia）下 5 组语义色互相混淆

- **定位**：同上 `:127-135`
- **问题**：经 deuteranopia 变换后 RGB 距离 <45 的配对：

| 配对 | 模拟后 RGB | 距离 |
|---|---|---|
| thmblue ↔ reviewblue | (81,81,255) / (90,90,212) | 45 |
| defgreen ↔ historybrown | (101,101,51) / (95,95,71) | **22** |
| defgreen ↔ tipgreen | (101,101,51) / (116,116,61) | **23** |
| cororange ↔ green60(框) | (124,124,0) / (108,108,3) | **23** |
| historybrown ↔ tipgreen | (95,95,71) / (116,116,61) | 31 |

- **影响**：约 8% 男性读者的红绿视觉下，「定义（绿）」与「学习提示（绿）」、「物理史话（棕）」三代同色；「推论（橙）」与「检查清单（绿60）」同色。
- **修复建议**：在蓝–橙轴上重排，避免绿/棕/红三者共存于同一亮度带：
```latex
\definecolor{defgreen}{RGB}{0, 121, 107}      % 改青绿（teal），脱离红绿轴
\definecolor{tipgreen}{RGB}{27, 94, 32}       % 与定义拉开明度
\definecolor{historybrown}{RGB}{93, 64, 55}   % 压暗，脱离 defgreen 亮度带
\definecolor{insightred}{RGB}{176, 0, 32}     % 偏品红，与 cororange 分离
% 验证工具：\usepackage{colorblind} 或在 PDF 上跑 deuteranopia 模拟
```
  配合 P3-1 的线型/图标冗余，即使色相仍接近也可辨识。

### P3-3　页面上下留白不对称（正文约占页高 78%）

- **定位**：`src/plasma_physics_textbook_v2.tex:7-8`
- **问题**：`top=2.1cm`、`bottom=2.1cm`、`headsep=0.5cm`、`footskip=1.0cm`。页眉占据 top 边距内的 0.5cm+15pt，实际正文起点被推低；底部 2.1cm 中 1.0cm 是未使用的 footskip ⇒ 视觉上**顶部紧凑、底部空旷**。
- **复现证据**：第 20 页实际测量——正文首行基线约在距页顶 2.95cm 处，末行基线约在距页底 1.35cm 处；89 页页脚带墨迹为 0。
- **修复建议**：
```latex
\geometry{left=2.1cm,right=2.6cm,top=2.2cm,bottom=1.7cm,
          footskip=0.5cm,headsep=0.5cm}
% right 由 2.3 -> 2.6cm 同时缓解 P2-1 的边注越界（20.62cm 落回 20.80cm 盒内边界内）
```

### P3-4　TOC 层级与末级页码拥挤

- **定位**：`src/plasma_physics_textbook_v2.tex:332`（`\tableofcontents`）；PDF 第 2–7 页（6 页目录）
- **问题**：TOC 占 6 页，`\subsection` 级条目（如 `11.1.1 气体电离与汤森放电机制`）与页码之间几乎没有引线缓冲，`\cftsetindents` 未调整，末级缩进偏小。
- **复现证据**：`pdftotext -layout -f 5 -l 5`：
  ```
  11.1.1 气体电离与汤森放电机制 . . . . . . . . . . . . . . . . . . . . . . . .    128
  11.1.2 帕邢定律：击穿电压 . . . . . . . . . . . . . .                 128
  ```
  同行内第二个条目的点引线明显短于第一个（`\@dottedtocline` 在长标题下挤压）。
- **修复建议**（统一点引线并增加缩进）：
```latex
\usepackage{tocloft}   % 或在 ctex 下直接调 \cftsetindents
\cftsetindents{subsection}{3.8em}{3.2em}
\cftsetindents{section}{1.5em}{2.4em}
\renewcommand{\cftdotsep}{4.5}          % 点引线密度
\setlength{\cftbeforesubsecskip}{1pt}
% 若目录仍过长，可只收 tocdepth=2（章 + 节），subsection 不进目录
\setcounter{tocdepth}{2}
```

### P3-5　页眉左标记为西文斜体，与中文正文风格断裂

- **定位**：`src/plasma_physics_textbook_v2.tex:107`
- **问题**：`\fancyhead[L]{\nouppercase{\small\itshape\leftmark}}`。`\leftmark` 内容是中文章名（如「等离子体基础诊断」），`\itshape` 对中文触发合成倾斜（同 P2-6），且 auto-`\MakeUppercase` 已被 `\nouppercase` 抑制但斜体保留。
- **复现证据**：第 20 页页眉 300dpi 实际观察：「Chapter 1. 数学与统计力学预备知识」中中文部分呈机械倾斜，与正文宋体正体不一致。
- **修复建议**：
```latex
\fancyhead[L]{\nouppercase{\small\heiti\leftmark}}   % 中文用黑体，不用斜体
% 或 \small\kaishu 与章前引语呼应
```

### P3-6　hyperref 书签锚点警告（符号说明章）

- **定位**：日志 `Package hyperref Warning: The anchor of a bookmark and its parent's must not be the same. Added a new anchor on input line 347.`；对应 `src/plasma_physics_textbook_v2.tex:337-338`
- **问题**：`\chapter*{符号说明}` + `\addcontentsline{toc}{chapter}{符号说明}` 使书签父子锚点相同，hyperref 自动补锚点。共 2 处同类警告。
- **复现证据**：日志原文（见上）；`pdftotext -f 8 -l 8` 确认第 8 页为「符号说明」章。
- **修复建议**（用 `\phantomsection` 显式建锚）：
```latex
\chapter*{符号说明}
\phantomsection
\addcontentsline{toc}{chapter}{符号说明}
% 前言、结语等所有 \chapter* + \addcontentsline 组合处同样处理
```
  另建议核对 PDF 书签树：`\usepackage{bookmark}` 替代 hyperref 内建 bookmark 可消除此类警告并生成更规范的书签层级。

---

## 5. 已验证为**良好**的部分（避免误改）

以下经实际观察确认无问题，评审建议保持：

1. **封面（`src/plasma_physics_textbook_v2.tex:220-330`）质量高**。第 1 页 200/300dpi 观察：三层径向 fade 消除了硬边，环向辉光 + 三条缠绕放电丝层次分明，四层同心环面线给出正确透视；标题块（46pt 黑体「等离子体物理」+ 21pt 西文副标 + letterspace 260 的 `A SELF-STUDY TEXT` + 分隔线 + 作者块）层级清晰、对齐居中、无溢出。**唯一可议**是封面纯视觉主图无任何文字标注（`tex:313` 注释亦确认这是有意设计），对自学教材可能缺少「这是什么」的提示——但属设计取向，非缺陷。
2. **`tcbset` 统一色块样式（`:192-202`）执行正确**。注释明确指出「必须在 15 个定义之后，否则被各自的 `left=4pt` 反压」，实测所有 15 类框的 `left/right/top/bottom` 一致为 5pt/5pt/3pt/3pt，`sharp corners=southwest` 统一生效，框内文字与框线间距均匀。
3. **`titlesec` 章标题设计有效**。第 19/29/39 页观察：64pt `coverMid!30` 大号章号 + 黑体章名 + `titlerule[1.2pt]` 分隔线，视觉冲击力与层级感良好。
4. **全书无表格水平溢出**（除 P2-7 的单元格换行问题）。像素扫描未发现任何表格线越过正文右边界；booktabs 三线表使用规范。
5. **数学公式排版质量高**。`amsmath`/`mathtools`/`bm`/`esint` 组合下，多行对齐、矩阵、矢量粗体、`siunitx` 单位（`\mathrm{eV}`、`\mathrm{m/s}`）均正确；`\abovedisplayskip`/`\belowdisplayskip` 收紧到 5pt 有效减少了公式周边留白。
6. **`microtype` 已启用**且加载了 lmroman/lmtt 的 protrusion 配置（日志确认 `mt-cmr.cfg`、`mt-msa.cfg`、`mt-msb.cfg` 加载，`lmtt` 使用通用 protrusion 设置）。
7. **PDF 元数据完整**。`pdfinfo` 确认 Title/Subject/Keywords/Author/Creator/Producer/PageSize 全部正确，`pdfcreationdate` 固定为 `D:20260913000000Z` 以保证可复现构建。

---

## 6. 修复优先级建议

| 顺序 | 条目 | 理由 |
|---|---|---|
| 1 | P1-2 + P1-3 | 同一处 `\appendix` 机制引发，一次改动修复 25 页页眉 + 全部附录编号，收益最大 |
| 2 | P1-1 | 9 处缺字导致核心教学图不可读，改动局部（`ch10_supplement.tex` 4 行） |
| 3 | P1-4 | 原始 TeX 泄漏属正确性问题，需全书图内标注排查 |
| 4 | P2-3 + P2-4 | 图无编号 + 难度标记错位，影响全书导航体验，但机械改动量大 |
| 5 | P2-1 + P3-3 | 同属边距问题，`right: 2.3cm → 2.6cm` 一并缓解 |
| 6 | P2-2 | 3 处文案重复，一行修复 |
| 7 | P2-5 | `plain` 页式 5 行修复，感知明显 |
| 8 | P2-6 + P3-5 | 中文合成斜体，需统一 CJK 字体声明 |
| 9 | P2-7 | 单表列宽 |
| 10 | P3-1 + P3-2 | 配色冗余编码，改动面大但可增量推进（先做线型） |
| 11 | P3-4 + P3-6 | 打磨项 |

---

## 附录：审查脚本与原始证据

本轮生成的中间证据（均在 `reviews_tmp/`，非交付物）：

| 文件 | 内容 |
|---|---|
| `full.txt` | `pdftotext -layout` 全文（211 段/分页符） |
| `scan_out.txt` | 210 页每页前 5 行摘要（用于定位章首页与附录页） |
| `header_out.txt` | 每页页眉带/页脚带墨迹像素数 + 空页眉清单（12 页） |
| `bbox_out.txt` | 每页墨迹左右边界 + 越界清单（50 页） |
| `edge_out.txt` | 越界分类：20 页越过边注盒 20.50cm（均为 20.62cm） |
| `rightink_out.txt` | 越界墨迹的 y 区间（统一 15px 高 → 判定为 `\rule` 竖线） |
| `color_out.txt` | 9 语义色 WCAG 对比度、灰度亮度、26 组灰度贴近配对、5 组色盲混淆 |
| `log_warn.txt` | 日志中 Overfull/Underfull/Missing character/警告提取 |
| `hdr/g-*.png` | 210 页 100dpi 灰度栅格 |
| `pg*.png` / `cover_*.png` / `p*_*.png` | 关键页与局部 200–400dpi 裁剪 |

复现命令：
```bash
pdftoppm -r 100 -gray -png pdf/plasma_physics_textbook_v2.pdf reviews_tmp/hdr/g
pdftotext -layout pdf/plasma_physics_textbook_v2.pdf reviews_tmp/full.txt
python reviews_tmp/colors.py
python reviews_tmp/edges.py        # bbox.py + edge.py + rightink.py 的结论来源
```
