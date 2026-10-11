# CLOSURE —— 计划评审 6 项问题的修复与验证记录（v1.8.1）

- 评审来源：`reviews/2026-10-11-scheduled/`（`findings.json` 基线 `b95117b3ee8e8e53b4c3d2fbfe1e5ef3a368ce8a`，即 v1.8.0 文档提交）。
- 待办清单：`reviews/2026-10-11-scheduled/TODO.md`（6 项，均为 P2；末行要求「修订完成后写 CLOSURE，逐项记录新版源码、实际 PDF 页码、反例行为与验收结果」）。
- 本文件即对该要求的回答。历史评审原文（`01-content.md`、`02-ui-reading.md`、`03-workflow.md`、`04-previous-round.md`、`checks/`、`evidence/`、`recheck_review.py`）**一律不回写**，保留原样以备对照。

## 0. 口径

1. **页码一律指新版 PDF 的物理页码**（`pdfinfo` 的 Pages 计数，从封面第 1 页算起），由 PyMuPDF 在最终交付的 PDF 上重新搜索定位（`tmp/locate_v181.py`、`tmp/locate_v181b.py`），不沿用评审报告里的旧页码——正文在两轮之间变长了（附录 I/J/K 新增），旧页码只能作对照。
2. **「修复」的判据是可证伪的机械证据**：源码断言（`check/physics_regression.py`）、夹具植入回退后必须 FAIL（`check/fixtures/gate_physics_negative.py`）、最终 PDF 上的实测定位。三者缺一不算关闭。
3. **不为了压 WARN 而删内容**：可读性 WARN、书签 −1 偏移等如实记录在 §7，不做粉饰。
4. 版本号：本轮为 **v1.8.1**（内容勘误 + 门禁加固），不改设计、只做增补与勘误。

## 1. 交付物、环境与干净构建

| 项目 | 值 |
| --- | --- |
| 主文件 | `src/plasma_physics_textbook_v2.tex`（6253 行） |
| 新增附录源 | `src/appendix_i_learning.tex`、`src/appendix_j_textbook_map.tex`、`src/appendix_k_problems.tex`（v1.8.0 引入，本轮勘误） |
| 最终 PDF | `pdf/plasma_physics_textbook_v2.pdf` = `src/plasma_physics_textbook_v2.pdf`（逐字节一致，`Copy-Item -Force` 同步） |
| 页数 | **241 页** |
| 字节数 | **1 816 654 B** |
| SHA256 | `7011AA9A9FDDA2192CA6A90695274D44A37654DB058049664A544C3DDDE9F08F` |
| PDF 标题 | `等离子体物理自学教材 v1.8.1`（`pdftitle`，主文件 `:85`；封面版本串在主文件 `:373`） |
| 编译 | `latexmk -xelatex -interaction=nonstopmode`，首跑 exit 0、次跑 `All targets are up-to-date` exit 0 |
| TeX 错误 | `^! ` 计数 **0** |
| 缺字 | 仅白名单豁免项（`U+000A`） |
| Overfull | 仅 2 处历史基线：主文件 `2048--2049` 0.84378 pt、`ch7_supplement.tex:155` 6.76945 pt（本轮新增的三处溢出已全部消除，见 §6） |
| 环境 | XeTeX 3.141592653-2.6-0.999998（TeX Live 2026）、pdfinfo 25.02.0、Python 3.14.7、numpy 2.5.3、scipy 1.18.1、PyMuPDF 1.28.2、matplotlib 3.11.1 |

## 2. 六套检查与五个可证伪夹具（均在最终 PDF 上复跑）

| 套件 | 判定 | 说明 |
| --- | --- | --- |
| `check/full_check.py` | **FAIL 0 / WARN 1 / 32 项** | WARN 1 即上表两处历史基线 Overfull；v1.8.0 曾报「README 最新版本未出现在源码封面」的 FAIL，抬版本号后消除 |
| `check/ui_readability.py` | **FAIL 0 / WARN 17 / 35 项** | 可读性类告警（36 句 >100 字、句长 p90 = 65、2 条超长段 369/351 字），未为压数删文献条目 |
| `check/consistency_audit.py` | **FAIL 0 / 14 项** | 标签唯一、章节连贯、附录字母连续、README 计数与正文一致 |
| `check/physics_regression.py` | **FAIL 0 / 156 项**（v1.8.0 为 122 项） | 本轮重写 P 类：逐题绑定 + 有符号 + 单位绑定，并新增 P05/P06 |
| `check/figure_layout.py` | **FAIL 0 / 13 项** | 图 9.1/9.2/9.3 的节点框与文字排版 |
| `check/pdf_outline.py` | **FAIL 0 / WARN 2 / 11 项** | O1 页数 241、书签 220（1 级 26 / 2 级 137 / 3 级 57）、链接 323、印刷目录 7 页（第 2--8 页，动态探测）；两条 WARN 见 §7 |

| 夹具 | 判定 |
| --- | --- |
| `check/fixtures/gate_compile_log.py` | **PASS**（7 例；其已知前提：若 README 页数与正文不符会连带 FAIL） |
| `check/fixtures/gate_physics_negative.py` | **PASS**（植入 **9** 处回退 → exit 1 / FAIL **11**；恢复后 exit 0 / FAIL 0，且文件按原字节恢复） |
| `check/fixtures/gate_figure_layout.py` | **PASS**（5 例合成 PDF，含 `textcut`） |
| `check/fixtures/gate_pdf_outline.py` | **PASS** |
| `check/fixtures/gate_g2_listing.py` | **PASS**（含篡改数值、删行两个反例） |

报告文件：`tmp/final_<suite>_v181.txt`、`tmp/final_gate_*_v181.txt`、汇总 `tmp/final_all_v181.txt`。

## 3. 逐项关闭表

| 编号 | 问题 | 源码位置（新版） | 实际 PDF 页 | 反例行为 | 判定 |
| --- | --- | --- | --- | --- | --- |
| P-C01 | 曲率漂移方向与曲率矢量混淆 | `src/appendix_k_problems.tex:79-87` | **p.237** | `check_P05` 缺 `\unitvec{b}\times\kappa` 或出现旧句即 FAIL | ✅ |
| P-C02 | CL 鞘厚与 Bohm 通量密度不一致 | `src/appendix_k_problems.tex:181-205` | **p.239--240** | 厚度或通量单侧变异即 FAIL（夹具 2 条） | ✅ |
| P-C03 | 用非相对论性推出磁冻结 | `src/appendix_k_problems.tex:108-114`、`:122` | **p.238** | 结论反转（夹具第 8 条）即 FAIL | ✅ |
| P-W01 | 门禁忽略答案符号且不绑定题目 | `check/physics_regression.py:736-1063` | —（工具） | 只翻符号 → FAIL；量级/单位/幻章变异 → FAIL（合计 11 项） | ✅ |
| PP-C04 | EEPF 与归一化 EEDF 仍混淆 | `src/plasma_physics_textbook_v2.tex:4069-4101`、`:4108` | **p.164--165** | 三种记法缺一、丢 √ε、缺 J/eV 雅可比 → FAIL | ✅ |
| PP-C06 | 散射判据符号表反向 | 主文件 `:6047`（符号表）、`:6200`（术语表） | **p.221**（符号表）、**p.223**（术语表） | 表格行方向反转 → FAIL（P06） | ✅ |

## 4. 逐项明细

### P-C01 曲率漂移方向 —— 改为 $\unitvec{b}\times\kappa$ ✅

- **评审所指**：K.3（`$\grad B$ 漂移的方向判别与电荷分离`）把「曲率漂移的方向」写成指向曲率中心，与正文柱坐标叉乘结论矛盾。
- **新版源码**（`src/appendix_k_problems.tex:79-87`）：
  - `曲率矢量 $\kappa=(\unitvec{b}\cdot\grad)\unitvec{b}$ 确实指向曲率中心（环向场下 $\kappa=-\unitvec{R}/R_c$），但曲率漂移沿 $\unitvec{b}\times\kappa$，\emph{垂直于}曲率矢量`；
  - 显示式 `\vect{v}_R=\frac{mv_\parallel^2}{qB}\,(\unitvec{b}\times\kappa)` 与 `\unitvec{b}\times\kappa=\unitvec{\varphi}\times(-\unitvec{R}/R_c)=+\unitvec{z}/R_c`；
  - 结论句 `漂移方向是竖直的（正离子向上、电子向下），\textbf{不是}指向曲率中心`。
- **常见错误**同步补一条：`把曲率漂移的方向当作曲率矢量的方向（漂移 $\propto\unitvec{b}\times\kappa$，与 $\kappa$ 垂直，两者不可互换）`（`:91`）。
- **门禁**：`check_P05_appendix_k_claims()` 断言 K.3 必须含 `\unitvec{b}\times\kappa`、`\textbf{不是}指向曲率中心`，且**不得**再出现旧句 `曲率漂移沿磁力线的曲率中心方向`。
- **实测页码**：p.237（搜「不是指向曲率中心」；同页亦有「曲率中心」上下文）。

### P-C02 Child--Langmuir 鞘厚与 Bohm 通量 —— 统一到鞘边密度并互相回代 ✅

- **评审所指**：K.7 用体密度德拜长度估厚度，却用鞘边密度算通量，两者不守同一个电流边界条件。
- **新版源码**（`src/appendix_k_problems.tex:181-205`）：
  - 题干区分三问：`（b）求\emph{体}等离子体的德拜长度`、`（c）用 $s\approx\frac{\sqrt2}{3}\lambda_D(n_s)(2V_w/T_e)^{3/4}$ 估算 Child--Langmuir 鞘层厚度（$\lambda_D(n_s)$ 为鞘边密度处的德拜长度）`、`（d）……说明它与 (c) 的厚度是否守同一个电流边界条件`；
  - 数值：`λ_D(n_0)=1.29e-4 m=0.129 mm` → `λ_D(n_s)=λ_D(n_0)/√0.61=1.65e-4 m` → `s=0.471×1.65e-4×23.3=1.81×10^{-3} m≈1.8 mm`；
  - 明确写出误用体密度的后果：`若误用体密度算出的 $\lambda_D(n_0)$，厚度会被低估 $\sqrt{n_0/n_s}=1.28$ 倍（得到 1.4 mm），并使厚度与 (d) 的通量不守同一个电流边界条件`；
  - (d) 回代：`Γ=n_s c_s=6.1×10^{15}×2.68×10^{3}=1.6×10^{19} m^{-2}s^{-1}`、`j=eΓ≈2.62 A m^{-2}`；再用 `j=(4/9)ε_0√(2e/M)V_w^{3/2}/s^2` 反解得 `1.813 mm`，`与估算式相差 0.2\%`；结论句 `厚度与通量守同一个电流边界条件`。
- **门禁**：`check_P01` 逐题绑定 K.7 的 6 个数值（含 `1.6478e-4 m`、`1.8099e-3 m`、`1.8 mm`、`Γ`、`j`、`1.8132e-3 m`，容差 0.5%--3% 且**绑定单位**）；`check_P05` 还要求字面出现 `\lambda_D(n_s)`、`Child--Langmuir`、`同一个电流边界条件`。
- **夹具**：`=1.81\times10^{-3}\ \mathrm{m}` → `=1.81\times10^{-2}\ \mathrm{m}`（量级）、`\approx1.8\ \mathrm{mm}` → `\approx1.8\ \mathrm{cm}`（单位）两处变异都必须 FAIL——后者只能靠 P05 的字面断言抓到（同块内还有 `1.813 mm`，纯数值检查会漏）。
- **实测页码**：p.239--240（「鞘边密度」p.239 起，「1.813」在 p.240）。

### P-C03 非相对论性 ≠ 磁冻结 —— 拆成独立条件 ✅

- **评审所指**：K.4 用 $v_A/c=2.6\%\ll1$ 直接推出磁冻结成立。
- **新版源码**（`src/appendix_k_problems.tex:108-114`、`:122`）：
  - `这里只能得出一个结论：\textbf{非相对论 Alfvén 速度公式适用}（位移电流可以忽略），\textbf{不能据此断定磁冻结成立}`；
  - 给出真正的条件：`磁冻结要求理想欧姆定律 $\vect{E}+\vect{v}\times\vect{B}=0$ 成立，至少要磁雷诺数 $R_m=\mu_0Lv/\eta\gg1$，还要检查霍尔项与电子压强项是否可丢（第 5 章「MHD 适用条件：三道独立的门槛」）`；
  - 强调输入不足不能下结论：`本题没有给温度、电阻率与流动尺度，用现有输入无法判断——\emph{同一个} $v_A/c$ 既可属于强碰撞的电阻性等离子体，也可属于理想 MHD 等离子体`；
  - 给出定量对照：`R_m=1.2566×10^{-6}×1×7.69×10^{6}/1.64×10^{-5}≈5.9×10^{5}\gg1` 与扩散时标 `τ_η=μ_0L^2/η≈0.076 s`；
  - 常见错误补一条：`把 $v_A/c\ll1$ 当成磁冻结的充分条件（它只保证非相对论公式适用，$R_m=\mu_0Lv/\eta\gg1$ 才是丢掉电阻项的条件）`。
- **门禁**：`check_P05` 要求 K.4 含 `不能据此断定磁冻结成立` 与 `R_m=\mu_0Lv/\eta`，且不得出现反转句 `磁冻结图像在非相对论等离子体中非常好`（夹具第 8 条即植入该反转）。
- **实测页码**：p.238（搜「磁冻结成立」）。

### P-W01 门禁忽略符号且不绑定题目 —— 重写 P 类检查 ✅（本条本身就是「检查是否可信」）

- **评审所指**：`check/physics_regression.py:672-703` 的 `values_in()`/`has_value()` 在全篇范围内按**绝对值**查数字是否出现过，既不看符号，也不要求答案出现在所属题目里。
- **新版实现**（`check/physics_regression.py`，`tmp/p01_new.py` 经 `tmp/apply_p01_patch.py` 拼接替换第 675--770 行，96 行 → 325 行；仍只依赖标准库）：
  - `NUM_PAT`（`:736`）：支持 `a\times10^{b}`、`a e b`、纯小数三种写法；
  - `_is_negative(text, pos)`（`:753`）：向前取 12 字符并右去空白，末字符为 `-` 且其前为 `=(,:$[{` 之一或以 `SIGN_HEADS`（`\approx`/`\sim`/`\simeq`/`\pm`）结尾才算**负号**——因此 `10^{-2}` 的指数负号不会误判；
  - `value_in(text, v, rel, unit='', signed=False)`（`:775`）：在**指定块**内找值；`unit` 必须出现在数值后 80 字符内且其后不是字母（`\mathrm{m}` 不会错配 `\mathrm{mm}`）；
  - `split_blocks(text, kind)`（`:811`）：`kind='example'` 切出附录 K 的 8 道例题（每块只含该 `example` 环境），`kind='section'` 切出附录 I 的小节；
  - `check_P01_new_appendix_numbers()`（`:835`）：**43** 条 case，元组带 `(tag, 题目序号, 值, 相对容差, 单位, 是否带符号, 说明)`；先自检 8 道题与 ≥4 个 I 小节都在，再逐条在**所属块**内核对；K.5 的 `γ/ω` 用 `signed=True`（负号必须在）。
- **新增两条方向性门禁**：
  - `check_P05_appendix_k_claims()`（`:955`）：三道判断题的结论方向（K.3 / K.4 / K.7）字面核对；
  - `check_P06_scattering_alpha()`（`:990`）：对全篇「含 `&` 且含 α 判据」的表格行按子句切分判定方向（`α≪1` 子句必须是非相干/非集体、`α≳1`/`α≫1` 子句必须是集体/相干），实测命中 2 行且方向全部正确；另加 4 条显式断言（符号表行、术语表行、12.4 正文、`ch12_supplement.tex` 自测清单）；正文中的同句并存只记 INFO。
- **反例行为**（夹具 `check/fixtures/gate_physics_negative.py`，`PATCHES` 9 条、`MIN_FAIL = 11`）：

  | # | 植入位置 | 变异 | 目的 |
  | --- | --- | --- | --- |
  | 1 | `ch9_supplement.tex` | `56.6 ps` → `5.66 ps` | 数值量级 |
  | 2 | 主文件 | `0.0543` → `0.0600` | 数值容差 |
  | 3 | `ch6_supplement.tex` | `ω_pe/\cosθ_0` → `ω_pe\cosθ_0` | 公式方向 |
  | 4 | 主文件 | `w^2≈-0.1101` → `-0.2202` | 主文件与补充材料交叉一致 |
  | 5 | 附录 K | `=-2.00\times10^{-2}` → `=+2.00\times10^{-2}` | **只翻符号**（绝对值不变） |
  | 6 | 附录 K | `1.81\times10^{-3} m` → `1.81\times10^{-2} m` | 量级 |
  | 7 | 附录 K | `≈1.8 mm` → `≈1.8 cm` | 单位绑定 |
  | 8 | 附录 K | 结论反转（磁冻结） | 结论方向 |
  | 9 | 附录 J | `第 10 章` → `第 13 章` | 幻章（README/正文无此章） |

  - 实测：植入后 `exit=1`、`FAIL=11`；恢复后 `exit=0`、`FAIL=0`，并且恢复是按 `io.open(path,'rb').read()` 的**原字节**比对确认的（仓库 `core.autocrlf=true`，不能用文本比较）。
  - 迭代记录：首轮只得 FAIL 10，用临时探针 `tmp/probe_neg.py` 定位到第 4 条旧补丁**此前从未命中任何检查**（旧 `check_P01` 只查附录 I/K 的 26 个值），补上主文件交叉对照后才达到 11——即该变异此前是**能悄悄通过门禁的**，这正是 P-W01 的攻击面。
- **结果**：`python check\physics_regression.py` → **FAIL 0 / WARN 0 / 156 项**；module docstring 已同步说明 P 类的四项加固。

### PP-C04 三种 EEDF 记法 —— 列出归一化、单位与 J/eV 雅可比 ✅

- **评审所指**：v1.8.0 的修订仍把 EEPF 与「归一化 EEDF」混为一谈。
- **新版源码**（主文件 `:4069-4108`）：
  - `\textbf{先分清三种互相冲突的记法}——它们相差\emph{两次不同的运算}：密度归一化、以及除去 $\sqrt{\eps}$。这两步不能合并，也不能靠「改单位」或乘一个常数互换`；
  - 三行表：`EEDF $F(\eps)$`（`∫F dε=n_e`，`m^{-3}(eV)^{-1}`）、`归一化 EEDF $p=F/n_e$`（`∫p dε=1`，`(eV)^{-1}`）、`EEPF $g=p/\sqrt{\eps}$`（`∫√ε g dε=1`，`(eV)^{-3/2}`）；
  - `\textbf{EEPF 不是「归一化 EEDF」的同义词}：它比归一化 EEDF 又多除了一个 $\sqrt{\eps}$` + 文献口径（Lieberman & Lichtenberg、Cantera PlasmaPhase）+ 麦克斯韦试金石 `p∝√ε e^{-ε/T_e}`、`g∝e^{-ε/T_e}`、`ε=0.25,1,4 eV` 处 `p/g=√ε=0.5,1,2`；
  - 「每焦耳」换算：`F_{\mathrm{eV}}(\eps_{\mathrm{eV}})=e\,F_J(e\,\eps_{\mathrm{eV}})`、`∫F_eV dε_eV=∫F_J dε_J`，并警告 `只改单位栏而不做这一步，整条曲线的数值会差 $e$ 倍`；
  - 雅可比段落：`F(\eps)=4πv^2f(v)(dv/dε)`、`d²I_e/dV²=e² d²I_e/dε²`；
  - Druyvesteyn 定理（`:4100`）补 $p$ 与 $g$ 的两种换算（`F=n_e√ε g`）；警告框（`:4108`）补 `\textbf{只有把 $\sqrt{\eps}$ 除掉的 EEPF $g=p/\sqrt{\eps}$ 才呈纯指数型}——这正是低气压放电文献偏爱 $g$ 的原因`；insight（`:4112`）的被积函数由 `f_\eps(\eps)` 改为 `F(\eps)`（全书已无 `f_\eps`/`\hat F` 残留）。
- **门禁**（`check_M07_eedf()` 扩展）：(d) 数值验证 `p/g=√ε`（0.25/1/4 eV）与 `∫√ε g dε=1`；(e) 验证 J↔eV 雅可比下 `∫F_eV dε_eV=∫F_J dε_J`（相对差 <1e-9）；源码断言三种记法表三行原文、`EEPF 不是「归一化 EEDF」的同义词`、`F=n_e\sqrt{\eps}\,g`、`F_{\mathrm{eV}}=e F_J`、探针通量积分用 `F(\eps)`，并 `lacks('EEPF 通常就是')`、`lacks('f_\\eps(')`。
- **实测页码**：p.164--165（`EEPF` 两页都有；「每焦耳」与「纯指数型」在 p.165）。

### PP-C06 附录 H 符号表的 α 判据 —— 同步为非集体 ✅

- **评审所指**：正文已改对，符号表仍反向（评审基线行 `:6028`）。
- **新版源码**：主文件 `:6047` `$\alpha$（散射参数） & 散射参数 & 无量纲 & \S 12.4，$\alpha = 1/(k\lambda_D)$；$\alpha\ll1$ 为\textbf{非集体（非相干）散射}，$\alpha\gtrsim1$ 才需集体散射理论`；术语表 `:6200` `集体散射 & collective scattering & $\alpha \gtrsim 1$ 时的散射 regime（探测尺度大于德拜长度）`。（行号由 6028 变为 6047，是 PP-C04 在同一轮新增 19 行所致。）
- **全文对照**：`check_P06_scattering_alpha()` 扫描全篇表格行，实测**恰好 2 行**含 α 判据且方向全部正确；正文另加 4 条显式断言（含 `ch12_supplement.tex:620`）；唯一 INFO 是 2 处正文把两种情形写在同一句里，属正常表述。
- **实测页码**：附录 H 符号表 p.221（搜「非集体」）、术语表 p.223（搜「collective scattering」）。

## 5. 目视检查记录

- 本轮为**文字与门禁勘误**，未改图、未改版式，故未做新的图像目视；上一轮（v1.8.0）对附录 I/J/K 的表格、长 URL、DOI 溢出的目视结论仍然有效。
- 已复核项目：新 PDF 的 241 页总数、`pdfinfo` 标题、六套件与五夹具的最终报告（`tmp/final_*_v181.txt`）。
- 需要目视时建议优先看：p.237（K.3 叉乘式）、p.239--240（K.7 三段口径）、p.164--165（三记法表）、p.221（符号表行）。

## 6. 评审范围外的相邻修复（本轮顺带）

1. `src/appendix_k_problems.tex:87` 打字错误：`上下交替，,短接了电荷分离` → `上下交替，短接了电荷分离`。
2. K.6 结论句（`:169`）收紧：明确 `\textbf{在给定尺度与流速下}`、`（注意 $R_m$ 还含 $L$ 与 $v$——扩散时标本身只由 $\eta$ 与尺度决定）`，避免把 $R_m$ 的温度依赖说成唯一因素。
3. 上一轮遗留的三处 Overfull 全部消除（主文件附录 E 长 URL 加 `{\small\url{...}}`、附录 I 出处表列宽调整、附录 J 资源条目分行），本轮复核只剩 2 处历史基线。

## 7. 未完成 / 待验证（如实声明）

1. `pdf_outline` 的 **O3b**：18/220 条书签的落点比印刷目录页码少 1（如「等离子体特征参数」印刷 10 / 书签 9）。成因是 `\appchapter`/`\frontchapter`/`\addcontentsline` 的锚点结构，偏移**一致且有界**（7.5%），插 `\needspace` 会在交付前重排页面，故不改，只记录。
2. `pdf_outline` 的 **O4**：附录 F 有 1/35 条链接的目标页与 `.aux` 记录差 1（`sec:ch7-landau`，86 vs 87），与上一条同源。
3. `ui_readability` 的 17 条 WARN：36 句 >100 字、2 条超长段（369/351 字），主要来自附录 I 的英文文献条目与附录 K 的「出处」行。**不为压 WARN 删出处**（注明来源是本轮的硬要求）。
4. `figure_layout` 目前只覆盖第 9 章的 3 张图；其余图（含附录 I/J/K 的表格）尚无版式守卫。
5. 2 处历史基线 Overfull（0.84 pt / 6.77 pt）未处理。
6. `reviews/2026-10-11-scheduled/recheck_review.py` 是评审方脚本，其断言针对评审基线（v1.8.0）的两个对照控制项，本轮按「历史评审不回写」原则未改；验收由六套件 + 五夹具承担。

## 8. 版本与发布

- 版本：**v1.8.1**（主文件 `pdftitle` 与封面版本串已同步；`:6238` 注释里的「v1.8.0 新增」保留，因为它在描述历史）。
- 内容提交：`49797dbfefaad66683f7cf5b16f420d46cae4494`（6 个文件：`README.md`、`check/physics_regression.py`、`check/fixtures/gate_physics_negative.py`、`pdf/plasma_physics_textbook_v2.pdf`、`src/appendix_k_problems.tex`、`src/plasma_physics_textbook_v2.tex`；+466 / −112）。
- 标注标签：**v1.8.1** = `239f2c65f3ae90d9ef3ad493284775dca4f63456` → 内容提交 `49797db`。
- 推送：`refs/heads/main` 与 `refs/tags/v1.8.1` 均已推送到 `origin`（`28ecdaa..49797db`，新标签）。
- Gitea 发布：release **id=86**，`.../releases/tag/v1.8.1`，正文即本轮 release notes；资产 `plasma_physics_textbook_v2.pdf` 上传 HTTP **201**、1 816 654 B。
- 独立复核（下载而非本地文件）：`curl.exe -L` 取回 HTTP **200**、**1 816 654 B**、SHA256 `7011AA9A…9F08F` 与本地一致、文件头 `%PDF-1.7`。至此仓库共有 **14** 个 release。
- 文档提交（本文件 + `processed.md` §13）紧随其后，内容见 Git 历史。

