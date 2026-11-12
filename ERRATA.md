# ERRATA.md —— 勘误总台账

本表记录 2026-11 统一编辑轮对本书源码的逐条核对与修改。八列固定为：
**ID｜位置｜原文/现象｜问题｜改法｜证据｜状态｜备注**。

证据分级：**A** = 已读源码并独立复算/量纲核对；**B** = 有第三方公开来源可查；**C** = 需作者拍板的编辑决定，不算物理错误。
状态取值：**已修**、**待作者**、**待核**、**不改**、**已复核无误**。
行号以本轮改后的源码为准；改前的行号写在「位置」的括号里。

| ID | 位置 | 原文/现象 | 问题 | 改法 | 证据 | 状态 | 备注 |
|---|---|---|---|---|---|---|---|
| PL-01 | `src/plasma_physics_textbook_v2.tex:2912`（原 2912）、`:3064`、`:6827` | insight 与章末回顾写 `ω²=ω_pe²+3k²v_te²`；附录 G 写 `ω²=ω_pe²+(3k_BT_e/m_e)k²=ω_pe²+3k²v_te²` | 同节定理 `:2898` 用单分量宽度 `σ_ve=√(k_BT_e/m_e)`，且全书 `v_te=√2σ_ve`；因此 `3k²v_te²=6k²σ_ve²`，热修正项大 2 倍 | 三处统一为 `3k²σ_ve²`；附录 G 补一句 `σ_ve` 的定义与 `v_te=√2σ_ve` | A：代数核对；`grep -rn "3k\^2v_{te}" src/*.tex` 无匹配；数值 `8.32e17` 未变 | 已修 | 对应 qa-verification B10；与 PL-02 同源 |
| PL-02 | `src/ch6_supplement.tex:41,46,50,58,77,412,423`；`src/appendix_l_problems.tex:145-146` | 补篇把 `v_te` 定义成 `√(k_BT_e/m_e)=8.39e5 m/s`，并以此写 `ω²=ω_pe²+3k²v_te²`；附录 L 题与答案也用 `3v_te²` | 该定义实际是 `σ_ve`，与全书 `v_te=√(2k_BT_e/m_e)` 冲突；按全书口径会把热项放大 2 倍 | 补篇 6 处变量改名 `σ_ve`，公式卡与符号说明同步；附录 L 题与答案改 `3k²σ_ve²`、`3σ_ve²` | A：数值 `8.32e17`、`2.62e-4`、`1.49e-5 m` 全部保持不变；`λ_De=σ_ve/ω_pe` 自洽 | 已修 | 局部约定已在 `ch7_supplement.tex:132`、`ch9_supplement.tex:207` 明说的两处未动 |
| PL-03 | `src/plasma_physics_textbook_v2.tex:5699`（原 5664） | 第 12 章量纲自检：`φ=\frac{e^2}{2ε_0m_ec^2}\int n_eB\,dl` | 多一个 `B`、少一个 `ω`；量纲不是 rad，与定理 `:5234` 的 `Δφ=e²/(2ε_0m_ecω)∫n_e dl` 不一致 | 改为 `φ=e²/(2ε_0m_ecω)∫n_e dl`；量纲句补 `r_eλ_0` 同量纲说明 | A：`r_eλ_0` 为 `m²`，乘 `n_e` 与 `dl` 后无量纲；与 `eq:phase-shift` 逐项一致 | 已修 | 对应 qa-verification B9 |
| PL-04 | `src/plasma_physics_textbook_v2.tex:3515`（目标）、`:3621`（新节） | 第 9 章目标承诺「瑞利--泰勒、开尔文--亥姆霍兹和双流不稳定性的推导和图像」，正文只有 RT 与双流；KH 仅在术语表 `:7228` | 目标与正文不符，读者按目录找不到 KH | 在双流之后新增 `§9.4 开尔文--亥姆霍兹不稳定性`：不可压剪切流色散、`v_0>2v_A` 磁稳定判据、磁层顶/边界层图像、SpacePhysics 第 6 章对应；更新学习路径与章末回顾 | A：色散关系与判据独立推导，三条极限与量纲核对通过；full_check 仍为 12 条量纲自检、18 只自测盒；consistency 14/14 | 已修 | 新增 6 个 label、2 处 `\ref`；未新增图、未新增盒、未改回归期望 |
| PL-05 | `src/appendix_p_lowtemp.tex:14,103,118,129,138,154,155,478,480,542` | 正文「帕申」与「帕邢」混用；改前全局帕邢 28 次、帕申 12 次 | 同一术语两种译名 | 正文统一为「帕邢定律」「帕邢曲线」；公式、数值、图注物理内容不动 | A：改后全局帕邢 39 次、帕申 1 次（仅标题，见 PL-06）；C：译名取通行译法 | 已修 | 对应 qa-verification B8；属编辑统一，不计物理错误 |
| PL-06 | `src/appendix_p_lowtemp.tex:77` | `\section{汤森雪崩与帕申定律}` | 与正文统一目标冲突；但 `check/physics_regression.py:1887,1985` 以该标题前缀作为固定期望键，改名会使 7195 项回归出现 FAIL | 暂保留标题；建议作者批准后改标题，并同步更新回归脚本的期望键 | A：`grep` 定位；C：改标题需作者拍板 | 待作者 | 若改标题，需把 `want_prefix` 两处一起改，再重跑 `physics_regression.py`；本轮不动期望值 |
| PL-07 | `src/ch2_supplement.tex:55,232,241`；主册 `§2.3` 定理 `plasma-criteria` | 例题 (iii) 写「准中性」；速查卡写 `ω_peτ_coll>1`；清单写「准中性」 | 三处措辞不同，容易被读成三个独立判据；物理上弱耦合 `Γ≪1` 与 `N_D≫1` 同源 | 以 `Γ≪1` 为主判据；速查卡注明「等价于常写的 `ω_peτ_coll>1`」；清单补 `Γ=(1/3)N_D^{-2/3}` 把三者连起来 | A：与主册关系式一致；consistency 14/14；物理判据未改 | 已修 | 对应 qa-verification B5；性质是「表述统一」，不是改判据 |
| PL-08 | `src/appendix_aa_alternates.tex:63,103,111,128,323`；`src/appendix_ad_pedestal.tex:162,197,201`；`src/appendix_ae_kinetics.tex:295,314,320,594,621` | AA 用 D1 的 9 MA、`q_95=3.00`；AD/AE 用 15 MA、`q_95=3.5` | 读者可能当成同一工况下的矛盾 | 不改数值；记为不同工作点；建议 AA `:103` 附近加一句「`q_95=3.00` 对应 9 MA 对比工况，AD/AE 的 15 MA 设计工况取 3.5」 | A：两处电流与 q95 均在文中声明；q95 是输入参数而非推导结果 | 不改 | 对应 qa-verification B7；若 AD 改用 3.0，`α_MHD≈2.41` 会跨 2.5，属敏感性问题，留作者判断 |
| PL-09 | `src/ch1_supplement.tex`（2 处）、`src/ch2_supplement.tex`（4 处）、`src/ch4_supplement.tex`（1 处）、`src/ch5_supplement.tex`（6 处） | 行首出现孤立的「。\item」 | 排版噪声，列表项前多一个句号 | 删除行首「。」，保留 `\item` | A：改后 `grep -rn "^\。" src/*.tex` 无匹配 | 已修 | 内容审查顺手修，不改语义 |
| PL-10 | `src/plasma_physics_textbook_v2.tex:6827` 附近 | 「端口输出」 | 应为「脚本输出」；与上下文（Python 脚本打印）不符 | 改为「脚本输出」 | A：上下文核对 | 已修 | 同一段已因 PL-01 改过 |
| PL-11 | `README.md`；`pdf/plasma_physics_textbook_v2.pdf` | 「610 页」 | 需核对是否与实际 PDF 一致 | 无（保留 610） | A：`pdfinfo` 四册为 610 / 207 / 403 / 47；分册和 610 自洽 | 复核无误 | 对应 qa-verification B12 |
| PL-12 | 主册 `§2.3`；`src/appendix_ac_numerics.tex:56`、`src/appendix_ae_kinetics.tex:37` | `N_D` 与 `Γ` 的定义及换算 | qa-verification B6 怀疑 `N_D` 口径不一致 | 无；三处都明确 `Γ=(1/3)N_D^{-2/3}`，且 `σ_ve`/`v_te` 关系写清 | A：按定义复算，`Γ=1/3 N_D^{-2/3}` 自洽 | 复核无误 | 不改正文 |
| PL-13 | `src/plasma_physics_textbook_v2.tex:5234-5245` | 定理 `eq:phase-shift`：`Δφ=(ω/c)∫(1-n_r)dl=e²/(2ε_0m_ecω)∫n_e dl` | 作为 PL-03 的正确锚点复核 | 无 | A：与均匀近似 `r_eλ_0∫n_e dl` 一致 | 复核无误 | 不改正文 |
| PL-14 | `src/ch2_supplement.tex:101`、`:127` | 两处 `\ref{ex:low-temp-params}`，全库无此 label | 编译会印 `??`；目标应为 `ex:typical-params`（主文件 `:1805` 例题，`newtcbtheorem` 自动加 `ex:` 前缀） | 两处改为 `\ref{ex:typical-params}` | A：`grep` 无 `ex:low-temp-params`；目标 label 存在；consistency 仍 14/14 | 已修 | 审计 P0-1；旧 PDF 显示「例题 2.3」是更早 label 版本，当前源原本会印 `??` |
| PL-15 | `src/plasma_physics_textbook_v2.tex:3069`、`:3264` | 两处量纲自检写 `k^2v_{Te}^2` | 大写 `v_{Te}` 与全书 `v_te=√2σ_ve` 口径不一致；同页其他公式已用 `σ_ve` | 两处改为 `k^2\sigma_{ve}^2` | A：`grep -n "v_{Te}"` 无匹配；两处均与色散关系 `k²σ_ve²` 同量纲 | 已修 | 审计 PL-008；不改数值 |
| PL-16 | `src/plasma_physics_textbook_v2.tex:3350` | 「$12\pi n_e\lambda_D^3$ ……与本书约定只差一个常数因子 $4\pi$」 | 算术错：$12\pi n_e\lambda_D^3$ 是 $9N_D$，与本书的 $3N_D$ 相差 3 倍，不是 $4\pi$ | 改为「与本书的 $3N_D$ 相差 3 倍；若以 $N_D$ 计，则相差 9 倍」 | A：代入 $N_D=\frac{4\pi}{3}n_e\lambda_D^3$ 复算 | 已修 | 审计 PL-017；不改 $\ln\Lambda$ 的定义与数值 |

## 计数

- 已修：11（PL-01、PL-02、PL-03、PL-04、PL-05、PL-07、PL-09、PL-10、PL-14、PL-15、PL-16）
- 待作者：1（PL-06）
- 待核：0
- 不改：1（PL-08）
- 已复核无误：3（PL-11、PL-12、PL-13）
- 合计：16 条

**计数口径**：按表体 `| PL-` 行计数，PL-01 … PL-16 共 16 行；已修 11 + 待作者 1 + 待核 0 + 不改 1 + 已复核无误 3 = 16。本文件只按表体计数；若他处仍引用 13/10 等旧数，以本段与表体为准。

**门禁口径**：文中 `7195` 指完整构建（有 `.aux`）的物理回归项数；`consistency 14/14` 指 14 个 PASS 项全过。本机静态（无 TeX）为 `7193` 项、一致性 15 项（14 PASS + 1 INFO），差额与原因见 [`STATUS.md`](STATUS.md) 与 [`README.md`](README.md) 的双栏表。

## 本轮没有改的事

- 没有升版本号；主文件封面、`pdftitle` 与 README 仍是 v1.33.9。
- 没有改 `\input` 路径、label 名、构建脚本、`check/physics_regression.py` 的期望值。
- 没有构建 PDF：本机无 XeLaTeX/latexmk，四册 PDF 仍是 **2026-09-13** 构建的 v1.33.9 旧产物，与当前源码不同步。
- 没有编造引用；KH 节只引用本书 `\ref` 与《SpacePhysics》第 6 章 §6.4 的定性对应。
