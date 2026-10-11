# -*- coding: utf-8 -*-
"""可证伪性夹具：故意在源码里植入本轮的典型回退，确认 check/physics_regression.py 会 FAIL。

流程：备份 → 植入 78 处回退 → 跑套件（期望 退出码≠0 且 FAIL≥MIN_FAIL）→ `finally` 恢复 → 再跑（期望 0）。其中四处打在 v1.8.0/v1.8.1 新增的附录 I/J/K 上（K.5 的阻尼率符号、K.7 的鞘厚量级与单位、K.4 的磁冻结结论、J 的逐章对照表章号），
三处打在 v1.10.0 新增的附录 M 上（束速量级、H 模阈值标度、欧姆/自加热倍数口径），
一处是 v1.10.0 的结构性回退（语义盒退回裸 `\begin{tip}`，正文首字会被吞、标签可能重复），
三处打在 v1.11.0 新增的附录 N 上（脱靶辐射功率量级、ELM 侵蚀倍数、表标签被改名），
三处打在 v1.12.0 新增的附录 O 上（电子磁化程度 ω_ceτ_e 的量级、可调比冲的取舍结论、损失清单表标签被改名），
三处打在 v1.13.0 新增的附录 P 上（§P.2 帕申外推与经验值的比值、§P.4 表内刻蚀速率量级、来源表标签被改名），
两处打在 v1.13.0 新增的附录 I「间隔与提取：把复习排进日历」一节上（遗忘曲线的幂律指数、排程表 D=84 那行的复习日）。
两处打在 v1.14.0 新增的附录 Q 上（探针判据表里的半径/德拜比、I--V 斜率结论里的低估比例、RC 积分输出电压、单发散射光子数、来源表标签被改名，
共五处），
一处打在 v1.14.0 新增的附录 I「十种学习技术的评级」一节上（把「分散练习」从高效用降档，档次计数 2/3/5 被打破），
一处是 v1.14.0 修好的真排版缺陷回退（把 `\textbf{}` 改回 Markdown `**` 简写——v1.13.0 已发布 PDF 里确实印出了字面星号），
七处打在 v1.15.0 新增的附录 R 上（表 R.2 的半径/德拜比、Z_d 量级、§R.4 表内 Γ_eff、Havnes 警戒线、
尘埃声速量级、辐射压判据 β(0.1 μm)、来源表标签被改名），
三处打在 v1.15.0 新增的附录 I「双编码与多媒体学习」一节与「动作六」上（规则表标签被改名、
动作六小节被改名、两本期刊被拼成一个刊名的真错回退），
五处打在 v1.15.1 勘误的附录 R 量级口径上（表 R.2 的 C4 半径/德拜比十倍错、密度量级跨度写法、
尘埃等离子体频率量级、尘埃与离子质量比区间、氚库存的克数与摩尔数——这五条取自 v1.15.0
发出去的口径错（正文共修七处），回退它们用来证明新增的字面绑定与 `lacks` 守卫都能红）。
十二处打在 v1.16.0 上：附录 S 的十一处（IPB98 的 P 幂次、κ 因子、回旋半径的 2T 约定、
回旋玻姆写成 ρ*²、玻姆系数 1/16、表 S.4 的错式 τ_E、玻姆型外推 τ_E、回旋玻姆型标度、
Greenwald 极限、测地声模频率、四条对照表的一行）——其中 2T 约定与 ρ*² 两处同时触发
字面绑定与 `lacks` 守卫，各产生 2 个 FAIL；另一处打在第 10 章公式卡上（P 幂次退回 −0.58），
它会同时打掉「公式卡含 P^{-0.69}」「公式卡不得含 P^{-0.58}」与「旧错式只许出现在附录 S」
三条断言，是跨文件检查 P20 的可证伪性证据。
十九处打在 v1.17.0 上：附录 T 的十五处（q 剖面幂次、剪切恒等式、平行电子热导的电子质量、
电子碰撞时间的 3/2 次幂、岛宽系数 4、反转半径的根号、Ψ_h'' 与剪切的对应、表 T.2 的 L_q、
表 T.3 的 χ_∥、表 T.4 的岛宽、表 T.6 的 j_bs、§T.7 的体积份额、一条「出处：」被改名、
一道练习的解答被删、表 T.4 标签被改名）、主文件一处（附录 T 的 \appchapter 被改成 U，
用来证明附录字母连续性检查会红）与附录 I 新节三处（认知层级表标签被改名、
康奈尔笔记的复习动作被翻转、参考文献表里的 Bloom 条目被改写）。其中剪切恒等式、岛宽系数、
反转半径、Ψ_h''、电子质量与 3/2 次幂这六处同时被字面绑定与 `lacks` 守卫覆盖，
是 P22 双层防线的可证伪性证据。
一处打在 v1.17.1 勘误上（§T.4 正文印出的 τ_flat 从 3.128e-15 退回 3.125e-15——
v1.17.0 发布时表与脚本已对齐、正文漏改，两者之差落在任何数值用例的容差之内，
只能靠逐字绑定抓；这条补丁同时触发字面绑定与 `lacks`）。
用来证明 P 类新检查确实会失败——不是「从来没红过」的装饰：K.5 那条只翻符号，绝对值一字不变。
恢复用 `newline=''` 写回，保证不把 LF 改写成 CRLF（否则整个文件变成「已修改」）。

对应验收：reviews/2026-10-10-independent/TODO.md 完成条件「复算本轮负例与正对照……
使『修复通过』不依赖旧错误仍存在」。

用法：python check/fixtures/gate_physics_negative.py
"""
import io
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAIN = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.tex')
CH9 = os.path.join(ROOT, 'src', 'ch9_supplement.tex')
CH6 = os.path.join(ROOT, 'src', 'ch6_supplement.tex')
APPJ = os.path.join(ROOT, 'src', 'appendix_j_textbook_map.tex')
APPI = os.path.join(ROOT, 'src', 'appendix_i_learning.tex')
APPK = os.path.join(ROOT, 'src', 'appendix_k_problems.tex')
APPM = os.path.join(ROOT, 'src', 'appendix_m_heating.tex')
APPN = os.path.join(ROOT, 'src', 'appendix_n_boundary.tex')
APPO = os.path.join(ROOT, 'src', 'appendix_o_propulsion.tex')
APPP = os.path.join(ROOT, 'src', 'appendix_p_lowtemp.tex')
APPQ = os.path.join(ROOT, 'src', 'appendix_q_diagnostics.tex')
APPR = os.path.join(ROOT, 'src', 'appendix_r_dusty.tex')
APPS = os.path.join(ROOT, 'src', 'appendix_s_confinement.tex')
APPT = os.path.join(ROOT, 'src', 'appendix_t_topology.tex')
CH10 = os.path.join(ROOT, 'src', 'ch10_supplement.tex')
SUITE = os.path.join(ROOT, 'check', 'physics_regression.py')

PATCHES = [
    (CH9, r'56.6\,\mathrm{ps}', r'5.66\,\mathrm{ps}'),                     # 例题数值回退
    (MAIN, '0.0543', '0.0600'),                                            # Buneman 标度回退
    (CH6, r'\omega=\frac{\omega_{pe}}{\cos\theta_0}', r'\omega=\omega_{pe}\cos\theta_0'),
    (MAIN, r'w^2\approx-0.1101', r'w^2\approx-0.2202'),                    # 双流复核例回退
    # 附录 K.5：只翻符号、量级不变——旧的「按绝对值查存在性」检查抓不到，只有带符号匹配能抓
    (APPK, r'=-2.00\times10^{-2}$', r'=+2.00\times10^{-2}$'),
    (APPK, r'=1.81\times10^{-3}\ \mathrm{m}', r'=1.81\times10^{-2}\ \mathrm{m}'),  # 鞘厚量级回退
    (APPK, r'\approx1.8\ \mathrm{mm}', r'\approx1.8\ \mathrm{cm}'),         # 鞘厚单位回退
    # 附录 K.4：把「v_A/c 小」重新当成磁冻结的充分条件（旧错误结论回归）
    (APPK, '不能据此断定磁冻结成立', '磁冻结图像在非相对论等离子体中非常好'),
    (APPJ, '第 10 章 & 约束与受控核聚变', '第 13 章 & 约束与受控核聚变'),   # 附录 J 幻章/漏行
    # 附录 M（v1.10.0）：三个不同层面的回退——数值量级、经验标度、结论口径
    (APPM, r'\approx3.10\times10^{6}\ \mathrm{m/s}', r'\approx3.10\times10^{5}\ \mathrm{m/s}'),
    (APPM, r'P_{\mathrm{LH}}\approx86.5\ \mathrm{MW}', r'P_{\mathrm{LH}}\approx8.65\ \mathrm{MW}'),
    (APPM, '差 22 倍', '差 220 倍'),
    # v1.10.0 结构性回退：两参数语义盒退回裸调用（会吞掉正文首字、并可能造成标签重复）
    (MAIN, r'\begin{tip}{磁应力记忆口诀}{tip-magstress-mnemonic}', r'\begin{tip}'),
    # 附录 N（v1.11.0）：三处不同层面——数值量级、结论里的倍数、表标签被改名（结构性）
    (APPN, r'\approx0.861', r'\approx0.761'),                              # 脱靶辐射份额量级回退
    (APPN, r'ELM 造成的侵蚀比稳态高约 $390$ 倍', r'ELM 造成的侵蚀比稳态高约 $39$ 倍'),
    (APPN, r'\label{tab:app-n-impurity}', r'\label{tab:app-n-impuritx}'),
    # 附录 O（v1.12.0）：三处不同层面——磁化程度的量级、可调比冲的取舍结论、损失清单表标签被改名
    (APPO, r'\omega_{ce}\tau_e\approx663', r'\omega_{ce}\tau_e\approx66.3'),
    (APPO, '比冲调高一倍，推力就减半', '比冲调高一倍，推力就翻倍'),
    (APPO, r'\label{tab:app-o-loss}', r'\label{tab:app-o-lost}'),
    # 附录 P（v1.13.0）：三处不同层面——帕申外推与经验的比值（口径数值）、
    # 刻蚀速率量级、来源表标签被改名（结构性）
    (APPP, '经验值的 $1.68$ 倍', '经验值的 $1.86$ 倍'),          # §P.2 帕申外推与经验的比值
    (APPP, '& 1.0 & 323.4 &', '& 1.0 & 32.34 &'),               # §P.4 表内刻蚀速率量级
    (APPP, r'\label{tab:app-p-source}', r'\label{tab:app-p-sourcx}'),
    # 附录 I（v1.13.0）：两处不同层面——遗忘曲线的幂律指数（结论层的借用数据）、
    # 复习排程表里 D=84 那一行的复习日（规则输出的数据层）
    (APPI, r'$b\approx0.13$', r'$b\approx0.31$'),
    (APPI, r'$1,\,3,\,7,\,15,\,27,\,39,\,51,\,63,\,75$',
     r'$1,\,3,\,7,\,15,\,27,\,39,\,51,\,63,\,73$'),
    # 附录 Q（v1.14.0）：五处不同层面——探针判据表内的半径/德拜比（数据层）、
    # I--V 斜率结论里的低估比例（结论层）、RC 积分时间常数（数值层）、
    # 汤姆孙散射单发光子数（量级层）、来源表标签被改名（结构性）
    (APPQ, '40.72 & 12.28', '40.72 & 1.228'),                    # §Q.2 探针判据表
    (APPQ, r'低估 $50.63\%$', r'低估 $5.063\%$'),                # §Q.3 不减离子电流的代价
    (APPQ, r'50.0\ \mathrm{\mu s}', r'5.00\ \mathrm{\mu s}'),    # §Q.6 RC 积分时间常数
    (APPQ, '1781.6', '178.16'),                                  # §Q.8 单发散射光子数
    (APPQ, r'\label{tab:app-q-source}', r'\label{tab:app-q-sourcx}'),
    # 附录 I（v1.14.0）：Dunlosky 评级表里把「分散练习」从高效用降档（档次计数 2/3/5 被打破）
    (APPI, '高效用 & 分散练习', '低效用 & 分散练习'),
    # v1.14.0 修好的真排版缺陷回退：v1.13.0 已发布 PDF 里 `**` 被原样印成星号，
    # 这一处把 \textbf 改回 Markdown 简写，用来证明 check_P16_no_markdown_bold 真的会红
    (APPP, r'\textbf{不能把低压帕申公式外推到大气压}',
     r'**不能把低压帕申公式外推到大气压**'),
    # 附录 R（v1.15.0）：七处不同层面——表 R.2 的 a/λ_D（数据层）、Z_d 量级、
    # §R.4 表内 Γ_eff（结论层，同一数值在 §R.10 练习 3 里仍正确，用来验证
    # 「数值必须落在自己那一节」这条断言）、Havnes 警戒线、尘埃声速量级、
    # 辐射压判据 β(0.1 μm)、来源表标签被改名（结构性）
    (APPR, r'$40.72$ & $0.0123$', r'$40.72$ & $0.123$'),
    (APPR, r'$2.46\times10^{3}$', r'$2.46\times10^{4}$'),
    (APPR, r'4.35\times10^{3}', r'4.35\times10^{2}'),
    (APPR, r'4.07\times10^{13}', r'4.07\times10^{12}'),
    (APPR, r'$c_{DA}=1.5455\ \mathrm{cm/s}$', r'$c_{DA}=15.455\ \mathrm{cm/s}$'),
    (APPR, r'$\beta=2.87>1$', r'$\beta=0.287>1$'),
    (APPR, r'\label{tab:app-r-source}', r'\label{tab:app-r-sourcx}'),
    # 附录 I（v1.15.0）：三处——双编码规则表标签被改名（结构性）、
    # 「动作六」小节被改名（动作清单退回五个）、
    # 以及 v1.14.0 查出的真错回退（把两本期刊拼成一个刊名，DOI 反证）
    (APPI, r'\label{tab:app-i-dualcode}', r'\label{tab:app-i-dualcodx}'),
    (APPI, r'\subsection{动作六：把图讲成句子}', r'\subsection{动作五之二：把图讲成句子}'),
    (APPI, r'\textit{Psychological Science in the Public Interest} 14(1), 4--58；DOI: 10.1177/1529100612453266',
     r'\textit{Perspectives on Psychological Science in the Public Interest} 14(1), 4--58；DOI: 10.1177/1529100612453266'),
    # 附录 R（v1.15.1）：五处——把 v1.15.0 发出去的四类量级口径错逐一回退，
    # 证明新增的数值断言、字面绑定与 lacks 守卫都能抓住（其中第 46 处就是
    # v1.15.0 已发布 PDF 里表 R.2 的十倍错：C4 的 a/λ_D 应为 1.3e-8）
    (APPR, r'$1.3\times10^{-8}$', r'$1.3\times10^{-7}$'),
    (APPR, r'电子密度跨十三个数量级、', r'电子密度跨十四个数量级、'),
    (APPR, r'$\omega_{pd}$ 也降到 $4.1\times10^{3}\ \mathrm{rad/s}$',
     r'$\omega_{pd}$ 也降到 $4.1\times10^{2}\ \mathrm{rad/s}$'),
    (APPR, r'$10^{5}$--$10^{14}$ 倍', r'$10^{5}$--$10^{13}$ 倍'),
    (APPR, r'$0.75\ \mathrm{g}$，', r'$0.075\ \mathrm{g}$，'),
    # v1.16.0：附录 S（能量约束与输运）十处 + 第 10 章公式卡一处
    (APPS, r'P^{-0.69}', r'P^{-0.68}'),                                    # IPB98 的 P 幂次
    (APPS, r'\kappa^{0.78}', r'\kappa^{0.77}'),                            # 漏掉 κ 因子的等价写法
    (APPS, r'\rho_i=\frac{\sqrt{2m_iT_i}}{eB}', r'\rho_i=\frac{\sqrt{m_iT_i}}{eB}'),
    (APPS, r'\chi_{\mathrm{gB}}=\chi_B\,\rho_\ast', r'\chi_{\mathrm{gB}}=\chi_B\,\rho_\ast^2'),
    (APPS, r'\chi_B=\frac{T_e[\mathrm{eV}]}{16\,B[\mathrm{T}]}',
     r'\chi_B=\frac{T_e[\mathrm{eV}]}{8\,B[\mathrm{T}]}'),                 # 玻姆系数口径
    (APPS, '$7.4446$', '$7.4464$'),                                        # 错式 τ_E（表 S.4）
    (APPS, r'\tau_E^{\mathrm{Bohm}}=0.426', r'\tau_E^{\mathrm{Bohm}}=0.462'),
    (APPS, r'\tau_E\propto\frac{a^3B^2}{T^{1.5}}', r'\tau_E\propto\frac{a^3B}{T^{1.5}}'),
    (APPS, r'1.064\times10^{20}', r'1.64\times10^{20}'),                   # Greenwald 极限
    (APPS, r'f_{\mathrm{GAM}}\simeq\frac{c_s}{2\pi R}',
     r'f_{\mathrm{GAM}}\simeq\frac{c_s}{4\pi R}'),                         # 测地声模频率
    (CH10, r'P^{-0.69}', r'P^{-0.58}'),                                    # 公式卡退回旧错式
    # v1.16.0：附录 I「专家与新手」一节
    (APPI, '表征习惯', '表示习惯'),                                         # 四条对照之一被改
    # v1.17.0：附录 T（磁面、磁岛与磁通量管）
    # ① 公式形式层：q 剖面的幂次、剪切恒等式、平行电子热导的质量、电子碰撞时间的 3/2 次幂、
    #    岛宽系数、反转半径的根号、Ψ_h'' 与剪切的对应——七处都由 P22 的字面绑定 / lacks 守卫
    (APPT, r'q(\rho)=q_0+(q_{95}-q_0)\,\rho^2', r'q(\rho)=q_0+(q_{95}-q_0)\,\rho'),
    (APPT, r's(r)=2\,\frac{q(r)-q_0}{q(r)}=2-\frac{2q_0}{q(r)}', r's(r)=2-\frac{q_0}{q(r)}'),
    (APPT, r'\chi_\parallel=\frac{3.16\,T_e\,\tau_e}{m_e}',
     r'\chi_\parallel=\frac{3.16\,T_e\,\tau_e}{m_i}'),
    (APPT, r'\tau_e=\frac{3.44\times10^{5}\,T_e^{3/2}}{n_e\ln\Lambda}',
     r'\tau_e=\frac{3.44\times10^{5}\,T_e^{1/2}}{n_e\ln\Lambda}'),
    (APPT, r'w=4\sqrt{\frac{\tilde\Psi}{B_\theta\,s}}',
     r'w=2\sqrt{\frac{\tilde\Psi}{B_\theta\,s}}'),
    (APPT, r'r_{\mathrm{inv}}=a\sqrt{\frac{1-q_0}{q_{95}-q_0}}',
     r'r_{\mathrm{inv}}=a\frac{1-q_0}{q_{95}-q_0}'),
    (APPT, "\\Psi_h''(r_{mn})=B_\\theta\\,s", "\\Psi_h''(r_{mn})=B_\\theta\\,s^2"),
    # ② 数据层：表 T.2 的 L_q、表 T.3 的 χ_∥、表 T.4 的岛宽、表 T.6 的 j_bs、§T.7 的体积份额
    (APPT, r'1.6013', r'1.6103'),
    (APPT, r'9.023\times10^{11}', r'9.023\times10^{10}'),
    (APPT, r'0.0698', r'0.0689'),
    (APPT, r'1.755\times10^{5}', r'1.775\times10^{5}'),
    (APPT, r'4.8\%', r'48\%'),
    # ②b 正文印出的 τ_flat 与表/脚本不一致（v1.17.1 修掉的旧写法，只能逐字钉住）
    (APPT, r'=3.128\times10^{-15}', r'=3.125\times10^{-15}'),
    # ③ 结构层：一条「出处：」被改名（计数跌破 7）、一道练习的解答被删、表标签被改名
    (APPT, r'\textbf{出处：}Spitzer', r'\textbf{来源：}Spitzer'),
    (APPT, r'\textbf{解答}：$\rho=\sqrt{', r'解答：$\rho=\sqrt{'),
    (APPT, r'\label{tab:app-t-w}', r'\label{tab:app-t-wx}'),
    # ④ 主文件：附录字母连续性（\appchapter 的 T 被改成 U）
    (MAIN, r'\appchapter{附录T：磁面、磁岛与磁通量管：拓扑如何决定约束（注明出处）}',
     r'\appchapter{附录U：磁面、磁岛与磁通量管：拓扑如何决定约束（注明出处）}'),
    # v1.17.0：附录 I「学习目标与笔记结构」——表标签、康奈尔笔记的复习动作、参考文献表条目
    (APPI, r'\label{tab:app-i-bloom}', r'\label{tab:app-i-bloomx}'),
    (APPI, r'遮住右栏，只看左栏', r'遮住左栏，只看右栏'),
    (APPI, r'Bloom 等 (1956) & \textit{Taxonomy of Educational Objectives',
     r'Bloom 等 (1956) & \textit{Taxonomies of Educational Objectives'),
]

MIN_FAIL = 74        # 实测值（78 处植入回退合并计入 82，留 8 的余量）：见下
# 逐条实测（tmp/probe_neg_v1100.py … tmp/probe_neg_v1151.py、tmp/probe_neg_v1160.py）
# 每处补丁单独植入都至少产生 1 个 FAIL，v1.16.0 单独计共 69 个；合并植入时降为 56，因为 P01（附录 K 数值）、
# P08（附录 M 数值）、P10（附录 N 数值）、P11（附录 O 数值）与 P12/P13/P14/P17/P19 各自把同一附录/同一小节内的多处不符汇总成一行。
# 故 56 是「合并植入」的实测值，MIN_FAIL 取 50 留余量，per-patch 可证伪性由上面的探针单独保证
# （v1.16.0 实测：58/58 全部命中，无落空，未植入时 exit=0 FAIL=0）。
# v1.15.0 的探针还抓出一个检查侧的洞：§R.6 的 1.5455 cm/s 改成 15.455 cm/s 后，
# 同节 1.55 Hz 落在 0.5% 相对容差窗内，数值断言仍通过 → P17 因此补了 9 条小节内字面绑定。
# v1.15.1 的探针又抓出第二个洞：P17 只对 C1 复算了 a/λ_D，于是表 R.2 里 C4 的十倍错
# （1.3e-7 应为 1.3e-8）在 419 项全绿的情况下发出去 → 四个工况的 a/λ_D 与尘埃/离子
# 等离子体频率全部逐行复算，量级口径另加 5 条字面绑定与 4 条 lacks 守卫。
# v1.16.0 的探针一次抓出四个洞（附录 S 侧）：① `P^{-0.69}` 与 `\kappa^{0.78}` 在正文、
# 表说明与练习解答里各出现 2–3 次，只改其中一处的补丁仍能让松散的字面绑定通过 → 字面绑定
# 改成钉住整条式子 `n_{19}^{0.41}P^{-0.69}` 与 `\epsilon^{0.58}\kappa^{0.78}M^{0.19}`；
# ② 表 S.4 的 4 位小数值用 0.1% 容差，把 7.4446 抄成 7.4464 的转写错落在容差内 → 六个数
# 全部收紧到 1e-4；③ 夹具锚点 `$0.426` 先命中正文句而不是公式 → 锚点改为 `\tau_E^{\mathrm{Bohm}}=0.426`
# 并补同款字面绑定；④ 锚点 `$1.064\times10^{20}$` 在原稿里根本不存在（真串带 `\ \mathrm{m^{-3}}`）。
# 另记一条流程教训：批量探针偶发漏报（第 54 号一度报 FAIL=0），凡 `[X]` 必须用
# `tmp/probe_one_v1160.py <序号>` 单条复核后才可判定为洞。
# v1.17.0（附录 T + 附录 I 学习目标节）植入扩到 77 处，逐条探针（`tmp/probe_neg_v1170.py`）
# 抓出两个新洞：第 67 号（表 T.3 的 `9.023\times10^{11}`，全篇出现 2 次）与第 68 号
# （表 T.4 的 `0.0698`，全篇出现 4 次）——数值用例只要求「这个数在文中某处出现」，
# 补丁只改第一处（表里那一处）时用例仍由正文/练习里的同值满足。修法不是在夹具里改成全局替换，
# 而是在 P22 里把**每一处印出的位置**都钉住：两个整行字面绑定（表 T.3 的 D1 行、表 T.4 的
# D1 3/2 行）＋四处整句绑定（`\chi_\parallel=9.02\times10^{11}`、`\chi_\parallel=9.023\times10^{11}`、
# `w_{2/1}=0.0531\ \mathrm{m}<w_{3/2}=0.0698\ \mathrm{m}`、`w=0.0698\times\sqrt4=0.1396\ \mathrm{m}`、
# 练习 4 的给定值）。修后单条复核两处均为 `FAIL=1`，整轮 77/77 命中、单独计 94、合并计 81。
# v1.17.1 勘误轮植入扩到 78 处：新增第 78 号把 §T.4 正文里的 τ_flat 从 3.128e-15 退回
# 3.125e-15（v1.17.0 发布时表与脚本已对齐、正文漏改；两者之差落在任何数值用例容差内，
# 只有逐字绑定与 `lacks` 能抓）。单条复核 `FAIL=1`，合并植入实测 **82**，MIN_FAIL 取 74。
# 同一轮的勘误还包含附录 T 七张表的列宽（数学在窄 `p` 列里会断成两行，观感问题不是门禁项，
# 由渲染 PNG 目视发现——与 v1.15.1 的附录 R 表 R.2 同类）。


def run_suite():
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'          # 控制台是 GBK，直接跑会 UnicodeEncodeError
    p = subprocess.run([sys.executable, SUITE], capture_output=True, env=env)
    t = p.stdout.decode('utf-8', 'replace')
    nf = 0
    for line in t.split('\n'):
        if line.startswith('汇总：'):
            nf = int(line.split('FAIL')[1].split('／')[0].strip())
    if nf == 0 and '汇总' not in t:
        err = p.stderr.decode('utf-8', 'replace').strip().split('\n')[-1:]
        print('    （套件未输出汇总，疑似异常：%s）' % (err or ['?']))
    return p.returncode, nf, t


def write(path, text):
    """newline='' 保证不把 LF 改写成 CRLF（否则整个文件都变「已修改」）。"""
    io.open(path, 'w', encoding='utf-8', newline='').write(text)


def main():
    backup, original = {}, {}
    for path, _, _ in PATCHES:
        if path not in backup:
            original[path] = io.open(path, 'rb').read()
            backup[path] = original[path].decode('utf-8')
    ok = True
    try:
        for path, old, new in PATCHES:
            t = io.open(path, encoding='utf-8').read()
            if old not in t:
                print('[X] 夹具锚点未找到：%s -> %s' % (os.path.basename(path), old))
                ok = False
                continue
            io.open(path, 'w', encoding='utf-8', newline='').write(t.replace(old, new, 1))
        code, nf, _ = run_suite()
        print('[%s] 植入回退后：退出码=%d FAIL=%d（期望 !=0 且 >=%d）'
              % ('OK' if (code != 0 and nf >= MIN_FAIL) else 'X', code, nf, MIN_FAIL))
        ok = ok and code != 0 and nf >= MIN_FAIL
    finally:
        for path, t in backup.items():
            write(path, t)
    code, nf, _ = run_suite()
    print('[%s] 恢复后：退出码=%d FAIL=%d（期望 0 且 0）'
          % ('OK' if (code == 0 and nf == 0) else 'X', code, nf))
    ok = ok and code == 0 and nf == 0
    # 按**原字节**比对恢复结果：既不允许多出 CRLF，也不允许任何行尾/字节漂移。
    # （不能只查「有没有 CRLF」——core.autocrlf=true 的检出本身可能就是 CRLF。）
    for path in backup:
        b = io.open(path, 'rb').read()
        if b != original[path]:
            d = abs(len(b) - len(original[path]))
            print('[X] %s 未按原字节恢复（CRLF 变化 %+d 处，长度差 %d 字节）'
                  % (os.path.basename(path), b.count(b'\r\n') - original[path].count(b'\r\n'), d))
            ok = False
    print('夹具总判定：%s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
