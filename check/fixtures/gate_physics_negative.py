# -*- coding: utf-8 -*-
"""可证伪性夹具：故意在源码里植入本轮的典型回退，确认 check/physics_regression.py 会 FAIL。

流程：备份 → 植入 317 处回退 → 跑套件（期望 退出码≠0 且 FAIL≥MIN_FAIL）→ `finally` 恢复 → 再跑（期望 0）。其中四处打在 v1.8.0/v1.8.1 新增的附录 I/J/K 上（K.5 的阻尼率符号、K.7 的鞘厚量级与单位、K.4 的磁冻结结论、J 的逐章对照表章号），
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
三十处打在 v1.26.0 的收尾轮上（附录 I 第 17 节「读一篇论文」，门禁 P33）：结构层 11
（节标题、节 label、两张表的 label、语义盒的标题与标签、五个小节标题各改一处）、
口径层 13（表 1 第一遍的「五个 C」退回三个 C、第三遍的「能复现其中一个数」改成两个数、
第二遍行的产出与判据两格整行改字、表 2 的脚本接口从 M--AA 退回 M--X、读者/同行方向翻转、
第一问被换成「它引用了什么」、「读懂」定义句的「一个数」改成「一些数」、笔记模板的
「代入」改成「参考」、四条出处的年份/卷期页/版次/书名各改一处、三条常见错误之一改字、
「一次真正的收获」改成「一次意外的收获」）、`lacks` 回流注入 6（感觉懂了、复现不上就算了、
引用次数越多越正确、第二遍读完就算读懂，以及把本节核心口径「可以随时停止」反写成
「必须读完三遍」）。其中 227/229/231/245 四条同时触发字面绑定与 `lacks`，各产生 2 个 FAIL。
这一轮的逐条探针先抓出 P33 的四个洞——语义盒只钉了标签没钉标题、表 1 漏了第二遍行、
表 2 漏了第一问、以及「三遍法可以随时停止」这句核心口径完全没有绑定（反写成「必须读完
三遍」时六套件全绿）；按 v1.17.0 与 v1.22.0 的做法，修法是回去给 P33 补 4 条字面绑定
与 1 条 `lacks` 守卫（28 -> 33 条断言，门禁 1347 -> 1352 项），而不是改夹具。
补完后 30/30 全部可证伪、无落空。
另记一条夹具自身的坑：语义盒**标题**补丁（220）必须排在**标签**补丁（221）之前——
后者的锚点是前者锚点的子串，而夹具按列表顺序 replace(old, new, 1)，顺序错了会报
「夹具锚点未找到」而不是产生 FAIL。
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
APPU = os.path.join(ROOT, 'src', 'appendix_u_gain.tex')
APPV = os.path.join(ROOT, 'src', 'appendix_v_disruption.tex')
APPW = os.path.join(ROOT, 'src', 'appendix_w_fastparticle.tex')
APPX = os.path.join(ROOT, 'src', 'appendix_x_radiation.tex')
APPY = os.path.join(ROOT, 'src', 'appendix_y_strongfield.tex')
APPZ = os.path.join(ROOT, 'src', 'appendix_z_reactor.tex')
APPAA = os.path.join(ROOT, 'src', 'appendix_aa_alternates.tex')
APPAB = os.path.join(ROOT, 'src', 'appendix_ab_space_astro.tex')
APPAC = os.path.join(ROOT, 'src', 'appendix_ac_numerics.tex')
APPAD = os.path.join(ROOT, 'src', 'appendix_ad_pedestal.tex')
APPAE = os.path.join(ROOT, 'src', 'appendix_ae_kinetics.tex')
APPAF = os.path.join(ROOT, 'src', 'appendix_af_turbulence.tex')
APPAG = os.path.join(ROOT, 'src', 'appendix_ag_equilibrium.tex')
APPAH = os.path.join(ROOT, 'src', 'appendix_ah_atomic.tex')
APPAI = os.path.join(ROOT, 'src', 'appendix_ai_icf.tex')
CH12 = os.path.join(ROOT, 'src', 'ch12_supplement.tex')
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
    # v1.18.0：附录 U（聚变增益 Q、自加热与点火）
    # ① 公式形式层：alpha 份额、f_alpha 的分母、回路线性指数、Q 的单调性、工程增益的括号、
    #    盈亏平衡的 −1、三重积阈值——七处各自同时触发字面绑定与 lacks 守卫（第 7 处只触发字面绑定）
    (APPU, r'P_\alpha = \frac{1}{5}P_{\mathrm{fus}}, \qquad',
     r'P_\alpha = \frac{1}{4}P_{\mathrm{fus}}, \qquad'),
    (APPU, r'= \frac{Q}{Q+5} = \frac{1}{1+5/Q}', r'= \frac{Q}{Q+4} = \frac{1}{1+4/Q}'),
    (APPU, r'P_{\mathrm{tot}}^{\,0.62} < 1', r'P_{\mathrm{tot}}^{\,1.62} < 1'),
    (APPU, r'Q \propto P_{\mathrm{aux}}^{\,-0.38}', r'Q \propto P_{\mathrm{aux}}^{\,0.38}'),
    (APPU, r'= \eta_{\mathrm{th}}\eta_{\mathrm{aux}}\left(Q+1\right)',
     r'= \eta_{\mathrm{th}}\eta_{\mathrm{aux}}\,Q+1'),
    (APPU, r'\eta_{\mathrm{aux}}} - 1 = 5.250', r'\eta_{\mathrm{aux}}} = 6.250'),
    (APPU, r'nT\tau_E > 3.0992\times10^{21}', r'nT\tau_E > 3.0992\times10^{20}'),
    # ② 数据层：整行绑定——三个工作点表、f_alpha 基准行、闭环收敛步、Q(P_aux) 扫描行、工程账
    (APPU, r'D1 ITER 量级    & 150.0 & 2.4762 & 9.194 & 371.4 & 556.5 & 3.710 \\',
     r'D1 ITER 量级    & 150.0 & 2.4762 & 9.194 & 371.4 & 556.5 & 3.170 \\'),
    (APPU, r'3.710  & 0.4259 & 附录 S 的 D1 外部加热点（本书基准） \\',
     r'3.710  & 0.4529 & 附录 S 的 D1 外部加热点（本书基准） \\'),
    (APPU, r'12 & 332.20 & 1.4306 & 11.764 & 911.00 & 182.20 & 6.073 \\',
     r'12 & 332.20 & 1.4306 & 11.764 & 911.00 & 182.20 & 6.037 \\'),
    (APPU, r'5.0   & 25.8825 & 3.203 & 129.4 & 67.5  & 13.510 \\',
     r'5.0   & 25.8825 & 3.203 & 129.4 & 67.5  & 13.150 \\'),
    (APPU, r'D1 外部加热点    & 556.5 & 150.0 & 706.5  & 282.6 & 375.0 & 0.754 \\',
     r'D1 外部加热点    & 556.5 & 150.0 & 706.5  & 282.6 & 375.0 & 0.745 \\'),
    # ②b 正文与练习印出的 nT tau_E 不一致（0.0004e21 落在任何数值用例的容差内，
    #    与 v1.17.1 的 tau_flat 同类：只有逐字钉住整句才能抓）
    (APPU, r'$nT\tau_E = 2.299\times10^{21}\ \mathrm{keV\,s\,m^{-3}}$，是阈值的 $74.2\%$——',
     r'$nT\tau_E = 2.2993\times10^{21}\ \mathrm{keV\,s\,m^{-3}}$，是阈值的 $74.2\%$——'),
    # ③ 结构层：一条「出处：」被改名（计数跌破 5）、一道练习的解答被删、表标签被改名
    (APPU, r'\textbf{出处：} 功率平衡与', r'\textbf{来源：} 功率平衡与'),
    (APPU, r'\textbf{解答}：直接代入式', r'解答：直接代入式'),
    (APPU, r'\label{tab:app-u-qeng}', r'\label{tab:app-u-qengx}'),

    # ── v1.19.0（附录 V：破裂、误差场与逃逸电子）16 处 ──
    # 95. 式 V.5 的 RWM 增长率退回 4 位有效数字的 18.95
    (APPV, r'\gamma_{\mathrm{RWM}} \approx \frac{1}{\tau_w} = 18.947\ \mathrm{s^{-1}}',
          r'\gamma_{\mathrm{RWM}} \approx \frac{1}{\tau_w} = 18.95\ \mathrm{s^{-1}}'),
    # 96. V.3 正文的壁时间退回 52.78 ms
    (APPV, r'给 $\tau_w=52.7788\ \mathrm{ms}$，铜（',
          r'给 $\tau_w=52.78\ \mathrm{ms}$，铜（'),
    # 97. V.5 正文的环向场/Dreicer 场比值退回 0.0151
    (APPV, r'$E_{\mathrm{tor}}/E_D=0.01507$，主体电子仍被碰撞拖住',
          r'$E_{\mathrm{tor}}/E_D=0.0151$，主体电子仍被碰撞拖住'),
    # 98. V.5 的雪崩放大倍数退回粗略的 e^92≈1.0e40
    (APPV, r'$e^{92.12}\approx1.021\times10^{40}$',
          r'$e^{92}\approx1.0\times10^{40}$'),
    # 99. 表 V.2 热猝灭行的功率错一个量级（371.4 GW→37.14 GW）
    (APPV, r'$371.4\ \mathrm{MJ}$ 在 $1\ \mathrm{ms}$ 内是 $371.4\ \mathrm{GW}$ \\',
          r'$371.4\ \mathrm{MJ}$ 在 $1\ \mathrm{ms}$ 内是 $37.14\ \mathrm{GW}$ \\'),
    # 100. 表 V.3 m=3 行的壁面允许场误差退回 3 位有效数字
    (APPV, r'$1.7778\times10^{-4}$',
          r'$1.78\times10^{-4}$'),
    # 101. 表 V.4 工作点行的 beta 被转写错（2.6359→2.6539）
    (APPV, r'本书工作点 1.8627 & 2.6359 &',
          r'本书工作点 1.8627 & 2.6539 &'),
    # 102. 表 V.5 力账行的晕电流力错一个量级（79.5 MN→7.95 MN）
    (APPV, r'力账 & $I_{\mathrm{halo}}=4.5\ \mathrm{MA}$；力约 $79.5\ \mathrm{MN}$',
          r'力账 & $I_{\mathrm{halo}}=4.5\ \mathrm{MA}$；力约 $7.95\ \mathrm{MN}$'),
    # 103. 表 V.6 临界场行的比值小数点错位（769.9→76.99）
    (APPV, r'$E_c$（临界场） & 0.08668 & $E_{\mathrm{tor}}/E_c=769.9$',
          r'$E_c$（临界场） & 0.08668 & $E_{\mathrm{tor}}/E_c=76.99$'),
    # 104. 表 V.7 末行的 E_tor/E_c 被转写错（5.921→5.291）
    (APPV, r'1.00 & 5.133 & 5.921 & 还要慢六倍 \\',
          r'1.00 & 5.133 & 5.291 & 还要慢六倍 \\'),
    # 105. 表 V.8 杂质份额被转写错（1.488%→1.848%）
    (APPV, r'$1.488\%$（$L_z=3\times10^{-32}$） \\',
          r'$1.848\%$（$L_z=3\times10^{-32}$） \\'),
    # 106. 表 V.9 来源行的 Troyon 年份被改错（1984→1985）
    (APPV, r'V.3 & Troyon 等（1984）',
          r'V.3 & Troyon 等（1985）'),
    # 107. 表 V.1 密度极限行的 Greenwald 密度错一个量级
    (APPV, r'附录 M、附录 U（$n_G=1.1937\times10^{20}$） \\',
          r'附录 M、附录 U（$n_G=1.1937\times10^{21}$） \\'),
    # 108. V.1 的「出处：」被改名，出处计数跌破 6
    (APPV, r'\textbf{出处：} 三条运行极限（密度、$\beta$、误差场）的分类与破裂的时间线参照 Wesson',
          r'\textbf{来源：} 三条运行极限（密度、$\beta$、误差场）的分类与破裂的时间线参照 Wesson'),
    # 109. 练习 1 的解答标记被删，解答计数跌破 6
    (APPV, r'\textbf{解答}：由式 \eqref{eq:app-v-tauw}，不锈钢给',
          r'解答：由式 \eqref{eq:app-v-tauw}，不锈钢给'),
    # 110. 表 V.7 的 label 被改名
    (APPV, r'\label{tab:app-v-cq}',
          r'\label{tab:app-v-cqx}'),
    # 111. 第 4 只语义盒（tip-app-v-runaway）的标签被改名——P26 原本只钉 3 只盒，
    #      本轮把 runaway 补进守卫后，这一处由新守卫抓出（补洞的可证伪性证据）
    (APPV, r'{tip-app-v-runaway}', r'{tip-app-v-runawayx}'),
    # 112. 表 W.1 的 label 被改名
    (APPW, r'\label{tab:app-w-life}',
          r'\label{tab:app-w-lifex}'),
    # 113. 第 6 只语义盒的标签被改名
    (APPW, r'{tip-app-w-honest}',
          r'{tip-app-w-honestx}'),
    # 114. 公式 label eq:app-w-tauv 被改名
    (APPW, r'\label{eq:app-w-tauv}',
          r'\label{eq:app-w-tauvx}'),
    # 115. 公式 label eq:app-w-ratios 被改名
    (APPW, r'\label{eq:app-w-ratios}',
          r'\label{eq:app-w-ratiosx}'),
    # 116. alpha 出生能量 5.6076e-13 J 被改成 6.6076e-13 J
    (APPW, r'5.6076\times10^{-13}',
          r'6.6076\times10^{-13}'),
    # 117. W.2 的 v_alpha0/v_Te 由 0.2284 改成 0.3284
    (APPW, r'$0.2284$',
          r'$0.3284$'),
    # 118. 首轨损失份额由 0.4999% 退回 0.499%
    (APPW, r'0.4999\%',
          r'0.499\%'),
    # 119. 快 alpha 密度由 2.1231e17 退回 2.13e17
    (APPW, r'2.1231\times10^{17}',
          r'2.13\times10^{17}'),
    # 120. 能量 e 折时间由 0.1627 s 退回 ~0.16 s
    (APPW, r'$0.1627\ \mathrm{s}$',
          r'$\sim0.16\ \mathrm{s}$'),
    # 121. TAE 频率由 47.09 kHz 改成 48.09 kHz
    (APPW, r'47.09\ \mathrm{kHz}',
          r'48.09\ \mathrm{kHz}'),
    # 122. 临界能量由 295.0 keV 改成 295.1 keV
    (APPW, r'295.0\ \mathrm{keV}',
          r'295.1\ \mathrm{keV}'),
    # 123. 出处标记被换成「来源」，出处计数跌破 6
    (APPW, r'\textbf{出处：}',
          r'\textbf{来源：}'),
    # 124. W.2 的「两条完全不同的推导路径」被改弱
    (APPW, r'两条完全不同的推导路径',
          r'另一条推导路径'),
    # 125. W.4 的尾巴结论句被改弱
    (APPW, r'四成的粒子带着七成半的能量',
          r'一部分粒子带着大部分能量'),
    # 126. W.6 的结论句被改弱
    (APPW, r'快粒子损失不是总能量问题，是分布问题',
          r'快粒子损失是个很小的问题'),
    # 127. 表 W.4 题注的阈值声明被改成给了阈值
    (APPW, r'本表不给阈值',
          r'本表列出各模式的阈值'),
    # 128. W.5 出处段的「不给出任何阈值数字」声明被改成给了数字
    (APPW, r'不给出任何阈值数字',
          r'给出各模式的阈值数字'),
    # 129. W.8 的第三条不复制约定被改
    (APPW, r'不复制任何教材或论文的表格数值',
          r'表格数值参照教材'),
    # 130. tip W.6 的诚实口径声明被改
    (APPW, r'把做不到的地方标出来',
          r'把能算的都算出来'),
    # ---- v1.21.0：带符号物理量（check_P28_signed_quantities）----
    # 131. 附录 S：IPB98(y,2) 里 P 的指数被翻成正号（负号就是结论本身）
    (APPS, r'$P$ & $-0.69$ & $0.68$ &',
           r'$P$ & $0.69$ & $0.68$ &'),
    # 132. 附录 U：同一条定标律的指数被翻成正号
    (APPU, r'$P$ 的指数是 $-0.69$：',
           r'$P$ 的指数是 $0.69$：'),
    # 133. 附录 T：磁岛 Δ′ 的符号被翻正（符号决定稳定性）
    (APPT, r'$-7.680\times10^{-3}$',
           r'$7.680\times10^{-3}$'),
    # 134. 正文：悬浮电位的负号被去掉（物理上必须为负）
    (MAIN, r'氩约为 $-4.7\,T_e$',
           r'氩约为 $4.7\,T_e$'),
    # 135. 第 12 章补充：Langmuir 探针 I-V 表的离子电流被翻成正号
    (CH12, r'$-20$ & $-0.51$',
           r'$-20$ & $0.51$'),
    # 136. 附录 Q：探针曲线表的偏压列负号被去掉
    (APPQ, r'$-1.0$ & 7.1653 & 6.1653',
           r'$1.0$ & 7.1653 & 6.1653'),
    # ---- v1.22.0：附录 X（辐射输运与辐射损失，check_P29）----
    # 137. X.4 表格 label 被改名
    (APPX, r'\label{tab:app-x-cyclotron}',
          r'\label{tab:app-x-cyclo}'),
    # 138. tip X.5 的语义盒标签被改名
    (APPX, r'{tip-app-x-honest}',
          r'{tip-app-x-modest}'),
    # 139. X.6 公式 label 被改名
    (APPX, r'\label{eq:app-x-taurad}',
          r'\label{eq:app-x-taur}'),
    # 140. X.1 的「出处：」被写成「来源：」
    (APPX, r'\textbf{出处：}三种机制的划分',
          r'\textbf{来源：}三种机制的划分'),
    # 141. 解答 4 被写成「答案 4」
    (APPX, r'\textbf{解答 4}',
          r'\textbf{答案 4}'),
    # 142. 表 X.4 题注的诚实口径被改
    (APPX, r'本附录不给输运方程的解',
          r'本附录给出输运方程的解'),
    # 143. 表 X.5 题注的定性声明被改
    (APPX, r'不给具体数值，因为本附录没有引入原子物理数据库',
          r'不给具体数值，因为本附录引入了原子物理数据库'),
    # 144. 表 X.1 回旋辐射行：890.01 被四舍五入成 890.0
    (APPX, r'$890.01\ \mathrm{MW}$',
          r'$890.0\ \mathrm{MW}$'),
    # 145. 表 X.4 光学薄上限行：1.5993 被写成 1.6
    (APPX, r'$=1.5993\,P_{\mathrm{fus}}$，不可能',
          r'$=1.6\,P_{\mathrm{fus}}$，不可能'),
    # 146. 表 X.3 光子波长被四舍五入成 0.13 nm
    (APPX, r'$0.134853\ \mathrm{nm}$',
          r'$0.13\ \mathrm{nm}$'),
    # 147. 表 X.3 相对论修正行被四舍五入
    (APPX, r'$27.9716\ \mathrm{MW}$ & 提高 $3.598\%$',
          r'$28.0\ \mathrm{MW}$ & 提高 $3.6\%$'),
    # 148. 表 X.4 再吸收倍数被四舍五入成 5400
    (APPX, r'再吸收压低倍数 & $5406.41$',
          r'再吸收压低倍数 & $5400$'),
    # 149. 表 X.6 自加热余量被四舍五入成 4.0
    (APPX, r'$3.99897$ & --- & 约 4 倍，不是 100 倍',
          r'$4.0$ & --- & 约 4 倍，不是 100 倍'),
    # 150. 表 X.6 辐射冷却时间被四舍五入成 13.3 s
    (APPX, r'$13.3443\ \mathrm{s}$ & --- & 比 $\tau_{\mathrm{CQ}}$ 慢 $173.511$ 倍',
          r'$13.3\ \mathrm{s}$ & --- & 比 $\tau_{\mathrm{CQ}}$ 慢 $173.5$ 倍'),
    # 151. 表 X.6 壁面轫致热流被四舍五入成 55 kW/m^2
    (APPX, r'$55.1547\ \mathrm{kW\,m^{-2}}$ & --- & X 射线',
          r'$55\ \mathrm{kW\,m^{-2}}$ & --- & X 射线'),
    # 152. 解答 5 的平均自由程被写成 6.8e5 m
    (APPX, r'6.8301\times10^{5}',
          r'6.8\times10^{5}'),
    # 153. 解答 3 的钨浓度份额被写成 0.032%
    (APPX, r'$n_W/n_e = 0.0324441\%$',
          r'$n_W/n_e = 0.032\%$'),
    # 154. 解答 4 的亮度温度被取整成 1508 keV
    (APPX, r'$T = 1507.93\ \mathrm{keV}$',
          r'$T = 1508\ \mathrm{keV}$'),
    # 155. X.1 的磁能密度被四舍五入成 1.118e7
    (APPX, r'1.11767\times10^{7}',
          r'1.118\times10^{7}'),
    # 156. X.1 的 beta_perp^2 被写成 0.036
    (APPX, r'=0.035984$',
          r'=0.036$'),
    # 157. 解答 2 的功率密度被四舍五入成 47.85
    (APPX, r'47.8519\ \mathrm{kW\,m^{-3}}',
          r'47.85\ \mathrm{kW\,m^{-3}}'),
    # 158. 表 X.4 黑体上限份额被写成 0.61%
    (APPX, r'黑体上限占轫致辐射 & $0.60971\%$',
          r'黑体上限占轫致辐射 & $0.61\%$'),
    # 159. 解答 3 的钨密度被四舍五入成 3.24e16
    (APPX, r'3.24441\times10^{16}',
          r'3.24\times10^{16}'),
    # 160. 解答 3 引用的崩塌阈值被写成 0.13%
    (APPX, r'$0.133742\%$ 低四倍多',
          r'$0.13\%$ 低四倍多'),
    # 161. 表 X.1 轫致辐射行整行数值被微调（容差抓不到，靠整行绑定）
    (APPX, r'轫致辐射 & $32.444\ \mathrm{kW\,m^{-3}}$ & $27.00\ \mathrm{MW}$ &',
          r'轫致辐射 & $32.45\ \mathrm{kW\,m^{-3}}$ & $27.00\ \mathrm{MW}$ &'),
    # ---- v1.23.0（附录 Y 强场与激光等离子体）新增 18 处 ----
    # 162. 表 Y.2 的 label 被改名
    (APPY, r'\label{tab:app-y-wake}', r'\label{tab:app-y-waves}'),
    # 163. tip Y.3「透明不等于穿透」的语义盒标签被改名
    (APPY, r'{tip-app-y-transparency}', r'{tip-app-y-transparent}'),
    # 164. 式 Y.9 的 label 被改名
    (APPY, r'\label{eq:app-y-ewb}', r'\label{eq:app-y-wavebreak}'),
    # 165. 第一处「出处：」被写成「来源：」，出处计数跌破 6
    (APPY, r'\textbf{出处：}', r'\textbf{来源：}'),
    # 166. 「解答 4」被写成「答案 4」，解答计数跌破 6
    (APPY, r'\textbf{解答 4}', r'\textbf{答案 4}'),
    # 167. tip Y.1 的标题被改字（整行绑定抓）
    (APPY, r'\begin{tip}{场的强度要用无量纲数说话}{tip-app-y-field}',
          r'\begin{tip}{场的强度要用无量纲数字说话}{tip-app-y-field}'),
    # 168. tip Y.2 的标题被改字（整行绑定抓）
    (APPY, r'\begin{tip}{失相不是误差，是原理}{tip-app-y-dephasing}',
          r'\begin{tip}{失相不是误差而是原理}{tip-app-y-dephasing}'),
    # 169. 表 Y.1 激光电场：2.74309 -> 2.7431（差值 3.6e-6，落在 rel=1e-3 内）
    (APPY, r'激光电场 & $2.74309\times10^{14}\ \mathrm{V\,m^{-1}}$',
          r'激光电场 & $2.7431\times10^{14}\ \mathrm{V\,m^{-1}}$'),
    # 170. 表 Y.2 波破场：96.1592 -> 96.159（差值 2.6e-6，落在容差内）
    (APPY, r'波破场 $E_{WB}$ & $96.1592\ \mathrm{GV\,m^{-1}}$',
          r'波破场 $E_{WB}$ & $96.159\ \mathrm{GV\,m^{-1}}$'),
    # 171. 表 Y.3 理想增益：5.5929 -> 5.5928（差值 1.8e-5，落在容差内）
    (APPY, r'理想增益 $\Delta W$ & $5.5929\ \mathrm{GeV}$',
          r'理想增益 $\Delta W$ & $5.5928\ \mathrm{GeV}$'),
    # 172. 表 Y.1 归一化振幅：68.3486 -> 68.3487（差值 1.5e-6，落在容差内）
    (APPY, r'归一化振幅 & $a_0 = 68.3486$',
          r'归一化振幅 & $a_0 = 68.3487$'),
    # 173. 表 Y.3 脉冲能量：235.619 -> 235.62（差值 4.2e-6，落在容差内）
    (APPY, r'脉冲能量 & $235.619\ \mathrm{J}$',
          r'脉冲能量 & $235.62\ \mathrm{J}$'),
    # 174. 表 Y.3 自聚焦功率：30.31 -> 30.30（差值 3.3e-4，落在容差内）
    (APPY, r'自聚焦功率 $P_c$ & $30.31\ \mathrm{TW}$',
          r'自聚焦功率 $P_c$ & $30.30\ \mathrm{TW}$'),
    # 175. 表 Y.6 聚变能：960.2 -> 906.2 MJ（数字易位，容差外，数值用例抓）
    (APPY, r'聚变能 $E_{\mathrm{fus}}$ & $960.2\ \mathrm{MJ}$',
          r'聚变能 $E_{\mathrm{fus}}$ & $906.2\ \mathrm{MJ}$'),
    # 176. 表 Y.2 临界密度：1.742 -> 1.74（差 1.1e-3，刚出容差，数值用例抓）
    (APPY, r'临界密度 $n_c$ & $1.742\times10^{21}\ \mathrm{cm^{-3}}$',
          r'临界密度 $n_c$ & $1.74\times10^{21}\ \mathrm{cm^{-3}}$'),
    # 177. 表 Y.4 固体折射率：57.40661 -> 57.40671（差值 1.8e-6，落在容差内）
    (APPY, r'固体 $n_e/n_c$ & $57.40661$',
          r'固体 $n_e/n_c$ & $57.40671$'),
    # 178. 表 Y.5 辐射压强：3.33564e17 -> 3.34e17（数值用例 + 整行绑定）
    (APPY, r'辐射压强 $I/c$ & $3.33564\times10^{17}\ \mathrm{Pa}$',
          r'辐射压强 $I/c$ & $3.34\times10^{17}\ \mathrm{Pa}$'),
    # 179. 表 Y.4 阈值强度：7.05447 -> 7.05（落在容差内，靠 lacks 守卫 + 整行绑定）
    (APPY, r'阈值强度 & $7.05447\times10^{21}\ \mathrm{W\,cm^{-2}}$',
          r'阈值强度 & $7.05\times10^{21}\ \mathrm{W\,cm^{-2}}$'),
    # ---- v1.24.0（附录 Z 聚变堆工程）新增 18 处 ----
    # 180. 表 Z.3 的 label 被改名
    (APPZ, r'\label{tab:app-z-tritium}', r'\label{tab:app-z-tritiumcycle}'),
    # 181. tip Z.6 的语义盒标签被改名
    (APPZ, r'{tip-app-z-availability}', r'{tip-app-z-avail}'),
    # 182. 式 Z.12 的 label 被改名
    (APPZ, r'\label{eq:app-z-tbrbreak}', r'\label{eq:app-z-tbrlimit}'),
    # 183. 第一处「出处：」被写成「来源：」，出处计数跌破 6
    (APPZ, r'\textbf{出处：}', r'\textbf{来源：}'),
    # 184. 「解答 3」被写成「答案 3」，解答计数跌破 6
    (APPZ, r'\textbf{解答 3}', r'\textbf{答案 3}'),
    # 185. tip Z.1 的标题被改字（整行绑定抓）
    (APPZ, r'\begin{tip}{工程账是乘法，不是加法}{tip-app-z-chain}',
          r'\begin{tip}{工程账是相乘，不是相加}{tip-app-z-chain}'),
    # 186. tip Z.2 的标题被改字（整行绑定抓）
    (APPZ, r'\begin{tip}{难受的不是热，是 dpa}{tip-app-z-neutron}',
          r'\begin{tip}{难受的不是热，而是 dpa}{tip-app-z-neutron}'),
    # 187. 表 Z.1 包层热功率：639.97 -> 639.9（差值 1.1e-4，落在 rel=1e-3 内）
    (APPZ, r'包层热功率 & $639.97\ \mathrm{MW}$', r'包层热功率 & $639.9\ \mathrm{MW}$'),
    # 188. 表 Z.1 毛电功率：255.99 -> 255.9（差值 3.5e-4）
    (APPZ, r'毛电功率 & $255.99\ \mathrm{MW}$', r'毛电功率 & $255.9\ \mathrm{MW}$'),
    # 189. 表 Z.1 净电功率：217.59 -> 217.6（差值 4.6e-5）
    (APPZ, r'净电功率 & $217.59\ \mathrm{MW}$', r'净电功率 & $217.6\ \mathrm{MW}$'),
    # 190. 表 Z.2 中子壁负载：0.90944 -> 0.9094（差值 4.4e-5）
    (APPZ, r'中子壁负载 & $0.90944\ \mathrm{MW\,m^{-2}}$',
          r'中子壁负载 & $0.9094\ \mathrm{MW\,m^{-2}}$'),
    # 191. 表 Z.2 包层寿命：10.996 -> 10.99（差值 5.5e-4）
    (APPZ, r'包层寿命 & $10.996$ 年', r'包层寿命 & $10.99$ 年'),
    # 192. 表 Z.3 氚燃耗率：31.191 -> 31.19（差值 3.2e-5）
    (APPZ, r'氚燃耗率 & $31.191\ \mathrm{kg\,year^{-1}}$',
          r'氚燃耗率 & $31.19\ \mathrm{kg\,year^{-1}}$'),
    # 193. 表 Z.3 自持 TBR 下限：1.0036 -> 1.004（差值 4.0e-4）
    (APPZ, r'自持 TBR 下限 & $1.0036$', r'自持 TBR 下限 & $1.004$'),
    # 194. 表 Z.4 屏蔽厚度：2.0723 -> 2.072（差值 1.4e-4）
    (APPZ, r'屏蔽厚度（$10^{-6}$） & $2.0723\ \mathrm{m}$',
          r'屏蔽厚度（$10^{-6}$） & $2.072\ \mathrm{m}$'),
    # 195. 表 Z.5 质量流量：615.36 -> 615.4（差值 6.5e-5）
    (APPZ, r'质量流量 & $615.36\ \mathrm{kg\,s^{-1}}$',
          r'质量流量 & $615.4\ \mathrm{kg\,s^{-1}}$'),
    # 196. 表 Z.6 年发电量：1429.6 -> 1429（差值 4.2e-4）
    (APPZ, r'年发电量 & $1429.6\ \mathrm{GWh}$', r'年发电量 & $1429\ \mathrm{GWh}$'),
    # 197. 表 Z.6 折合标煤：0.48606 -> 0.49（容差外，且触发 lacks 守卫）
    (APPZ, r'折合标煤 & $0.48606$ 百万吨', r'折合标煤 & $0.49$ 百万吨'),
    # ---- v1.25.0（附录 AA 替代位形）新增 18 处 ----
    # 198. 表 AA.3 的 label 被改名
    (APPAA, r'\label{tab:app-aa-transport}', r'\label{tab:app-aa-transport2}'),
    # 199. tip AA.4 的语义盒标签被改名
    (APPAA, r'{tip-app-aa-st}', r'{tip-app-aa-st2}'),
    # 200. 式 AA.8 的 label 被改名
    (APPAA, r'\label{eq:app-aa-eff32}', r'\label{eq:app-aa-eff32b}'),
    # 201. 第一处「出处：」被写成「来源：」，出处计数跌破 6
    (APPAA, r'\textbf{出处：}', r'\textbf{来源：}'),
    # 202. 「解答 5」被写成「答案 5」，解答计数跌破 6
    (APPAA, r'\textbf{解答 5}', r'\textbf{答案 5}'),
    # 203. tip AA.1 的标题被改字（整行绑定抓）
    (APPAA, r'\begin{tip}{先问「它拿什么换什么」}{tip-app-aa-map}',
          r'\begin{tip}{先问「拿什么换什么」}{tip-app-aa-map}'),
    # 204. tip AA.3 的标题被改字（整行绑定抓）
    (APPAA, r'\begin{tip}{不同区间的数不能直接相除}{tip-app-aa-eps}',
          r'\begin{tip}{不同区间的数不能相除}{tip-app-aa-eps}'),
    # 205. 表 AA.1 D1 压强：1.473e5 -> 1.474e5（差值 6.8e-4，落在容差内）
    (APPAA, r'压强 $p$ & $1.473\times10^{5}\ \mathrm{Pa}$',
          r'压强 $p$ & $1.474\times10^{5}\ \mathrm{Pa}$'),
    # 206. 表 AA.1 D1 逆环径比：0.32258 -> 0.3226（差值 6.2e-5）
    (APPAA, r'逆环径比 $\varepsilon$ & $0.32258$ & $0.66667$ & $0.096364$ \\',
          r'逆环径比 $\varepsilon$ & $0.3226$ & $0.66667$ & $0.096364$ \\'),
    # 207. 表 AA.1 D1 归一化比压：1.5523 -> 1.5524（差值 6.4e-5）
    (APPAA, r'归一化比压 $\beta_N$ & $1.5523$ & $1.2483$ & $1.3659$ \\',
          r'归一化比压 $\beta_N$ & $1.5524$ & $1.2483$ & $1.3659$ \\'),
    # 208. 表 AA.3 香蕉区因子：0.56796 -> 0.5680（差值 7.0e-5）
    (APPAA, r'$\sqrt{\varepsilon_{D1}}$ & $0.56796$ & 托卡马克几何 & 香蕉区放大因子 \\',
          r'$\sqrt{\varepsilon_{D1}}$ & $0.5680$ & 托卡马克几何 & 香蕉区放大因子 \\'),
    # 209. 表 AA.3 优化品质：0.13491 -> 0.1349（差值 7.4e-5）
    (APPAA, r'$\varepsilon_{\mathrm{eff}}/\varepsilon$ & $0.13491$ & 优化品质 & 只有几何值的 $13.5\%$ \\',
          r'$\varepsilon_{\mathrm{eff}}/\varepsilon$ & $0.1349$ & 优化品质 & 只有几何值的 $13.5\%$ \\'),
    # 210. 表 AA.4 边缘极向场：0.64516 -> 0.6452 T（差值 6.2e-5）
    (APPAA, r'边缘极向场 $B_p$ & --- & $0.64516\ \mathrm{T}$ & $\mu_0I_p/(2\pi a)$ \\',
          r'边缘极向场 $B_p$ & --- & $0.6452\ \mathrm{T}$ & $\mu_0I_p/(2\pi a)$ \\'),
    # 211. 表 AA.4 自举份额：0.078989 -> 0.07899（差值 1.3e-5）
    (APPAA, r'自举份额 $f_{bs}$ & 几个 $\%$ & $0.078989$ & $C_{bs}\sqrt{\varepsilon}\beta_p$ \\',
          r'自举份额 $f_{bs}$ & 几个 $\%$ & $0.07899$ & $C_{bs}\sqrt{\varepsilon}\beta_p$ \\'),
    # 212. 表 AA.4 D1 中心柱：4.2 -> 4.21 m（差值 2.4e-3，容差外，数值用例抓）
    (APPAA, r'中心柱空间 $R-a$ & $4.2\ \mathrm{m}$ & $0.31\ \mathrm{m}$ & 屏蔽需要 $2.0723\ \mathrm{m}$ \\',
          r'中心柱空间 $R-a$ & $4.21\ \mathrm{m}$ & $0.31\ \mathrm{m}$ & 屏蔽需要 $2.0723\ \mathrm{m}$ \\'),
    # 213. 表 AA.6 ISS04 的 tau_E：0.142983 -> 0.14298 s（差值 2.1e-5）
    (APPAA, r'$\tau_E^{\mathrm{ISS04}}$ & $0.142983\ \mathrm{s}$ & 式 \eqref{eq:app-aa-iss04} & 研究装置量级 \\',
          r'$\tau_E^{\mathrm{ISS04}}$ & $0.14298\ \mathrm{s}$ & 式 \eqref{eq:app-aa-iss04} & 研究装置量级 \\'),
    # 214. 表 AA.6 三重积：2.28774e19 -> 2.2877e19（差值 1.7e-5）
    (APPAA, r'三重积 & $2.28774\times10^{19}$ & $nT\tau_E$ & 差 $135.51$ 倍 \\',
          r'三重积 & $2.2877\times10^{19}$ & $nT\tau_E$ & 差 $135.51$ 倍 \\'),
    # 215. 表 AA.3 有效纹波：0.013 -> 0.0132（差值 1.5e-2，容差外，数值用例抓）
    (APPAA, r'$\varepsilon_{\mathrm{eff}}$ & $0.013$ & 位形优化 & 声明值（优化位形） \\',
          r'$\varepsilon_{\mathrm{eff}}$ & $0.0132$ & 位形优化 & 声明值（优化位形） \\'),
    # ---- v1.26.0 收尾轮（附录 I 第 17 节「读一篇论文」/ 门禁 P33）新增 30 处 ----
    # 216. 节标题「读一篇论文」被写成「读一篇文章」
    (APPI, r'\section{读一篇论文：三遍法、提问清单与复现实验}',
          r'\section{读一篇文章：三遍法、提问清单与复现实验}'),
    # 217. 节 label 被改名
    (APPI, r'\label{sec:app-i-paper}', r'\label{sec:app-i-papers}'),
    # 218. 三遍法表 label 被改名
    (APPI, r'\label{tab:app-i-paperpass}', r'\label{tab:app-i-paperpas}'),
    # 219. 六问表 label 被改名
    (APPI, r'\label{tab:app-i-paperquestions}', r'\label{tab:app-i-paperquestionx}'),
    # 220. 语义盒**标题**被改字（必须排在 221 之前：221 的锚点是本条锚点的子串，
    #      先改标签会让本条锚点找不到——夹具按列表顺序 replace(,1)，顺序即依赖）
    (APPI, r'\begin{tip}{「读懂」的可操作定义}{tip-app-i-paper}',
          r'\begin{tip}{「读懂」的主观定义}{tip-app-i-paper}'),
    # 221. 语义盒标签被改名
    (APPI, '{tip-app-i-paper}', '{tip-app-i-pape}'),
    # 222. 小节一标题改字
    (APPI, r'\subsection*{一、三遍法：每一遍都有「通过判据」}',
          r'\subsection*{一、三遍法：每一遍都有「大致判据」}'),
    # 223. 小节二标题改字（把图当成「主文」-> 「插图」）
    (APPI, r'\subsection*{二、第二遍的核心动作：把图当成「主文」}',
          r'\subsection*{二、第二遍的核心动作：把图当成「插图」}'),
    # 224. 小节三标题改字（第四遍 -> 第五遍）
    (APPI, r'\subsection*{三、第四遍其实是「算一遍」：复现实验}',
          r'\subsection*{三、第五遍其实是「算一遍」：复现实验}'),
    # 225. 小节四标题改字（三个数 -> 四个数）
    (APPI, r'\subsection*{四、文献笔记：一行问题、三个数、一个步骤、一句反对}',
          r'\subsection*{四、文献笔记：一行问题、四个数、一个步骤、一句反对}'),
    # 226. 小节五标题改字（三条 -> 两条常见错误）
    (APPI, r'\subsection*{五、三条常见错误}', r'\subsection*{五、两条常见错误}'),
    # 227. 表 1 第一遍：五个 C 退回三个 C（同时触发字面绑定与 lacks 守卫，单条 FAIL=2）
    (APPI, '五个 C：类别、背景、正确性、贡献、清晰度', '三个 C：类别、背景、正确性'),
    # 228. 表 1 第三遍：能复现其中一个数 -> 两个数
    (APPI, r'\textbf{能复现其中一个数}（见第三节）',
          r'\textbf{能复现其中两个数}（见第三节）'),
    # 229. 表 1 第二遍行注入「第二遍读完就算读懂」（同时触发字面绑定与 lacks，单条 FAIL=2）
    (APPI, r'第二遍 & 图、表、坐标轴、误差棒、关键公式的\textbf{结构}（不追推导）',
          r'第二遍 & 图、表、坐标轴、误差棒、关键公式的\textbf{结构}（第二遍读完就算读懂）'),
    # 230. 表 1 第二遍行的产出与判据两格整行改字（三到五个量级 -> 一到两个）
    (APPI, r'三到五个关键量级 + 一份「没看懂的术语」清单 & 能说出每个图的横纵轴是什么、结论靠哪张图支撑 \\',
          r'一到两个关键量级 + 一份「没看懂的术语」清单 & 能说出每个图的标题是什么 \\'),
    # 231. 表 2「我能算出同一个数吗」行：M--AA 退回 M--X（同时触发 lacks，单条 FAIL=2）
    (APPI, r'我能算出同一个数吗？ & \textbf{读懂的判据} & 附录 M--AA 的 \texttt{calc} 脚本 \\',
          r'我能算出同一个数吗？ & \textbf{读懂的判据} & 附录 M--X 的 \texttt{calc} 脚本 \\'),
    # 232. 表 2 最后一行：读者/同行 方向翻转
    (APPI, r'我不同意哪一点？ & 从「读者」变成「同行」 & 附录 K 的例题与解答对照 \\',
          r'我不同意哪一点？ & 从「同行」变成「读者」 & 附录 K 的例题与解答对照 \\'),
    # 233. 表 2 第一行：「它解决什么问题」被换成「它引用了什么」
    (APPI, r'它解决什么问题？ & 判断与自己的目标是否相关 & 附录 J 的逐章对照表 \\',
          r'它引用了什么？ & 判断参考文献是否够新 & 附录 J 的逐章对照表 \\'),
    # 234. 「读懂」的可操作定义句：一个数 -> 一些数
    (APPI, r'\textbf{能用本书的脚本复现出论文里的一个数，',
          r'\textbf{能用本书的脚本复现出论文里的一些数，'),
    # 235. 笔记模板的可复现步骤：代入 -> 参考
    (APPI, r'用 \texttt{tmp/app\_\_calc.py} 代入', r'用 \texttt{tmp/app\_\_calc.py} 参考'),
    # 236. 出处：Keshav 年份 2007 -> 2006
    (APPI, 'Keshav（2007）', 'Keshav（2006）'),
    # 237. 出处：Keshav 卷期页 83--84 -> 83--85
    (APPI, r'ACM SIGCOMM Computer Communication Review \textbf{37}(3), 83--84',
          r'ACM SIGCOMM Computer Communication Review \textbf{37}(3), 83--85'),
    # 238. 出处：Booth 等 The Craft of Research 第 3 版 -> 第 2 版
    (APPI, r'\emph{The Craft of Research}（第 3 版，2008）',
          r'\emph{The Craft of Research}（第 2 版，2008）'),
    # 239. 出处：Wilson 等书名 Best -> Good
    (APPI, r'\emph{Best Practices for Scientific Computing}',
          r'\emph{Good Practices for Scientific Computing}'),
    # 240. 常见错误之一：把摘要当结论 -> 把摘要当摘要
    (APPI, r'\textbf{把摘要当结论}', r'\textbf{把摘要当摘要}'),
    # 241. 复现不上时的口径句：一次真正的收获 -> 一次意外的收获
    (APPI, '一次真正的收获——本书自己的勘误', '一次意外的收获——本书自己的勘误'),
    # 242. 注入「感觉懂了」（改表 1 题注，触发 lacks）
    (APPI, '读没读懂不由「感觉」判断', '读没读懂感觉懂了就行'),
    # 243. 注入「复现不上就算了」（触发 lacks）
    (APPI, r'\textbf{复现不上怎么办？记录差异，而不是忽略差异。}',
          r'\textbf{复现不上就算了，不必记录差异。}'),
    # 244. 注入「引用次数越多越正确」（触发 lacks）
    (APPI, r'\textbf{把引用次数当正确性}', r'\textbf{引用次数越多越正确}'),
    # 245. 口径反转：「可以随时停止」被反写成「必须读完三遍」（同时触发字面绑定与 lacks，
    #      单条 FAIL=2）——本节的核心方法论口径，v1.26.0 发布时完全没有绑定，是探针抓出的洞
    (APPI, r'三遍法最重要的性质是\textbf{可以随时停止}：',
          r'三遍法最重要的性质是\textbf{必须读完三遍}：'),
    # 246. 结构：表 2 的 label 改名（tab:app-ab-solarwind -> wind），P34 的表标签检查应报缺
    (APPAB, r'\label{tab:app-ab-solarwind}',
           r'\label{tab:app-ab-wind}'),
    # 247. 结构：第 3 只语义盒的标签改名（tip-app-ab-collisionless -> shock）
    (APPAB, r'{tip-app-ab-collisionless}',
           r'{tip-app-ab-shock}'),
    # 248. 结构：式 (AB.17) 的 label 改名（eq:app-ab-petschek -> pet）
    (APPAB, r'\label{eq:app-ab-petschek}',
           r'\label{eq:app-ab-pet}'),
    # 249. 结构：第一处「出处：」声明被改成「来源：」，出处计数应跌破 6
    (APPAB, r'\textbf{出处：}本节的分区口径与尺度取值参照',
           r'\textbf{来源：}本节的分区口径与尺度取值参照'),
    # 250. 结构：练习解答 4 的编号词被改（解答 -> 答案）
    (APPAB, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 251. 结构：第 1 只语义盒标题改字（数量级 -> 能标）
    (APPAB, r'四个尺度不是四个数量级，是四种物理',
           r'四个尺度不是四个能标，是四种物理'),
    # 252. 结构：第 6 只语义盒标题改字（带 -> 环）
    (APPAB, r'辐射带不是「带」，是磁矩守恒的必然结果',
           r'辐射带不是「环」，是磁矩守恒的必然结果'),
    # 253. 数据：表 2 电子等离子体频率 1.26147e5 -> 1.26148e5
    (APPAB, r'电子等离子体频率 & $1.26147\times10^{5}$',
           r'电子等离子体频率 & $1.26148\times10^{5}$'),
    # 254. 数据：表 2 离子等离子体频率 2943.89 -> 2943.9
    (APPAB, r'离子等离子体频率 & $2943.89$',
           r'离子等离子体频率 & $2943.9$'),
    # 255. 数据：表 2 声速 43.769 -> 43.77
    (APPAB, r'声速 & $43.769$',
           r'声速 & $43.77$'),
    # 256. 数据：表 2 等离子体 beta 0.80534 -> 0.8053
    (APPAB, r'等离子体 $\beta$ & $0.80534$',
           r'等离子体 $\beta$ & $0.8053$'),
    # 257. 数据：表 3 太阳风动压 1.3381e-9 -> 1.3382e-9
    (APPAB, r'太阳风动压 $\rho v^2$ & $1.3381\times10^{-9}$',
           r'太阳风动压 $\rho v^2$ & $1.3382\times10^{-9}$'),
    # 258. 数据：表 3 磁层顶距离 8.1158 -> 8.1159
    (APPAB, r'磁层顶距离 & $8.1158$',
           r'磁层顶距离 & $8.1159$'),
    # 259. 数据：表 3 太阳风输入功率 34.126 -> 34.127 GW
    (APPAB, r'太阳风输入功率 & $34.126$',
           r'太阳风输入功率 & $34.127$'),
    # 260. 数据：表 4 阿尔芬电导 16.316 -> 16.317 S
    (APPAB, r'$16.316\ \mathrm{S}$',
           r'$16.317\ \mathrm{S}$'),
    # 261. 数据：表 5 伦德奎斯特数 1.86362e10 -> 1.86363e10
    (APPAB, r'伦德奎斯特数 $S$ & $1.86362\times10^{10}$',
           r'伦德奎斯特数 $S$ & $1.86363\times10^{10}$'),
    # 262. 数据：表 5 两模型之比 2266.9 -> 2267.0（同时触发 lacks 守卫，单条 FAIL=2）
    (APPAB, r'两者之比 & $2266.9$',
           r'两者之比 & $2267.0$'),
    # 263. 数据：表 6 外辐射带电子回旋半径 47.432 -> 47.433 km
    (APPAB, r'$r_g=47.432\ \mathrm{km}$',
           r'$r_g=47.433\ \mathrm{km}$'),
    # 264. 结构：表 1 的 label 改名（tab:app-ac-resolution -> grid），P35 的表标签检查应报缺
    (APPAC, r'\label{tab:app-ac-resolution}',
           r'\label{tab:app-ac-grid}'),
    # 265. 结构：第 2 只语义盒的标签改名（tip-app-ac-conserve -> energy）
    (APPAC, r'{tip-app-ac-conserve}',
           r'{tip-app-ac-energy}'),
    # 266. 结构：式 (AC.17) 的 label 改名（eq:app-ac-trunc -> cut）
    (APPAC, r'\label{eq:app-ac-trunc}',
           r'\label{eq:app-ac-cut}'),
    # 267. 结构：第一处「出处：」声明被改成「来源：」
    (APPAC, r'\textbf{出处：}本节的分辨率判据',
           r'\textbf{来源：}本节的分辨率判据'),
    # 268. 结构：练习解答 4 的编号词被改（解答 -> 答案），解答计数应跌破练习数
    (APPAC, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 269. 结构：第 1 只语义盒标题改字（分辨率 -> 精度）
    (APPAC, r'先算分辨率，再谈算法',
           r'先算精度，再谈算法'),
    # 270. 结构：第 2 只语义盒标题改字（算得好 -> 算得准）
    (APPAC, r'守恒量不是「算得好」，是格式的性质',
           r'守恒量不是「算得准」，是格式的性质'),
    # 271. 数据：式 (AC.1) 德拜长度 7.12806e-5 -> 7.12807e-5
    (APPAC, r'=7.12806\times10^{-5}\ \mathrm{m}=71.28',
           r'=7.12807\times10^{-5}\ \mathrm{m}=71.28'),
    # 272. 数据：式 (AC.2) 等离子体频率 5.64146e11 -> 5.64147e11
    (APPAC, r'5.64146\times10^{11}',
           r'5.64147\times10^{11}'),
    # 273. 数据：式 (AC.2) 回旋频率 9.32175e11 -> 9.32176e11
    (APPAC, r'9.32175\times10^{11}',
           r'9.32176\times10^{11}'),
    # 274. 数据：式 (AC.4) 德拜球内粒子数 3.62171e7 -> 3.62172e7
    (APPAC, r'3.62171\times10^{7}',
           r'3.62172\times10^{7}'),
    # 275. 数据：式 (AC.5) 碰撞时间 8.92144e-5 -> 8.92145e-5
    (APPAC, r'8.92144\times10^{-5}',
           r'8.92145\times10^{-5}'),
    # 276. 数据：式 (AC.9) 每步相位误差 0.00166667 -> 0.00166668
    (APPAC, r'=0.00166667\ \mathrm{rad},',
           r'=0.00166668\ \mathrm{rad},'),
    # 277. 数据：式 (AC.14) 粒子步总数 1.65236e11 -> 1.65237e11
    (APPAC, r'\times8261.82=1.65236\times10^{11}',
           r'\times8261.82=1.65237\times10^{11}'),
    # 278. 数据：表 4 每格 20 粒子的相对噪声 22.361% -> 22.362%（同时触发 lacks 守卫）
    (APPAC, r'$20$ & $22.361\%$ & $2.2361\%$',
           r'$20$ & $22.362\%$ & $2.2361\%$'),
    # 279. 数据：式 (AC.16) FFT 加速比 50171.7 -> 50171.8（同时触发 lacks 守卫）
    (APPAC, r'50171.7',
           r'50171.8'),
    # 280. 数据：式 (AC.21) 256^3 网格内存 1.07374e9 -> 1.07375e9
    (APPAC, r'1.07374\times10^{9}',
           r'1.07375\times10^{9}'),
    # 281. 数据：练习 1 的德拜长度 7.43394e-5 -> 7.43395e-5
    (APPAC, r'7.43394\times10^{-5}',
           r'7.43395\times10^{-5}'),
    # 282. 结构：表 2 的 label 改名（tab:app-ad-pedestal -> ped），P36 的表标签检查应报缺
    (APPAD, r'\label{tab:app-ad-pedestal}',
           r'\label{tab:app-ad-ped}'),
    # 283. 结构：第 3 只语义盒的标签改名（tip-app-ad-two-branches -> branches）
    (APPAD, r'{tip-app-ad-two-branches}',
           r'{tip-app-ad-branches}'),
    # 284. 结构：式 (AD.22) 的 label 改名（eq:app-ad-detach -> det）
    (APPAD, r'\label{eq:app-ad-detach}',
           r'\label{eq:app-ad-det}'),
    # 285. 结构：小节 AD.6 的 sec label 改名（sec:app-ad-sol -> scrapeoff）
    (APPAD, r'\label{sec:app-ad-sol}',
           r'\label{sec:app-ad-scrapeoff}'),
    # 286. 结构：第一处「出处：」声明被改成「来源：」
    (APPAD, r'\textbf{出处：}L-H 转换的实验发现参照',
           r'\textbf{来源：}L-H 转换的实验发现参照'),
    # 287. 结构：练习解答 4 的编号词被改（解答 -> 答案），解答计数应跌破练习数
    (APPAD, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 288. 结构：第 2 只语义盒标题改字（自由参数 -> 固定参数）
    (APPAD, r'台基宽度不是自由参数',
           r'台基宽度不是固定参数'),
    # 289. 数据：台基壳层面积 489.532 m^2 -> 489.533 m^2（P36 精确绑定应报 FAIL）
    (APPAD, r'489.532\ \mathrm{m^2}',
           r'489.533\ \mathrm{m^2}'),
    # 290. 数据：密度减半时的阈值 38.4821 MW -> 38.4822 MW
    (APPAD, r'38.4821\ \mathrm{MW}',
           r'38.4822\ \mathrm{MW}'),
    # 291. 数据：单次 ELM 损失 2.59980e7 J -> 2.59981e7 J
    (APPAD, r'2.59980\times10^{7}\ \mathrm{J}',
           r'2.59981\times10^{7}\ \mathrm{J}'),
    # 292. 数据：ELM 能量密度 6.67373e5 J/m^2 -> 6.67374e5 J/m^2
    (APPAD, r'6.67373\times10^{5}\ \mathrm{J\,m^{-2}}',
           r'6.67374\times10^{5}\ \mathrm{J\,m^{-2}}'),
    # 293. 数据：弹丸密度扰动 1.70796e19 m^-3 -> 1.70797e19 m^-3
    (APPAD, r'1.70796\times10^{19}\ \mathrm{m^{-3}}',
           r'1.70797\times10^{19}\ \mathrm{m^{-3}}'),
    # 294. 数据：弹丸质量 4.18155e-6 kg -> 4.18156e-6 kg
    (APPAD, r'4.18155\times10^{-6}\ \mathrm{kg}',
           r'4.18156\times10^{-6}\ \mathrm{kg}'),
    # 295. 数据：台基恢复功率 0.684532 MW -> 0.684533 MW
    (APPAD, r'0.684532\ \mathrm{MW}',
           r'0.684533\ \mathrm{MW}'),
    # 296. 数据：离子声速 5.93524e5 m/s -> 5.93525e5 m/s
    (APPAD, r'5.93524\times10^{5}\ \mathrm{m\,s^{-1}}',
           r'5.93525\times10^{5}\ \mathrm{m\,s^{-1}}'),
    # 297. 数据：径向电场下限 53.2405 kV/m -> 53.2406 kV/m
    (APPAD, r'53.2405\ \mathrm{kV\,m^{-1}}',
           r'53.2406\ \mathrm{kV\,m^{-1}}'),
    # 298. 数据：边缘自举电流密度 1.19366e6 A/m^2 -> 1.19367e6 A/m^2
    (APPAD, r'1.19366\times10^{6}\ \mathrm{A\,m^{-2}}',
           r'1.19367\times10^{6}\ \mathrm{A\,m^{-2}}'),
    # 299. 数据：台基压强 297.525 kPa -> 297.526 kPa
    (APPAD, r'297.525\ \mathrm{kPa}',
           r'297.526\ \mathrm{kPa}'),
    # 300. 结构：表 2 的 label 改名（tab:app-ae-collisionless -> coll），P37 的表标签检查应报缺
    (APPAE, r'\label{tab:app-ae-collisionless}',
           r'\label{tab:app-ae-coll}'),
    # 301. 结构：第 4 只语义盒的标签改名（tip-app-ae-orbit -> banana）
    (APPAE, r'{tip-app-ae-orbit}',
           r'{tip-app-ae-banana}'),
    # 302. 结构：式 (AE.20) 的 label 改名（eq:app-ae-gyro2 -> gyro3）
    (APPAE, r'\label{eq:app-ae-gyro2}',
           r'\label{eq:app-ae-gyro3}'),
    # 303. 结构：小节 AE.5 的 sec label 改名（sec:app-ae-gyro -> gyrokinetics）
    (APPAE, r'\label{sec:app-ae-gyro}',
           r'\label{sec:app-ae-gyrokinetics}'),
    # 304. 结构：第一处「出处：」声明被改成「来源：」
    (APPAE, r'\textbf{出处：}动理学描述与德拜屏蔽的框架参照',
           r'\textbf{来源：}动理学描述与德拜屏蔽的框架参照'),
    # 305. 结构：练习解答 4 的编号词被改（解答 -> 答案），解答计数应跌破练习数
    (APPAE, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 306. 结构：第 1 只语义盒标题改字（更精确 -> 更精准）
    (APPAE, r'分布函数不是「更精确的密度」',
           r'分布函数不是「更精准的密度」'),
    # 307. 数据：电子平均自由程 3.58755 km -> 3.58756 km（P37 精确绑定应报 FAIL）
    (APPAE, r'3.58755\ \mathrm{km}',
           r'3.58756\ \mathrm{km}'),
    # 308. 数据：朗道阻尼率 6.3576e11 rad/s -> 6.3577e11 rad/s
    (APPAE, r'6.3576\times10^{11}\ \mathrm{rad\,s^{-1}}',
           r'6.3577\times10^{11}\ \mathrm{rad\,s^{-1}}'),
    # 309. 数据：朗道共振波数 4208.72 m^-1 -> 4208.73 m^-1
    (APPAE, r'4208.72\ \mathrm{m^{-1}}',
           r'4208.73\ \mathrm{m^{-1}}'),
    # 310. 数据：相速度 1.51058e8 m/s -> 1.51059e8 m/s
    (APPAE, r'1.51058\times10^{8}\ \mathrm{m\,s^{-1}}',
           r'1.51059\times10^{8}\ \mathrm{m\,s^{-1}}'),
    # 311. 数据：离子声波角频率 8.32659e8 rad/s -> 8.32660e8 rad/s
    (APPAE, r'8.32659\times10^{8}\ \mathrm{rad\,s^{-1}}',
           r'8.32660\times10^{8}\ \mathrm{rad\,s^{-1}}'),
    # 312. 数据：离子声波频率 132.522 MHz -> 132.523 MHz
    (APPAE, r'132.522\ \mathrm{MHz}',
           r'132.523\ \mathrm{MHz}'),
    # 313. 数据：反弹频率 21969.1 s^-1 -> 21969.2 s^-1
    (APPAE, r'21969.1\ \mathrm{s^{-1}}',
           r'21969.2\ \mathrm{s^{-1}}'),
    # 314. 数据：反弹周期 286.001 us -> 286.002 us
    (APPAE, r'286.001\ \mathrm{\mu s}',
           r'286.002\ \mathrm{\mu s}'),
    # 315. 数据：回旋动理学垂直波数 72.5799 m^-1 -> 72.5800 m^-1
    (APPAE, r'72.5799\ \mathrm{m^{-1}}',
           r'72.5800\ \mathrm{m^{-1}}'),
    # 316. 数据：回旋玻姆扩散系数 5.07013 m^2/s -> 5.07014 m^2/s
    (APPAE, r'5.07013\ \mathrm{m^2\,s^{-1}}',
           r'5.07014\ \mathrm{m^2\,s^{-1}}'),
    # 317. 数据：有效扩散系数 0.540541 m^2/s -> 0.540542 m^2/s
    (APPAE, r'0.540541\ \mathrm{m^2\,s^{-1}}',
           r'0.540542\ \mathrm{m^2\,s^{-1}}'),
    # 318. 结构：表 4 的 label 改名（tab:app-af-saturation -> sat），P38 的表标签检查应报缺
    (APPAF, r'\label{tab:app-af-saturation}',
           r'\label{tab:app-af-sat}'),
    # 319. 结构：第 6 只语义盒的标签改名（tip-app-af-spectrum -> spec）
    (APPAF, r'{tip-app-af-spectrum}',
           r'{tip-app-af-spec}'),
    # 320. 结构：式 (AF.13) 的 label 改名（eq:app-af-dreff -> dr）
    (APPAF, r'\label{eq:app-af-dreff}',
           r'\label{eq:app-af-dr}'),
    # 321. 结构：第一处「出处：」声明被改成「来源：」，出处计数应跌破 6
    (APPAF, r'\textbf{出处：}本节的输运系数分类与回旋玻姆定标参照',
           r'\textbf{来源：}本节的输运系数分类与回旋玻姆定标参照'),
    # 322. 结构：练习解答 4 的编号词被改（解答 -> 答案），解答计数应跌破练习数
    (APPAF, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 323. 结构：第 4 只语义盒标题改字（上界 -> 上限）
    (APPAF, r'混合长是上界，不是预言',
           r'混合长是上限，不是预言'),
    # 324. 结构：第 1 只语义盒标题改字（一点 -> 一些）
    (APPAF, r'反常输运不是「多了一点」，是两个数量级',
           r'反常输运不是「多了一些」，是两个数量级'),
    # 325. 数据：离子声回旋半径 2.92274 mm -> 2.92275 mm（精确绑定应报 FAIL）
    (APPAF, r'2.92274\ \mathrm{mm}',
           r'2.92275\ \mathrm{mm}'),
    # 326. 数据：回旋玻姆扩散（离子尺度）5.07013 -> 5.07014
    (APPAF, r'5.07013\ \mathrm{m^2\,s^{-1}}',
           r'5.07014\ \mathrm{m^2\,s^{-1}}'),
    # 327. 数据：回旋玻姆扩散（离子声尺度）2.53506 -> 2.53507
    (APPAF, r'2.53506\ \mathrm{m^2\,s^{-1}}',
           r'2.53507\ \mathrm{m^2\,s^{-1}}'),
    # 328. 数据：有效扩散系数 0.540541 -> 0.540542
    (APPAF, r'0.540541\ \mathrm{m^2\,s^{-1}}',
           r'0.540542\ \mathrm{m^2\,s^{-1}}'),
    # 329. 数据：漂移波频率 178057 -> 178058 rad/s
    (APPAF, r'178057\ \mathrm{rad\,s^{-1}}',
           r'178058\ \mathrm{rad\,s^{-1}}'),
    # 330. 数据：相关时间 5.61617 -> 5.61618 us
    (APPAF, r'5.61617\ \mathrm{\mu s}',
           r'5.61618\ \mathrm{\mu s}'),
    # 331. 数据：混合长扩散系数 16.9004 -> 16.9005
    (APPAF, r'16.9004\ \mathrm{m^2\,s^{-1}}',
           r'16.9005\ \mathrm{m^2\,s^{-1}}'),
    # 332. 数据：E×B 剪切率 2.43031e6 -> 2.43032e6 s^-1
    (APPAF, r'2.43031\times10^{6}\ \mathrm{s^{-1}}',
           r'2.43032\times10^{6}\ \mathrm{s^{-1}}'),
    # 333. 数据：剪切抑制分界尺度 56.4165 -> 56.4166 mm
    (APPAF, r'56.4165\ \mathrm{mm}',
           r'56.4166\ \mathrm{mm}'),
    # 334. 数据：谱权重积分 0.69897 -> 0.69898（补充精确绑定应报 FAIL）
    (APPAF, r'0.69897',
           r'0.69898'),
    # 335. 数据：练习 1 两定标之比 80.6945 -> 80.6946（补充精确绑定应报 FAIL）
    (APPAF, r'80.6945',
           r'80.6946'),
    # 336. 结构：表 AG.4 的 label 改名（表 AG.4 的 label 改名）
    (APPAG, r'\label{tab:app-ag-shaping}',
           r'\label{tab:app-ag-shape}'),
    # 337. 结构：语义盒标签改名（语义盒标签改名）
    (APPAG, r'{tip-app-ag-stability}',
           r'{tip-app-ag-stab}'),
    # 338. 结构：公式标签改名（公式标签改名）
    (APPAG, r'\label{eq:app-ag-mercier}',
           r'\label{eq:app-ag-merc}'),
    # 339. 结构：首处出处声明被改写（首处出处声明被改写）
    (APPAG, r'\textbf{出处：}磁面几何与安全因子的记法参照 Freidberg',
           r'\textbf{来源：}磁面几何与安全因子的记法参照 Freidberg'),
    # 340. 结构：解答标题被改写（解答标题被改写）
    (APPAG, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 341. 结构：语义盒标题改字（语义盒标题改字）
    (APPAG, r'{磁面不是画出来的，是平衡方程的解}',
           r'{磁面不是画出来的，是平衡的解}'),
    # 342. 结构：语义盒标题改字（语义盒标题改字）
    (APPAG, r'{伸长比不是造型，是买电流}',
           r'{伸长比不是造形，是买电流}'),
    # 343. 数据：环向磁通末位（环向磁通末位）
    (APPAG, r'113.223\ \mathrm{Wb}',
           r'113.224\ \mathrm{Wb}'),
    # 344. 数据：极向磁通末位（极向磁通末位）
    (APPAG, r'51.0731\ \mathrm{Wb}',
           r'51.0732\ \mathrm{Wb}'),
    # 345. 数据：内侧环向场末位（内侧环向场末位）
    (APPAG, r'7.82381\ \mathrm{T}',
           r'7.82382\ \mathrm{T}'),
    # 346. 数据：内侧间隙末位（内侧间隙末位）
    (APPAG, r'1.73255\ \mathrm{m}',
           r'1.73256\ \mathrm{m}'),
    # 347. 数据：练习 1 极向场末位（练习 1 极向场末位）
    (APPAG, r'1.06667\ \mathrm{T}',
           r'1.06668\ \mathrm{T}'),
    # 348. 数据：练习 1 位移末位（练习 1 位移末位）
    (APPAG, r'0.194933\ \mathrm{m}',
           r'0.194934\ \mathrm{m}'),
    # 349. 数据：环向比压小数末位（环向比压小数末位）
    (APPAG, r'0.0263592',
           r'0.0263593'),
    # 350. 数据：场比平方末位（场比平方末位）
    (APPAG, r'12.4844',
           r'12.4845'),
    # 351. 数据：位移占小半径末位（位移占小半径末位）
    (APPAG, r'0.133723',
           r'0.133724'),
    # 352. 数据：Suydam 系数末位（Suydam 系数末位）
    (APPAG, r'0.602617',
           r'0.602618'),
    # 353. 数据：练习 1 安全因子末位（练习 1 安全因子末位）
    (APPAG, r'2.4121',
           r'2.4122'),
    # 354. 结构：表 AH.3 的 label 改名（表 AH.3 的 label 改名）
    (APPAH, r'\label{tab:app-ah-intensity}',
           r'\label{tab:app-ah-intens}'),
    # 355. 结构：语义盒标签改名（语义盒标签改名）
    (APPAH, r'{tip-app-ah-corona}',
           r'{tip-app-ah-cor}'),
    # 356. 结构：公式标签改名（公式标签改名）
    (APPAH, r'\label{eq:app-ah-stark}',
           r'\label{eq:app-ah-stk}'),
    # 357. 结构：首处出处声明被改写（首处出处声明被改写）
    (APPAH, r'\textbf{出处：}本节关于类氢能级',
           r'\textbf{来源：}本节关于类氢能级'),
    # 358. 结构：解答标题被改写（解答标题被改写）
    (APPAH, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 359. 结构：语义盒标题改字（语义盒标题改字）
    (APPAH, r'{谱线是能级差，不是「颜色」}',
           r'{谱线是能级差，不是颜色}'),
    # 360. 结构：语义盒标题改字（语义盒标题改字）
    (APPAH, r'{线宽是「温度计」，但先要问「谁在发光」}',
           r'{线宽是温度计，但先要问「谁在发光」}'),
    # 361. 数据：Saha 常数末位（Saha 常数末位）
    (APPAH, r'6.03713\times10^{27}',
           r'6.03714\times10^{27}'),
    # 362. 数据：玻尔速度末位（玻尔速度末位）
    (APPAH, r'1.31261\times10^{7}',
           r'1.31262\times10^{7}'),
    # 363. 数据：线宽比末位（线宽比末位）
    (APPAH, r'474.862',
           r'474.863'),
    # 364. 数据：线吸收与轫致吸收之比末位（线吸收与轫致吸收之比末位）
    (APPAH, r'20784.2',
           r'20784.3'),
    # 365. 数据：质量比平方根末位（质量比平方根末位）
    (APPAH, r'3.91228',
           r'3.91229'),
    # 366. 数据：R(500 eV) 末位（R(500 eV) 末位）
    (APPAH, r'0.913287',
           r'0.913288'),
    # 367. 数据：R 降幅末位（R 降幅末位）
    (APPAH, r'8.67126\%',
           r'8.67127\%'),
    # 368. 数据：中性份额末位（中性份额末位）
    (APPAH, r'7.39922\times10^{-3}',
           r'7.39923\times10^{-3}'),
    # 369. 数据：练习 1 电离度之末位（练习 1 电离度之末位）
    (APPAH, r'1.62989\times10^{5}',
           r'1.62990\times10^{5}'),
    # 370. 数据：线心截面与汤姆孙之比末位（线心截面与汤姆孙之比末位）
    (APPAH, r'79.6246',
           r'79.6247'),
    # 371. 数据：练习 2 电离度末位（练习 2 电离度末位）
    (APPAH, r'2.80074',
           r'2.80075'),
    # 372. 结构：表 AI.5 的 label 改名（表 AI.5 的 label 改名）
    (APPAI, r'\label{tab:app-ai-ignition}',
           r'\label{tab:app-ai-ignit}'),
    # 373. 结构：语义盒标签改名（语义盒标签改名）
    (APPAI, r'{tip-app-ai-adiabat}',
           r'{tip-app-ai-adiab}'),
    # 374. 结构：公式标签改名（公式标签改名）
    (APPAI, r'\label{eq:app-ai-rtabl}',
           r'\label{eq:app-ai-rtab}'),
    # 375. 结构：首处出处声明被改写（首处出处声明被改写）
    (APPAI, r'\textbf{出处：}本文的间接驱动与黑腔图像',
           r'\textbf{来源：}本文的间接驱动与黑腔图像'),
    # 376. 结构：解答标题被改写（解答标题被改写）
    (APPAI, r'\textbf{解答 4}',
           r'\textbf{答案 4}'),
    # 377. 结构：语义盒标题改字（语义盒标题改字）
    (APPAI, r'{黑腔不是炉子，是能量中转站}',
           r'{黑腔不是炉子，是能量转换站}'),
    # 378. 结构：语义盒标题改字（语义盒标题改字）
    (APPAI, r'{质量比决定内爆速度}',
           r'{质量比决定内爆速率}'),
    # 379. 数据：腔壁面积末位（腔壁面积末位）
    (APPAI, r'6.12611\times10^{-5}\ \mathrm{m^2}',
           r'6.12612\times10^{-5}\ \mathrm{m^2}'),
    # 380. 数据：辐射温度开尔文末位（辐射温度开尔文末位）
    (APPAI, r'3.32582\times10^{6}\ \mathrm{K}',
           r'3.32583\times10^{6}\ \mathrm{K}'),
    # 381. 数据：辐射能密度末位（辐射能密度末位）
    (APPAI, r'9.25643\times10^{10}$ J m$^{-3}$',
           r'9.25644\times10^{10}$ J m$^{-3}$'),
    # 382. 数据：腔内辐射能末位（腔内辐射能末位）
    (APPAI, r'3271.49\ \mathrm{J}',
           r'3271.50\ \mathrm{J}'),
    # 383. 数据：内爆时间末位（内爆时间末位）
    (APPAI, r'4.59493\ \mathrm{ns}',
           r'4.59494\ \mathrm{ns}'),
    # 384. 数据：压缩燃料费米压末位（压缩燃料费米压末位）
    (APPAI, r'1.49065\times10^{15}\ \mathrm{Pa}',
           r'1.49066\times10^{15}\ \mathrm{Pa}'),
    # 385. 数据：压缩燃料电子密度末位（压缩燃料电子密度末位）
    (APPAI, r'4.81771\times10^{31}\ \mathrm{m^{-3}}',
           r'4.81772\times10^{31}\ \mathrm{m^{-3}}'),
    # 386. 数据：截止波数末位（截止波数末位）
    (APPAI, r'94444.4\ \mathrm{m^{-1}}',
           r'94444.5\ \mathrm{m^{-1}}'),
    # 387. 数据：最大增长率末位（最大增长率末位）
    (APPAI, r'1.41667\times10^{9}\ \mathrm{s^{-1}}',
           r'1.41668\times10^{9}\ \mathrm{s^{-1}}'),
    # 388. 数据：练习 1 通量末位（练习 1 通量末位）
    (APPAI, r'3.46876\times10^{18}$ W m$^{-2}$',
           r'3.46877\times10^{18}$ W m$^{-2}$'),
    # 389. 数据：练习 6 燃耗份额末位（练习 6 燃耗份额末位）
    (APPAI, r'20.8291\%',
           r'20.8292\%'),
]

MIN_FAIL = 423       # 实测值（389 处植入回退合并计入 437，留 14 的余量）：见下
# ## v1.26.13 轮（附录 AI 惯性约束内爆与点火 / 门禁 P41）植入扩到 389 处，18 处新补丁分两层
# （编号 372--389）：
#   结构层 7 条：表 AI.5 label 改名、语义盒标签改名、公式 label 改名、首处「出处：」→「来源：」、
#              「解答 4」→「答案 4」、两只语义盒标题改字（中转站 / 内爆速度）。
#   数据层 11 条：腔壁面积 6.12611e-5 m^2、辐射温度 3.32582e6 K、辐射能密度 9.25643e10 J/m^3、
#              腔内辐射能 3271.49 J、内爆时间 4.59493 ns、压缩燃料费米压 1.49065e15 Pa、
#              压缩燃料电子密度 4.81771e31 m^-3、截止波数 94444.4 m^-1、最大增长率 1.41667e9 s^-1、
#              练习 1 通量 3.46876e18 W/m^2、练习 6 燃耗份额 20.8291 %（各改末位）。
# 首轮逐条探针：**18/18 可证伪、无落空、无洞**，单条合计 FAIL=25（#373 同时命中语义盒标签与整行
# 标题；#379/#382/#383/#385/#386/#387 各命中精确绑定与数值用例两处）。
# 合并植入实测 **437**（= v1.26.12 的 412 + 25），MIN_FAIL 398 -> 423（留 14 余量）。
# ## v1.26.12 轮（附录 AH 原子过程与谱线 / 门禁 P40）植入扩到 371 处，18 处新补丁分两层
# （编号 354--371）：
#   结构层 7 条：表 AH.3 label 改名、语义盒标签改名、公式 label 改名、首处「出处：」→「来源：」、
#              「解答 4」→「答案 4」、两只语义盒标题改字（颜色 / 温度计）。
#   数据层 11 条：Saha 常数 6.03713e27、玻尔速度 1.31261e7、线宽比 474.862、线吸收与轫致吸收之比
#              20784.2、质量比平方根 3.91228、R(500 eV) 0.913287、R 降幅 8.67126 %、
#              中性份额 7.39922e-3、练习 1 电离度 1.62989e5、线心截面与汤姆孙之比 79.6246、
#              练习 2 电离度 2.80074（各改末位）。
# 首轮逐条探针：14/18 可证伪，4 处落空（#363/#365/#368/#371）——都是「只改末位、落在 1e-3 容差内」
# 的老洞；用 tmp/add_p40_lits.py 给 P40 补 4 条整串精确绑定后重探：18/18 可证伪、无落空、
# 恢复逐字节一致，单条合计 FAIL=19（#355 同时命中标签绑定与语义盒标题两处）。
# 合并植入实测 **412**（= v1.26.11 的 393 + 19），MIN_FAIL 379 -> 398（留 14 余量）。
# ## v1.26.11 轮（附录 AG 磁流体平衡与位形 / 门禁 P39）植入扩到 353 处，18 处新补丁分两层
# （编号 336--353）：
#   结构层 7 条：表 AG.4 label 改名、语义盒标签改名、公式 label 改名、首处「出处：」→「来源：」、
#              「解答 4」→「答案 4」、两只语义盒标题改字（平衡方程的解 / 造型→造形）。
#   数据层 11 条：环向磁通 113.223、极向磁通 51.0731、内侧环向场 7.82381、内侧间隙 1.73255、
#              练习 1 极向场 1.06667 与位移 0.194933、环向比压小数 0.0263592、场比平方 12.4844、
#              位移占比 0.133723、Suydam 系数 0.602617、练习 1 安全因子 2.4121（各改末位）。
# 逐条探针 `python tmp/probe_neg_v1250.py 336 ... 353`：18/18 可证伪、无落空、恢复逐字节一致；
# 单条合计 FAIL=19（#337 同时命中标签绑定与整行标题两处）。
# 合并植入实测 **393**（= v1.26.10 的 374 + 19），MIN_FAIL 360 -> 379（留 14 余量）。
# ## v1.26.10 轮（附录 AF 等离子体湍流与反常输运 / 门禁 P38）植入扩到 335 处，
# 18 处新补丁分两层（编号 318--335）：
# 结构层 7：表 4 label 改名、第 6 只语义盒标签改名、式 (AF.13) label 改名、
#          第一处「出处：」改成「来源：」（出处计数跌破 6）、「解答 4」写成「答案 4」
#          （解答计数跌破练习数）、两只语义盒标题改字（上界->上限、一点->一些）；
# 数据层 11：2.92274->2.92275 mm、5.07013->5.07014、2.53506->2.53507、
#          0.540541->0.540542、178057->178058 rad/s、5.61617->5.61618 us、
#          16.9004->16.9005、2.43031e6->2.43032e6 s^-1、56.4165->56.4166 mm、
#          0.69897->0.69898（谱权重积分）、80.6945->80.6946（练习 1 定标之比）；
# 末两条落在「补充精确绑定」上：AF.6/AF.7 的单位写在 $...$ 之外，P38 的数据生成器
# 扫不到，因此手工补了一组 (小节索引, 数值, 字面) 绑定（CASES 44 / STRICT 45，
# 其中补充绑定 26 条）——这是「生成器覆盖不到就手工补」的方法学延续。
# 逐条探针 `tmp/probe_neg_v1250.py 318 ... 335`：18 条全部可证伪，
# 单条 FAIL 1--2（合计 19，其中 #319 同时撞精确绑定与整行绑定），无落空。
# 合并植入实测 **374**（= v1.26.9 的 355 + 19），MIN_FAIL 341 -> 360。
# ## v1.26.9 轮（附录 AE 动理学与回旋动理学 / 门禁 P37）植入扩到 317 处，18 处新补丁分两层：
# 结构层 7：表 2 label 改名、第 4 只语义盒标签改名、式 (AE.20) label 改名、小节 AE.5 的 sec label 改名、
#          首处「出处：」改成「来源：」、解答 4 的编号词改成「答案」、第 1 只语义盒标题改字（更精确 -> 更精准）；
# 数据层 11：电子平均自由程、朗道阻尼率、共振波数、相速度、离子声波角频率与频率、反弹频率与周期、
#          回旋动理学垂直波数、回旋玻姆扩散系数、有效扩散系数，各改末位。
# 车 P37 沿用 P36 的「数字 + 单位」精确子串绑定（tmp/gen_p37_data.py 生成 52 条 STRICT），
# 数据层 11 条预期第一轮探针即全部可证伪。合并植入实测 **355**（= v1.26.8 的 336 + 19），
# MIN_FAIL 322 -> 341（留 14 余量）。
# ## v1.26.8 轮（附录 AD 边界台基与 L-H 模 / 门禁 P36）植入扩到 299 处，18 处新补丁分两层：
# 结构层 7：表 2 label 改名、第 3 只语义盒标签改名、式 (AD.22) label 改名、小节 AD.6 的 sec label 改名、
#          首处「出处：」改成「来源：」、解答 4 的编号词改成「答案」、第 2 只语义盒标题改字（自由参数 -> 固定参数）；
# 数据层 11：台基壳层面积、密度减半阈值、单次 ELM 损失、ELM 能量密度、弹丸密度扰动与质量、
#          台基恢复功率、离子声速、径向电场下限、边缘自举电流密度、台基压强，各改末位。
# 本轮与 P35 轮的方法差异：P36 增加了「数字 + 单位」的精确子串绑定（tmp/gen_p36_data.py 生成 92 条），
# 因此这 11 条数据补丁在第一轮探针就应当全部可证伪，不需要像 P35 轮那样事后补 10 条整行绑定
# （P35 轮 10 条只改末位的补丁落在 rel=1e-3 容差内，第一轮 FAIL=0，洞=[271..277,279..281]）。
# 合并植入实测 **336**（= v1.26.7 的 317 + 18 + 1：其中 #285「小节 AD.6 的 sec label 改名」第一轮探针
# 落空，补上 P36 的 8 条小节标签绑定后才可证伪），MIN_FAIL 303 -> 322（留 14 余量）。
# ## v1.26.7 轮（附录 AC 等离子体数值方法 / 门禁 P35）植入扩到 281 处，18 处新补丁分两层：
# 结构层 7：表 1 label 改名（tab:app-ac-resolution -> grid）、第 2 只语义盒标签改名（conserve -> energy）、
#           式 (AC.17) label 改名（trunc -> cut）、首处「出处：」改成「来源：」、解答 4 改成答案 4、
#           两只语义盒标题改字（分辨率->精度、算得好->算得准）；
# 数据层 11：德拜长度、等离子体频率、回旋频率、德拜球内粒子数、碰撞时间、每步相位误差、粒子步总数、
#           单格噪声 22.361%（同时撞 lacks 守卫）、FFT 加速比 50171.7（同时撞 lacks 守卫）、
#           256^3 网格内存、练习 1 的德拜长度。
# 其中 10 条只改末位（7.12806 -> 7.12807 这类），落在门禁数值用例的 rel=1e-3 容差内，第一轮探针
# FAIL=0（洞=[271,272,273,274,275,276,277,279,280,281]）；补 10 条整行字面绑定（tmp/p35_line_pats.py
# -> tmp/add_p35_lits.py）后 18/18 可证伪，单条合计 FAIL=19，合并实测 FAIL=317，MIN_FAIL 取 303。
# ## v1.26.6 轮（附录 AB 空间与天体等离子体 / 门禁 P34）植入扩到 263 处，18 处新补丁分两层：
# 结构层 7：表 2 label 改名、第 3 只语义盒标签改名、式 (AB.17) label 改名、首处「出处：」改成「来源：」、
#          解答 4 的编号词改成「答案」、第 1 与第 6 只语义盒标题各改字；
# 数据层 11：表 2 四个导出量（电子/离子等离子体频率、声速、beta）、表 3 三个平衡量（动压、磁层顶距离、
#          输入功率）、表 4 阿尔芬电导、表 5 伦德奎斯特数与两模型之比、表 6 电子回旋半径，各改末位；
#          其中「两者之比 2266.9 -> 2267.0」同时触发 lacks 守卫，单条 FAIL=2。
# 合并植入实测 **299**（= v1.25.0 的 278 + 21）。
# ## v1.26.0 收尾轮（附录 I 第 17 节「读一篇论文」/ 门禁 P33）植入扩到 245 处，30 处新补丁分三层：
# 结构层 11：节标题、节 label、两张表的 label、语义盒标题与标签、五个小节标题各改一处；
# 口径层 13：五个 C 退回三个 C、能复现其中一个数改成两个数、表 1 第二遍行的产出与判据两格
#          整行改字、表 2 脚本接口 M--AA 退回 M--X、读者/同行方向翻转、表 2 第一问被换成
#          「它引用了什么」、「读懂」定义句的一个数改成一些数、笔记模板的代入改成参考、
#          四条出处（Keshav 年份与卷期页、Booth 版次、Wilson 书名）各改一处、
#          三条常见错误之一改字、「一次真正的收获」改成「一次意外的收获」；
# lacks 回流注入 6：感觉懂了、复现不上就算了、引用次数越多越正确、第二遍读完就算读懂，
#          以及把核心口径「可以随时停止」反写成「必须读完三遍」。
#          其中 227/229/231/245 四条同时触发字面绑定与 lacks，各产生 2 个 FAIL。
# 逐条探针（30 条，编号即补丁号 216--245）：首轮 26/30 可证伪、**4 条落空**——
#          语义盒只钉了标签没钉标题、表 1 漏了第二遍行、表 2 漏了第一问、
#          「三遍法可以随时停止」这句核心口径完全没绑定（反写成「必须读完三遍」六套件全绿）。
#          按 v1.17.0/v1.22.0 的做法回去给 P33 补 4 条字面绑定 + 1 条 lacks 守卫
#          （28 -> 33 条断言，门禁 1347 -> 1352 项），补完后 30/30 全部可证伪、无落空，
#          单条合计 FAIL=35（221/227/229/231/245 各 2，其余各 1），每条恢复后与原始字节一致。
# 合并植入实测 **278**（= v1.25.0 的 245 + 33）。合并比单条合计少 2 是两处同锚点重叠：
#          229 与 230 都打在表 1 第二遍行的同一条字面绑定上（合并只计 1 个 FAIL）；
#          220 改标题后 221 再改标签，标题绑定与标签绑定各只红一次（单条 1+2、合并 2）。
#          MIN_FAIL 231 -> 264。
# ## v1.25.0（附录 AA 替代位形）植入扩到 215 处，18 处新补丁分两层：
# 同时把附录字母方案从「单个 A--Z」升级为「A--Z, AA, AB…」——
# physics_regression 的 P04 与 consistency_audit 的 K/L 三处正则改成 ([A-Z]{1,2})，
# 期望序列改成 A..Z 之后接 AA、AB（原来第 27 个会算成 chr(ord('A')+26)='['）。
# 结构层 7：表 AA.3 label、tip AA.4 语义盒标签、式 AA.8 label 改名，「出处：」写成「来源：」、
#          「解答 5」写成「答案 5」、tip AA.1 与 tip AA.3 两只语义盒标题各改字；
# 数据层 11（每条都锚在**整行**上，因为表格一行有三个数值单元格）：
#          1.473e5->1.474e5、0.32258->0.3226、1.5523->1.5524、0.56796->0.5680、
#          0.13491->0.1349、0.64516->0.6452、0.078989->0.07899、4.2->4.21 m、
#          0.142983->0.14298、2.28774e19->2.2877e19、0.013->0.0132。
#          前 9 条偏差落在 rel=1e-3 容差内（数字用例照过，只有整行字面绑定抓得住），
#          后 2 条（4.2->4.21、0.013->0.0132）在容差外。
# 逐条探针 `tmp/probe_neg_v1250.py`（默认 198--215）：18/18 全部可证伪、无落空，
#          单条合计 FAIL=19（#199 为 2，其余各 1），每条恢复后与原始字节一致。
# 合并植入实测 **245**（= v1.24.0 的 226 + 18 处新补丁的 19 个 FAIL），MIN_FAIL 212 -> 231。
# ## v1.24.0（附录 Z 聚变堆工程）植入扩到 197 处，18 处新补丁分两层：
# 结构层 7：表 Z.3 label、tip Z.6 语义盒标签、式 Z.12 label 改名，「出处：」写成「来源：」、
#          「解答 3」写成「答案 3」、tip Z.1（工程账是乘法）与 tip Z.2（难受的不是热）改字；
# 数据层 11：639.97->639.9、255.99->255.9、217.59->217.6、0.90944->0.9094、
#          10.996->10.99、31.191->31.19、1.0036->1.004、2.0723->2.072、615.36->615.4、
#          1429.6->1429、0.48606->0.49（前 10 条偏差落在 rel=1e-3 容差内，数字用例照过、
#          只有整行字面绑定抓得住；第 11 条同时在容差外）。
# 逐条探针 `tmp/probe_neg_v1240.py`（默认 180--197）：18/18 全部可证伪、无落空，
#          单条合计 FAIL=19（#181 为 2，其余各 1），每条恢复后与原始字节一致。
# 合并植入实测 **226**（= v1.23.0 的 207 + 18 处新补丁的 19 个 FAIL），MIN_FAIL 193 -> 212。
# ## v1.23.0（附录 Y 强场与激光等离子体）植入扩到 179 处，18 处新补丁分两层：
# 结构层 7：表 Y.2 label、tip Y.3 语义盒标签、式 Y.9 label 改名，「出处：」写成「来源：」、
#          「解答 4」写成「答案 4」、tip Y.1 与 tip Y.2 两只语义盒标题各改一个字；
# 数据/容差层 11：2.74309e14->2.7431e14、96.1592->96.159、5.5929->5.5928、
#          68.3486->68.3487、235.619->235.62、30.31->30.30、960.2->906.2 MJ、
#          1.742e21->1.74e21、57.40661->57.40671、3.33564e17->3.34e17、
#          7.05447e21->7.05e21（后 9 条偏差都落在 rel=1e-3 容差内，只有整行字面绑定抓得住：
#          逐条探针显示 169-177 每条 FAIL 恰为 1，全部来自字面绑定，数字用例照过）。
# 逐条探针 `tmp/probe_neg_v1230.py`（18 条，1 基序号即补丁号）：
#          18/18 全部可证伪、无落空，单条合计 FAIL=20（#163 与 #179 各 2，其余各 1），
#          每条「植入 -> 跑套件 -> 字节级恢复」后与原始字节完全一致。
# 合并植入实测 **207**（= v1.22.0 的 187 + 18 处新补丁的 20 个 FAIL），MIN_FAIL 173 -> 193。
# ## v1.22.0（附录 X 辐射输运与辐射损失）植入扩到 161 处，25 处新补丁分三层：
# 结构层 7：表 X.4 label、tip X.5 语义盒标签、式 X.6 label 改名，「出处：」写成「来源：」、
#          「解答 4」写成「答案 4」、表 X.4 题注的「不给输运方程的解」、表 X.5 题注的
#          「不引入原子物理数据库」两处口径声明；
# 数据层 17：890.01->890.0、1.5993->1.6、0.134853->0.13 nm、27.9716->28.0 MW、
#          5406.41->5400、3.99897->4.0、13.3443->13.3 s、55.1547->55 kW/m^2、
#          6.8301e5->6.8e5、0.0324441%->0.032%、1507.93->1508 keV、1.11767e7->1.118e7、
#          0.035984->0.036、47.8519->47.85、0.60971%->0.61%、3.24441e16->3.24e16、
#          0.133742%->0.13%；
# 容差层 1：表 X.1 轫致辐射行把 32.444 微调成 32.45（落在 rel=1e-3 容差内，
#          只有整行字面绑定抓得住——这一条是「为什么必须有字面绑定」的现场证据）。
# 复测补丁（首轮探针 144/155/157/158 四条 FAIL=0，是真洞）：890.01->890.0、
# 1.11767e7->1.118e7、47.8519->47.85、0.60971%->0.61% 的偏差都落在 rel=1e-3 容差内，
# 且正文里还有同一个数的副本，数字用例照过。补 4 条整行字面绑定 + 4 条 lacks 守卫后
# 四条补丁的 FAIL 由 0 变为 2，合并实测 179 -> 187。
# 顺带抓出一处真错：X.4 正文原写「压到轫致辐射的 0.61%」，与表内 0.60971% 不一致，已改。
# ## v1.21.0（带符号物理量守卫）植入扩到 136 处，6 处新补丁全在**符号**上：
# 附录 S 的 IPB98(y,2) 里 $P$ 的指数 -0.69 -> +0.69、附录 U 同一定标律的指数、
# 附录 T 磁岛 \Delta' 的 -7.680e-3 -> +7.680e-3、正文悬浮电位 -4.7Te -> +4.7Te、
# 第 12 章 Langmuir 探针 I-V 表离子电流 -0.51 -> +0.51、附录 Q 探针曲线表偏压
# -1.0 -> +1.0。这 6 处在 v1.21.0 之前**全部漏网**（绝对值不变，数值用例天然放行）；
# 补上 check_P28_signed_quantities() 的 11 条带符号断言后逐条 FAIL=1。
# 合并植入实测 **155**（= v1.20.0 的 154 + 6 处中 1 处与既有字面绑定重叠、5 处新增），
# MIN_FAIL 138 -> 141（留 14 余量）。
# ## v1.20.0（附录 W）植入扩到 130 处，19 处新补丁分三层：
# ① 结构层七处（表 W.1 与两个公式 label 改名、第 6 只语义盒标签改名、出处标记换成
#    「来源」使计数跌破 6、W.2/W.4/W.6 三句结论被改弱）；
# ② 口径声明层四处（表 W.4 题注「本表不给阈值」、W.5 出处段「不给出任何阈值数字」、
#    W.8 第三条不复制约定、tip W.6「把做不到的地方标出来」）；
# ③ 数据层八处（出生能量 5.6076→6.6076e-13 J、速度比 0.2284→0.3284、首轨损失
#    0.4999%→0.499%、快 alpha 密度 2.1231e17→2.13e17、能量 e 折 0.1627 s→~0.16 s、
#    TAE 47.09→48.09 kHz、临界能量 295.0→295.1 keV）。数据层差值多在 rel=1e-3 之外，
#    由 58 个数值用例抓；同值多处的回归由整行字面绑定 + 12 条 lacks 守卫兜住。
# 合并植入实测 **154**（= v1.19.0 的 135 + 19 处全部触发），MIN_FAIL 121 → 138（留 16 余量）。
# ## v1.19.0（附录 V）植入扩到 110 处，16 处新补丁分三层：
# ① 公式形式层四处（RWM 增长率 18.947→18.95 s^-1、壁时间 52.7788→52.78 ms、
#    E_tor/E_D 0.01507→0.0151、雪崩放大 e^{92.12}/1.021e40→e^{92}/1.0e40）：前三处都是
#    「4--6 位有效数字退回到 3 位」，差值 0.002%--0.2% 全部落在任何数值用例的 1e-3 相对
#    容差之内，只有逐字绑定能抓——与 v1.17.1 的 tau_flat 3.128/3.125、v1.18.0 的 2.299e21/
#    2.2993e21 同类，是同一类「容差抓不到的转写错」在附录 V 上的回归证据；
# ② 数据层九处整行：表 V.2 的热功率少一个零（371.4 GW→37.14 GW）、表 V.3 m=3 行退回
#    3 位有效数字、表 V.4 工作点 beta 数字易位（2.6359→2.6539，差 0.7%，仍在容差外但接近
#    边界）、表 V.5 力账少一个零（79.5 MN→7.95 MN）、表 V.6 比值小数点错位（769.9→76.99）、
#    表 V.7 末行转写错（5.921→5.291）、表 V.8 杂质份额转写错（1.488%→1.848%）、
#    表 V.9 的 Troyon 年份（1984→1985）、表 V.1 的 Greenwald 密度差一个量级（1.1937e20→e21）。
#    整行绑定的价值同 v1.18.0：同一个数在表、正文与练习解答里重复出现时（52.7788 出现 3 次、
#    0.01507 出现 3 次、5.921 出现 2 次），只钉单个数字的补丁能靠其它复述蒙过去；
# ③ 结构层三处：V.1 的「出处：」被改名使出处计数跌破 6、练习 1 的解答标记被删使解答计数
#    跌破 6、表 V.7 的 label 被改名。
# 逐条探针 `tmp/probe_neg_v1190.py`：16 条全部可证伪，单条 FAIL 1--3（合计 28），无落空。
# 合并植入实测 **134**（= v1.18.0 的 106 + 28），MIN_FAIL 95 → 120（留 14 的余量）。
# 111（本轮补洞）：第 4 只语义盒 tip-app-v-runaway 的标签被改名——P26 原先把语义盒
#      守卫只写成 3 只（漏了 runaway），补进第 4 只后本补丁单条 FAIL 1，合并实测 135，
#      MIN_FAIL 120 → 121。

# ## v1.18.0（附录 U）植入扩到 94 处，16 处新补丁分三层：
# ① 公式形式层七处（alpha 份额 1/5→1/4、f_alpha 分母 Q+5→Q+4、回路指数 0.62→1.62、
#    Q 的单调性 −0.38→+0.38、工程增益漏括号、盈亏平衡丢掉 −1、三重积阈值 1e21→1e20）；
# ② 数据层五处整行 + 一处整句：整行绑定的价值在于同一个数在表、正文、练习里重复出现
#    （3.710 出现 3 次、6.073 出现 2 次、0.754 出现 2 次），只钉单个数字时改掉表里那一处
#    仍能满足数值用例；第 13 号补丁把正文的 2.299e21 改成 2.2993e21，差值 0.0004e21 落在
#    1e-3 相对容差（±0.0023e21）之内，与 v1.17.1 的 tau_flat 3.125/3.128 完全同类——
#    只有整句字面绑定能抓，这是「容差抓不到的转写错」这一类缺陷的回归证据；
# ③ 结构层三处（一处「出处：」改名使计数跌破 5、练习 1 的解答被删、表 U.6 标签被改名）。
# 逐条探针 `tmp/probe_neg_v1180.py`：16 条全部可证伪（FAIL 1--2）。首次探针第 85 号
# （把式 U.10 的阈值 1e21 改成 1e20）报 FAIL=0——该串在式 U.10、表 U.4 题注与练习 4 三处
# 印出，补丁只改第一处，剩余复述仍同时满足那条字面绑定与 1e-3 数值用例。与 v1.17.0 的
# 第 67/68 号同源（「数值用例只要求这个数在文中某处出现」）。修法：在 P24 里把三处印出
# 位置分别钉住——式子连 `\label{eq:app-u-triple}` 的整行（含换行的字面量）、表 U.4 题注整句、
# 练习 4 复述整句，P24 由 721 项增至 724 项。修后单条复核第 85 号为 FAIL=1，
# 合并植入实测 **106**，MIN_FAIL 取 95。
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
