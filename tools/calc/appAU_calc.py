# -*- coding: utf-8 -*-
"""附录 AU《测量、不确定度与数据分析》数值唯一来源（只 print，不写文件）。

口径：
  - 工作点沿用全书 D1：n_e0 = 1e20 m^-3、T_e = 9.194 keV、B = 5.3 T、a = 2.0 m、R = 6.2 m；
  - 跨附录对照取附录 AQ 的相位分辨 sigma_phi = 0.000361482 rad 与折合线积分密度 1.07979e15 m^-2；
  - 统计口径：GUM（JCGM 100:2008）A/B 类合成 + k=2 扩展；小样本用 Student t；
  - 例：汤姆孙散射 200 点、单点统计误差 5%（相对）、系统误差 2%（均匀分布取例）。

运行：python tmp\\appAU_calc.py > tmp\\appAU_calc.txt（PowerShell 用 Out-File -Encoding utf8）
"""
from __future__ import print_function

import math

# ---------------- 声明输入 ----------------
NE0 = 1.0e20          # m^-3
TE_EV = 9194.0        # eV
B0 = 5.3              # T
A_MINOR = 2.0         # m
R0 = 6.2              # m
N_PH = 200            # 每个空间点的光子数（统计量的来源）
REL_STAT = 0.05       # 单点（每点一次测量）统计误差（相对）
REL_SYS = 0.02        # 系统误差（相对，均匀分布取例；B 类）
N_REP = 20            # A 类重复测量次数
S_REL = 0.05          # 重复测量样本标准差（相对）
N_FIT = 20            # 拟合点数
SIGMA_Y = 0.5         # 单点测量标准差（counts，等权）
SIG_PHI = 0.000361482  # 附录 AQ 的相位分辨（rad）
NE_LINE_AQ = 1.07979e15  # 附录 AQ：0.001 rad 折合线积分密度（m^-2）= 1e-3 rad 对应
K_COVER = 2.0         # 扩展不确定度覆盖因子（约 95%）
F_SIG = 1000.0        # Hz，采样率（混叠例）
F_TRUE = 1300.0       # Hz，真实信号频率

print('==== 声明输入 ====')
print('n_e0 = %.6g m^-3, T_e = %.6g eV, B = %.6g T, a = %.6g m, R = %.6g m'
      % (NE0, TE_EV, B0, A_MINOR, R0))
print('N_PH = %d（每个空间点的光子数）, 单点统计误差 %.4g（相对）, B 类系统误差 %.4g（相对，均匀分布）'
      % (N_PH, REL_STAT, REL_SYS))
print('N_REP = %d, 重复测量样本标准差 %.4g（相对）' % (N_REP, S_REL))
print('N_FIT = %d, 等权单点 sigma_y = %.4g counts' % (N_FIT, SIGMA_Y))
print('sigma_phi(AQ) = %.6g rad；0.001 rad 折合线积分密度 = %.6g m^-2'
      % (SIG_PHI, NE_LINE_AQ))
print()

# ================= AU.1 从信号到误差棒 =================
print('==== AU.1 从信号到误差棒：一根棒子里的三笔账 ====')
u_stat = REL_STAT
u_sys = REL_SYS
u_c = math.hypot(u_stat, u_sys)
U_exp = K_COVER * u_c
print('单点统计项 u_stat = %.6g（相对）' % u_stat)
print('B 类系统项 u_sys = %.6g（相对，按均匀分布直接取半宽）' % u_sys)
print('合成相对不确定度 u_c = sqrt(u_stat^2 + u_sys^2) = %.6g -> %.5f%%'
      % (u_c, 100 * u_c))
print('扩展不确定度 U = %.1f u_c = %.6g -> %.5f%%（约 95%%）'
      % (K_COVER, U_exp, 100 * U_exp))
print('判据对照：U = %.5f%% -> 满足 15%% 判据；超过 10%% 判据 %.5f 个百分点'
      % (100 * U_exp, 100 * U_exp - 10.0))
# 要压到 U <= 5% 需要多少光子
U_target = 0.05
u_c_target = U_target / K_COVER
u_stat_target = math.sqrt(u_c_target ** 2 - u_sys ** 2)
n_ph_need = (u_stat / u_stat_target) ** 2 * N_PH
print('若要求 U <= 5%%（k=2）: u_c <= %.6g, 统计项需 <= %.6g（现为 %.6g）'
      % (u_c_target, u_stat_target, u_stat))
print('  每点光子数需 x%.6g -> %.6g 个 -> 取整 %d 个'
      % (n_ph_need / N_PH, n_ph_need, int(math.ceil(n_ph_need))))
# 单点绝对不确定度（n_e = 1e20）
abs_ne = NE0 * u_c
print('折合绝对不确定度：n_e = %.6g +/- %.6g m^-3' % (NE0, abs_ne))
print()

# ================= AU.2 A 类与 B 类 =================
print('==== AU.2 两条腿：A 类统计、B 类系统、合成与扩展 ====')
uA = S_REL / math.sqrt(N_REP)
uB = REL_SYS / math.sqrt(3.0)
u_c2 = math.hypot(uA, uB)
nu_eff = u_c2 ** 4 / (uA ** 4 / (N_REP - 1))
print('A 类：u_A = s/sqrt(n) = %.6g / sqrt(%d) = %.6g -> %.5f%%'
      % (S_REL, N_REP, uA, 100 * uA))
print('B 类：u_B = a/sqrt(3) = %.6g / sqrt(3) = %.6g -> %.5f%%'
      % (REL_SYS, uB, 100 * uB))
print('合成：u_c = %.6g -> %.5f%%；U = %.1f u_c = %.5f%%'
      % (u_c2, 100 * u_c2, K_COVER, 100 * K_COVER * u_c2))
print('Welch-Satterthwaite 有效自由度 nu_eff = %.6g（A 类 nu=19，B 类 nu=inf）' % nu_eff)
# t 值
def t_at(nu, p=0.975):
    # 简化：用表值插值的常见点（数值固定，便于复算）
    table = {1: 12.7062, 2: 4.30265, 3: 3.18245, 4: 2.77645, 5: 2.57058,
             10: 2.22814, 18: 2.10092, 19: 2.09302, 20: 2.08596, 30: 2.04227,
             60: 2.00030, 120: 1.97993}
    keys = sorted(table)
    if nu in table:
        return table[nu]
    if nu > max(keys):
        return 1.95996  # 正态极限
    lo = max([k for k in keys if k <= nu])
    hi = min([k for k in keys if k >= nu])
    if lo == hi:
        return table[lo]
    w = (nu - lo) / float(hi - lo)
    return table[lo] * (1 - w) + table[hi] * w
k_t = t_at(nu_eff)
print('按 nu_eff 取 t(95%%) = %.6g；用 t 而非 2 时 U = %.5f%%（差 %.4f 个百分点）'
      % (k_t, 100 * k_t * u_c2, 100 * (k_t - K_COVER) * u_c2))
print('若把系统项当均匀分布半宽 a 直接当标准差（少除 sqrt(3)）：u_c = %.5f%%（高估 %.4f 个百分点）'
      % (100 * math.hypot(uA, REL_SYS),
         100 * (math.hypot(uA, REL_SYS) - u_c2)))
print()

# ================= AU.3 误差传播 =================
print('==== AU.3 误差传播：加减、乘除、相关项 ====')
# 减法：n_e = n_tot - n_imp
n_tot, r_tot = 1.0e20, 0.05
n_imp, r_imp = 5.0e18, 0.05
s_tot, s_imp = n_tot * r_tot, n_imp * r_imp
n_net = n_tot - n_imp
s_sub = math.hypot(s_tot, s_imp)
print('减法：n_net = %.6g - %.6g = %.6g m^-3' % (n_tot, n_imp, n_net))
print('  sigma = sqrt(%.6g^2 + %.6g^2) = %.6g -> 相对 %.5f%%（原 5%%，放大 %.4f 倍）'
      % (s_tot, s_imp, s_sub, 100 * s_sub / n_net, (s_sub / n_net) / 0.05))
# 相关项：r = 0.9
for r in (0.0, 0.5, 0.9, 0.99, 1.0):
    s_corr = math.sqrt(s_tot ** 2 + s_imp ** 2 - 2 * r * s_tot * s_imp)
    print('  相关 r = %.2f: sigma = %.6g -> 相对 %.5f%%' % (r, s_corr, 100 * s_corr / n_net))
# 乘除：两个 5%
rel_prod = math.hypot(0.05, 0.05)
print('乘法：z = x y，两个 5%% -> 相对 %.6g -> %.5f%%' % (rel_prod, 100 * rel_prod))
# 幂律
print('幂律：z = x^2，5%% -> %.5f%%（翻倍）' % (100 * 2 * 0.05)) 
print('幂律：z = x^(1/2)，5%% -> %.5f%%（减半）' % (100 * 0.5 * 0.05))
# AQ 跨附录：相位 -> 线积分密度
print('跨附录（AQ 口径）：sigma_phi = %.6g rad -> 线积分密度误差 = %.6g m^-2'
      % (SIG_PHI, NE_LINE_AQ * SIG_PHI / 1e-3))
print()

# ================= AU.4 最小二乘与卡方 =================
print('==== AU.4 最小二乘与卡方：把直线拟合进带误差的点 ====')
xs = [i / float(N_FIT - 1) for i in range(N_FIT)]
xbar = sum(xs) / len(xs)
Sxx = sum((x - xbar) ** 2 for x in xs)
sig_b = SIGMA_Y / math.sqrt(Sxx)
sig_a = SIGMA_Y * math.sqrt(1.0 / N_FIT + xbar ** 2 / Sxx)
print('x_i = i/(N-1), i=0..%d；xbar = %.6g, S_xx = %.6g' % (N_FIT - 1, xbar, Sxx))
print('斜率不确定度 sigma_b = sigma_y/sqrt(S_xx) = %.6g counts/单位 x' % sig_b)
print('截距不确定度 sigma_a = sigma_y sqrt(1/N + xbar^2/S_xx) = %.6g counts' % sig_a)
print('相关系数 corr(a,b) = -xbar/sqrt(xbar^2 + S_xx/N) = %.6g'
      % (-xbar / math.sqrt(xbar ** 2 + Sxx / N_FIT)))
# 卡方
dof = N_FIT - 2
print('自由度 dof = N - 2 = %d；chi2/dof 的固有散布 sigma = sqrt(2/dof) = %.6g -> %.4f%%'
      % (dof, math.sqrt(2.0 / dof), 100 * math.sqrt(2.0 / dof)))
print('95%% 区间（dof=%d 附近）：chi2 小于 %.6g 或大于 %.6g 都说明模型/误差有问题'
      % (dof, 9.390, 28.869))
print('举例：chi2 = %.1f -> chi2/dof = %.6g（好）；chi2 = %.1f -> chi2/dof = %.6g（坏）'
      % (18.0, 18.0 / dof, 36.0, 36.0 / dof))
# 加权平均
y1, s1, y2, s2 = 10.0, 1.0, 12.0, 2.0
w1, w2 = 1.0 / s1 ** 2, 1.0 / s2 ** 2
ym = (w1 * y1 + w2 * y2) / (w1 + w2)
sm = 1.0 / math.sqrt(w1 + w2)
chi2 = w1 * (y1 - ym) ** 2 + w2 * (y2 - ym) ** 2
print('加权平均：10.0 +/- 1.0 与 12.0 +/- 2.0 -> %.6g +/- %.6g；chi2 = %.6g (dof=1)'
      % (ym, sm, chi2))
print()

# ================= AU.5 置信区间与小样本 =================
print('==== AU.5 置信区间：为什么小样本必须用 t ====')
for n in (5, 20, 200):
    nu = n - 1
    kk = t_at(nu)
    half = kk * S_REL / math.sqrt(n)
    print('n = %3d (nu=%3d): t = %.6g, 95%% 区间半宽 = %.5f%%'
          % (n, nu, kk, 100 * half))
print('样本量翻 4 倍才能把区间压窄一半：n=20 -> 半宽 %.5f%%，n=80 -> %.5f%%'
      % (100 * t_at(19) * S_REL / math.sqrt(20),
         100 * t_at(79) * S_REL / math.sqrt(80)))
print('AU.2 的有效自由度 %.6g 落在 19 与 inf 之间，故 2.0 是偏乐观的覆盖因子' % nu_eff)
print()

# ================= AU.6 实验设计与灵敏度 =================
print('==== AU.6 实验设计：把机时花在敏感的方向 ====')
# 两条弦：c_i 为误差系数（绝对）
c1, c2 = 5.0e18, 2.5e17
frac1 = c1 / (c1 + c2)
print('两条弦的误差系数 c1 = %.6g, c2 = %.6g -> 最优时间分配 %.5f%% / %.5f%%'
      % (c1, c2, 100 * frac1, 100 * (1 - frac1)))
u_opt = c1 + c2
u_eq = math.sqrt(2.0) * math.hypot(c1, c2)
print('最优分配（t_i ∝ c_i）：u ∝ c1 + c2 = %.6g / sqrt(T)' % u_opt)
print('平均分配：u ∝ sqrt(2(c1^2 + c2^2)) = %.6g / sqrt(T)，最优分配好 %.5f 倍'
      % (u_eq, u_eq / u_opt))
# 全因子 2^3
print('2^3 全因子：8 次运行 -> 3 个主效应 + 3 个二因子交互 + 1 个三因子交互')
print('  主效应方差 Var = sigma^2/2 = 0.5 sigma^2；一次一因子（每因子两水平各 2 次）同为 0.5 sigma^2，但拿不到交互项')
# 灵敏度系数例
V, dV = 1000.0, 10.0
sens = 2.0  # T_e 对 V 的幂次（例：Te ∝ V^2）
print('灵敏度：V = %.6g +/- %.6g（%.4f%%），若 Te ∝ V^%.1f，则贡献 = %.6g -> %.5f%%'
      % (V, dV, 100 * dV / V, sens, sens * dV / V, 100 * sens * dV / V))
print()

# ================= AU.7 蒙特卡洛与常见陷阱 =================
print('==== AU.7 蒙特卡洛、样本量的不确定度与三类陷阱 ====')
for n in (20, 200, 2000):
    print('N = %4d：标准差本身的相对不确定度 1/sqrt(2(N-1)) = %.5f%%'
          % (n, 100.0 / math.sqrt(2 * (n - 1))))
print('混叠：真实 f = %.6g Hz，采样 f_s = %.6g Hz -> 假峰 |f - f_s| = %.6g Hz'
      % (F_TRUE, F_SIG, abs(F_TRUE - F_SIG)))
print('  要避开混叠需 f_s >= 2 f = %.6g Hz（Nyquist）' % (2 * F_TRUE))
# 去卷积：仪器 1.0 nm，真实 0.5 nm
inst, true = 1.0, 0.5
obs = math.hypot(inst, true)
print('分辨率合成：仪器 %.4g nm 与真实线宽 %.4g nm -> 观测 %.6g nm（方和根）'
      % (inst, true, obs))
print('  若直接反解真实线宽：sqrt(obs^2 - inst^2) = %.6g nm（对 obs 的误差放大 %.4f 倍）'
      % (math.sqrt(obs ** 2 - inst ** 2), obs / true))
# 基线扣除
sig_s, sig_b = 10.0, 10.0
print('基线扣除：信号 %.4g counts、基线 %.4g counts、各自 sigma = %.4g -> 扣后 sigma = %.6g（相对 %.4f%%）'
      % (100.0, 1000.0, sig_s, math.hypot(sig_s, sig_b),
         100 * math.hypot(sig_s, sig_b) / 100.0))
print('3 sigma 检测限：3 x %.6g counts = %.6g counts；若灵敏度 100 counts 对应 1e16 m^-3，则检测限 %.6g m^-3'
      % (math.hypot(sig_s, sig_b), 3 * math.hypot(sig_s, sig_b),
         3 * math.hypot(sig_s, sig_b) / 100.0 * 1e16))
print()

# ================= AU.8 练习 =================
print('==== AU.8 练习与解答 ====')
# 1) 三个独立量各 1% 相乘
v1 = math.sqrt(3) * 0.01
print('[练1] 三个独立量各 1%% 相乘：相对不确定度 = sqrt(3) x 1%% = %.5f%%' % (100 * v1))
# 2) 减法
a1, ar1 = 1.0e20, 0.03
a2, ar2 = 1.0e19, 0.03
s1_, s2_ = a1 * ar1, a2 * ar2
net = a1 - a2
snet = math.hypot(s1_, s2_)
print('[练2] n = %.6g +/- %.1f%% 减 %.6g +/- %.1f%%：net = %.6g, sigma = %.6g, 相对 %.5f%%'
      % (a1, 100 * ar1, a2, 100 * ar2, net, snet, 100 * snet / net))
# 3) 5 点、s=2%
n5, s5 = 5, 0.02
print('[练3] n=%d, s=%.4g, t(95%%,nu=4)=2.77645 -> 95%% 半宽 = %.5f%%'
      % (n5, s5, 100 * 2.77645 * s5 / math.sqrt(n5)))
# 4) 点数
need2 = (0.05 / math.sqrt(0.025 ** 2 - 0.02 ** 2)) ** 2
print('[练4] U 从 %.5f%% 压到 5%%：统计项需 <= 1.5%%，每点光子数 = (5/1.5)^2 x %d = %.6g -> %d 个'
      % (100 * K_COVER * math.hypot(0.05, 0.02), N_PH, need2 * N_PH,
         int(math.ceil(need2 * N_PH))))
# 5) 加权平均一致性
print('[练5] 上例 chi2 = %.6g < 3.84（dof=1 的 95%% 上限），两组测量一致' % chi2)
# 6) 混叠
print('[练6] f = %.6g Hz 用 f_s = %.6g Hz 采样：假峰 %.6g Hz；改用 f_s = 2000 Hz 时假峰 %.6g Hz；'
      '最小不混叠采样率 %.6g Hz' % (F_TRUE, F_SIG, abs(F_TRUE - F_SIG),
                                     abs(F_TRUE - 2000.0), 2 * F_TRUE))

# ================= AU.9 补充派生量（正文引用，必须在此 print） =================
print()
print('==== AU.9 补充派生量 ====')
N_MC = 100000
print('蒙特卡洛抽样次数（声明）N_MC = %d' % N_MC)
s_of_s = 100.0 / math.sqrt(2 * (20 - 1))
print('N=20 时标准差的相对不确定度 = %.5f%% -> 5%% 的严格写法 5%% +/- %.6g%%'
      % (s_of_s, 5.0 * s_of_s / 100.0))
half_t4 = 2.77645 * 0.02 / math.sqrt(5)
half_k2 = 2.0 * 0.02 / math.sqrt(5)
print('[练3] t 半宽 = %.6g%%；误用 k=2 得 %.6g%%；比 t 口径窄 %.5f%%'
      % (100 * half_t4, 100 * half_k2, 100 * (1 - half_k2 / half_t4)))
print('小样本 t 与正态极限之比 = %.6g / %.6g = %.6g'
      % (2.77645, 1.95996, 2.77645 / 1.95996))
ratio_t = 2.77645 / 1.95996
print('用 k=2 代替 t 时区间报窄 = %.5f%%' % (100 * (ratio_t - 1)))
s_3pct_1 = 1.0e20 * 0.03
s_3pct_2 = 1.0e19 * 0.03
print('[练2] 两个输入的绝对误差 = %.6g 与 %.6g m^-3' % (s_3pct_1, s_3pct_2))
