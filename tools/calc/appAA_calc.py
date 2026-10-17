# -*- coding: utf-8 -*-
"""附录 AA（替代位形：仿星器、球形托卡马克与场反位形）的唯一数值来源。

口径：D1 托卡马克沿用全书口径；三种替代位形各声明一组「装置级」参数（量级取自
W7-X / NSTX-U / FRC 实验装置），工程侧与定标系数一律声明。所有打印值即正文可引用值。
"""
import math

mu0 = 1.25663706212e-6
e = 1.602176634e-19
keV = 1e3 * e

# ---------------- DECLARED ----------------
# D1 托卡马克（全书口径）
R1, a1, B1 = 6.2, 2.0, 5.3
n1, T1 = 1e20, 9.194 * keV
Ip1 = 9.0            # MA（声明）
q95_1, q0_1 = 3.00, 1.05
tauE1 = 2.4762       # s，附录 S 的 IPB98(y,2) 值
# 球形托卡马克（NSTX-U 量级）
R2, a2, B2 = 0.93, 0.62, 1.0
n2, T2 = 1e20, 1.0 * keV
Ip2 = 2.0            # MA
# 仿星器（W7-X 量级）
R3, a3, B3 = 5.5, 0.53, 2.5
n3, T3 = 8e19, 2.0 * keV
iota3 = 0.85         # iota/2pi
P3 = 10.0            # MW 吸收功率
eps_eff3 = 0.013     # 优化位形的有效纹波（声明）
# 场反位形（FRC）
Bext, rs, Lfrc = 0.5, 0.4, 8.0
Tfrc = 1.0 * keV
# 定标与工程系数（声明）
ISS04 = (0.134, 2.28, 0.64, -0.61, 0.54, 0.84, 0.41)   # C, a, R, P, n, B, iota
shield_Z = 2.0723    # m，附录 Z 算出的中子屏蔽厚度

DAO = '《等离子体物理自学教材》附录 AA 数值脚本'
print('=' * 72)
print(DAO, '—— 全部数字的唯一来源')
print('=' * 72)

# ---------------- AA.1 三种位形的一张地图 ----------------
def beta_of(n, T, B):
    p = n * T
    return 2 * mu0 * p / (B * B)


def betan_of(beta, a, B, Ip):
    return beta * 100 * a * B / Ip


A1, A2, A3 = R1 / a1, R2 / a2, R3 / a3
e1, e2, e3 = a1 / R1, a2 / R2, a3 / R3
b1 = beta_of(n1, T1, B1)
b2 = beta_of(n2, T2, B2)
b3 = beta_of(n3, T3, B3)
bn1 = betan_of(b1, a1, B1, Ip1)
bn2 = betan_of(b2, a2, B2, Ip2)
bn3 = betan_of(b3, a3, B3, 1.0)      # 仿星器无环向电流，用 1 MA 只为给出量级对照
p1, p2, p3 = n1 * T1, n2 * T2, n3 * T3
print('\n[AA.1] 三种位形的地图')
print('  D1 托卡马克 : A=R/a=%.5g, eps=a/R=%.5g, p=%.5g Pa, beta=%.5g (%.4g %%), beta_N=%.5g'
      % (A1, e1, p1, b1, b1 * 100, bn1))
print('  球形托卡马克: A=%.5g, eps=%.5g, p=%.5g Pa, beta=%.5g (%.4g %%), beta_N=%.5g'
      % (A2, e2, p2, b2, b2 * 100, bn2))
print('  仿星器      : A=%.5g, eps=%.5g, p=%.5g Pa, beta=%.5g (%.4g %%), beta_N=%.5g（1 MA 对照）'
      % (A3, e3, p3, b3, b3 * 100, bn3))
print('  beta(ST)/beta(D1) = %.5g' % (b2 / b1))
print('  1/beta(D1) = %.5g（FRC 的 beta=1 相对 D1 的倍数，只在口径上成立）' % (1 / b1))
print('  eps(ST)/eps(D1) = %.5g, eps(D1)/eps(仿星器) = %.5g'
      % (e2 / e1, e1 / e3))

# ---------------- AA.2 旋转变换、安全因子与磁剪切 ----------------
q3 = 1.0 / iota3
shear1 = 2 - 2 * q0_1 / q95_1
print('\n[AA.2] 旋转变换与磁剪切')
print('  仿星器 iota/2pi = %.4g -> q = 1/iota = %.5g' % (iota3, q3))
print('  D1 的 q0=%.4g, q95=%.4g -> 剪切 s = 2-2q0/q95 = %.5g' % (q0_1, q95_1, shear1))
print('  仿星器的 iota 由线圈决定（无等离子体电流），托卡马克的 q 由 I_p 决定')

# ---------------- AA.3 新经典输运：sqrt(eps) 与 eps_eff^(3/2) ----------------
sq_e1 = math.sqrt(e1)
eff32 = eps_eff3 ** 1.5
ratio_eff = eps_eff3 / e3
print('\n[AA.3] 新经典输运的几何因子')
print('  sqrt(eps_D1) = %.5g' % sq_e1)
print('  eps_eff = %.5g, eps_eff^(3/2) = %.6g' % (eps_eff3, eff32))
print('  eps_eff/eps(仿星器) = %.5g = %.5g %%' % (ratio_eff, ratio_eff * 100))
print('  sqrt(eps_D1)/eps_eff^(3/2) = %.6g' % (sq_e1 / eff32))
print('  eps_eff 减半 -> eps_eff^(3/2) 变为 %.6g 倍' % (0.5 ** 1.5))

# ---------------- AA.4 球形托卡马克：中心柱的空间代价 ----------------
space2 = R2 - a2
space1 = R1 - a1
print('\n[AA.4] 球形托卡马克的空间代价')
print('  D1 的中心柱空间 R-a = %.5g m（附录 Z 的屏蔽需要 %.5g m）' % (space1, shield_Z))
print('  ST 的中心柱空间 R-a = %.5g m' % space2)
print('  比屏蔽需求 2.0723 m 少 %.5g 倍' % (shield_Z / space2))
print('  B_p(ST) = mu0 I_p/(2 pi a) = %.5g T' % (mu0 * Ip2 * 1e6 / (2 * math.pi * a2)))
bp2 = mu0 * Ip2 * 1e6 / (2 * math.pi * a2)
print('  beta_p(ST) = 2 mu0 p / B_p^2 = %.5g' % (2 * mu0 * p2 / bp2 ** 2))
fbs2 = math.sqrt(e2) * 2 * mu0 * p2 / bp2 ** 2
print('  自举份额量级 C_bs sqrt(eps) beta_p = %.5g = %.5g %%（C_bs 取 1）'
      % (fbs2, fbs2 * 100))

# ---------------- AA.5 场反位形：beta=1 意味着什么 ----------------
p_frc = Bext ** 2 / (2 * mu0)
n_frc = p_frc / Tfrc
I_frc = 2 * math.pi * rs * Bext / mu0
elong = Lfrc / (2 * rs)
print('\n[AA.5] 场反位形')
print('  B_ext = %.2f T -> beta=1 时的最大压力 p = B^2/(2 mu0) = %.5g Pa' % (Bext, p_frc))
print('  T = 1.0 keV 时对应密度 n = %.5g m^-3 = %.5g cm^-3' % (n_frc, n_frc / 1e6))
print('  反平行电流 I = 2 pi r_s B/mu0 = %.6g A = %.5g MA' % (I_frc, I_frc / 1e6))
print('  长细比 L/(2 r_s) = %.5g' % elong)
print('  若 B_ext 加倍：p 变为 %.5g 倍、n 变为 %.5g 倍' % (4.0, 4.0))

# ---------------- AA.6 定标：ISS04 与三重积 ----------------
C, ea, eR, eP, en, eB, ei = ISS04
tauE3 = (C * a3 ** ea * R3 ** eR * P3 ** eP * (n3 / 1e19) ** en
         * B3 ** eB * iota3 ** ei)
triple3 = n3 * (T3 / keV) * tauE3          # n[1e20? 用 m^-3] * T[keV] * tau
triple_thr = 3.1e21                        # 附录 U 的点火三重积阈值
print('\n[AA.6] 定标与三重积')
print('  ISS04: C=%.4g, 指数 a=%.4g R=%.4g P=%.4g n=%.4g B=%.4g iota=%.4g' % ISS04)
print('  tau_E(仿星器) = %.6g s' % tauE3)
print('  tau_E(D1, IPB98) = %.6g s -> 比值 = %.5g' % (tauE1, tauE3 / tauE1))
print('  三重积 n T tau_E = %.6g（阈值 %.4g）-> 差 %.5g 倍'
      % (triple3, triple_thr, triple_thr / triple3))
print('  iota 从 0.85 提到 1.00: tau_E 变为 %.5g 倍' % ((1.0 / iota3) ** ei))

# ---------------- AA.7 练习解答 ----------------
print('\n[AA.7] 练习解答')
print('  1) ST 中心柱 0.31 m，屏蔽需要 2.0723 m -> 差 %.5g 倍（无法容纳）' % (shield_Z / space2))
print('  2) iota 0.85->1.00: tau_E x %.5g（+%.5g %%）'
      % ((1.0 / iota3) ** ei, ((1.0 / iota3) ** ei - 1) * 100))
print('  3) B_ext 0.5->1.0 T: p x 4 = %.5g Pa, n x 4 = %.5g m^-3' % (4 * p_frc, 4 * n_frc))
print('  4) beta(ST)/beta(D1) = %.5g' % (b2 / b1))
print('  5) eps_eff 0.013->0.0065: eps_eff^(3/2) x %.6g' % (0.5 ** 1.5))
print('  6) ST 的 B 1.0->1.5 T: 因子 1.5^2 = %.5g, beta = %.5g (%.4g %%), beta_N = %.5g'
      % (1.5 ** 2, b2 / 1.5 ** 2, b2 / 1.5 ** 2 * 100, bn2 * 1.5 / 1.5 ** 2))

# ---------------- 对照 ----------------
print('\n[对照]')
print('  三种位形的 A: %.4g / %.4g / %.4g' % (A1, A2, A3))
print('  三种位形的 eps: %.4g / %.4g / %.4g' % (e1, e2, e3))
print('  三种位形的 beta[%%]: %.4g / %.4g / %.4g' % (b1 * 100, b2 * 100, b3 * 100))
print('  beta_N: %.4g / %.4g / %.4g' % (bn1, bn2, bn3))
print('  sqrt(eps_D1) = %.4g, eps_eff^(3/2) = %.4g' % (sq_e1, eff32))
print('  ISS04 tau_E = %.4g s, D1 IPB98 tau_E = %.4g s' % (tauE3, tauE1))
print('  FRC: p = %.4g Pa, n = %.4g m^-3, I = %.4g MA' % (p_frc, n_frc, I_frc / 1e6))
print('\n完成：以上即正文可以引用的全部数值。')
