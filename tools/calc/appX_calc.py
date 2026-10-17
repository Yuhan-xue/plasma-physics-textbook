# -*- coding: utf-8 -*-
"""附录 X（辐射输运与辐射损失）的数值唯一来源。

输出 tmp/appX_calc.txt；正文 src/appendix_x_radiation.tex 与门禁
check_P29_appendix_x_radiation() 只允许引用本脚本打印的数字。

工作点 D1 与附录 S/U/V/W 同口径：
  n_e = 1.0e20 m^-3, T_e = 9.194 keV, B = 5.3 T, a = 2.0 m, R = 6.2 m,
  I_p = 15 MA, V = 832.2 m^3, P_fus = 556.5 MW, P_alpha = 111.3 MW,
  W_th = 371.4 MJ, Z_eff = 2, L_W(10 keV) 取 1e-32 W m^3 的量级（声明见正文）。

简化声明：只算光学薄极限与黑体上限两个边界，不给输运方程的数值解；
轫致辐射用非相对论 + 一阶相对论修正；线辐射只取一个声明的量级值。
本脚本不引用任何外部结论数值，全部由常数与工作点现算。
"""
import io
import math
import os

out = []


def p(fmt, *args):
    s = fmt % args if args else fmt
    out.append(s)
    print(s)


def head(t):
    p('')
    p('=' * 72)
    p(t)
    p('=' * 72)


# ---------------------------------------------------------------- 常数
e = 1.602176634e-19          # C
eps0 = 8.8541878128e-12      # F/m
m_e = 9.1093837015e-31       # kg
c = 2.99792458e8             # m/s
kB = 1.380649e-23            # J/K
mu0 = 4.0e-7 * math.pi       # H/m
sigma_T = 6.6524587321e-29   # m^2
h = 6.62607015e-34           # J s
keV = 1.0e3 * e

# ---------------------------------------------------------------- 工作点 D1
n_e = 1.0e20                 # m^-3
T_keV = 9.194                # keV
B = 5.3                      # T
a = 2.0                      # m
R = 6.2                      # m
Ip_MA = 15.0
V = 832.2                    # m^3
Z_eff = 2.0
P_fus = 556.5e6              # W
P_alpha = 111.3e6            # W
W_th = 371.4e6               # J

T = T_keV * keV              # J
T_eV = T_keV * 1.0e3         # eV
L_W = 1.0e-32                # W m^3，钨在 10 keV 量级的冷却率系数（声明）
f_W = 1.0e-5                 # n_W / n_e

head('0. 工作点 D1')
p('n_e = %.4g m^-3, T_e = %.4g keV, B = %.4g T, a = %.4g m, R = %.4g m' %
  (n_e, T_keV, B, a, R))
p('V = %.5g m^3, Z_eff = %.4g, P_fus = %.5g W, P_alpha = %.5g W, W_th = %.5g J'
  % (V, Z_eff, P_fus, P_alpha, W_th))

head('1. 轫致辐射（光学薄，主要辐射损失）')
C_BR = 5.35e-37              # P_br = C_BR * Z_eff * n_e^2 * sqrt(T[keV])  W/m^3
P_br = C_BR * Z_eff * n_e ** 2 * math.sqrt(T_keV)
P_BR = P_br * V
p('P_br = %.5g * %.4g * (%.4g)^2 * sqrt(%.5g) = %.5g W/m^3 = %.5g kW/m^3'
  % (C_BR, Z_eff, n_e, T_keV, P_br, P_br / 1.0e3))
p('全等离子体 P_BR = P_br * V = %.5g W = %.5g MW' % (P_BR, P_BR / 1.0e6))
p('P_BR / P_alpha = %.5g%%' % (100.0 * P_BR / P_alpha))
p('P_BR / P_fus   = %.5g%%' % (100.0 * P_BR / P_fus))

# 一阶相对论修正（声明为一阶，只报量级）
rel = 1.0 + 2.0 * T / (m_e * c * c)
p('T_e / (m_e c^2) = %.5g；一阶相对论修正因子 1 + 2T/(m_e c^2) = %.5g'
  % (T / (m_e * c * c), rel))
P_br_rel = P_br * rel
P_BR_rel = P_br_rel * V
p('修正后 P_br = %.5g W/m^3（提高 %.4g%%），P_BR = %.5g MW'
  % (P_br_rel, 100.0 * (rel - 1.0), P_BR_rel / 1.0e6))

# 平均光子能量与波长
E_ph = T_eV                  # 量级：平均光子能量 ~ T_e
lam_ph = h * c / (E_ph * e)
p('平均光子能量量级 ~ T_e = %.5g eV；波长 = %.5g m = %.5g nm'
  % (E_ph, lam_ph, lam_ph * 1.0e9))

head('2. 回旋辐射：光学薄上限')
nu_ce = e * B / (2.0 * math.pi * m_e)
u_B = B * B / (2.0 * mu0)
beta2 = 2.0 * T / (m_e * c * c)
P1 = (4.0 / 3.0) * sigma_T * c * beta2 * u_B
P_cyc_thin = n_e * P1
P_CYC_thin = P_cyc_thin * V
p('f_ce = eB/(2 pi m_e) = %.5g Hz = %.5g GHz' % (nu_ce, nu_ce / 1.0e9))
p('磁能密度 u_B = B^2/(2 mu0) = %.5g J/m^3' % u_B)
p('beta_perp^2 = 2T/(m_e c^2) = %.5g' % beta2)
p('单电子 Larmor 功率 P_1 = (4/3) sigma_T c beta^2 u_B = %.5g W' % P1)
p('光学薄 P_cyc = n_e P_1 = %.5g W/m^3 = %.5g MW/m^3' % (P_cyc_thin, P_cyc_thin / 1.0e6))
p('全等离子体 P_CYC(光学薄) = %.5g W = %.5g MW' % (P_CYC_thin, P_CYC_thin / 1.0e6))
p('P_CYC(光学薄) / P_BR = %.5g' % (P_CYC_thin / P_BR))
p('P_CYC(光学薄) / P_fus = %.5g  <-- 大于 1，物理上不可能' % (P_CYC_thin / P_fus))

head('3. 回旋辐射：黑体（Rayleigh-Jeans）上限')
A_wall = 4.0 * math.pi ** 2 * R * a
F_RJ = 2.0 * math.pi * nu_ce ** 3 * T / (c * c)   # T 已是能量[J]，k_B T = T
P_RJ = F_RJ * A_wall
p('第一壁面积 A = 4 pi^2 R a = %.5g m^2' % A_wall)
p('RJ 通量（按一个谐波带宽 ~ f_ce 估）F = 2 pi nu^3 kB T / c^2 = %.5g W/m^2' % F_RJ)
p('P_CYC(黑体上限) = %.5g W = %.5g MW' % (P_RJ, P_RJ / 1.0e6))
p('再吸收把光学薄上限压低 %.5g 倍' % (P_CYC_thin / P_RJ))
p('P_CYC(黑体上限) / P_BR = %.5g%%' % (100.0 * P_RJ / P_BR))

head('4. 自由-自由吸收与光学深度（轫致辐射为什么逃得掉）')
# cgs: alpha_ff = 3.7e8 Z^2 n_e n_i T^-1/2 nu^-3 (1 - exp(-h nu / kT))  [cm^-1]
n_e_cgs = n_e * 1.0e-6
Z2_ne_ni = Z_eff * n_e_cgs * n_e_cgs
T_K = T / kB
x = h * nu_ce / T            # T 已是能量[J]，k_B T = T
alpha_cgs = 3.7e8 * Z2_ne_ni * T_K ** -0.5 * nu_ce ** -3 * (1.0 - math.exp(-x))
alpha_SI = alpha_cgs * 100.0
L_chord = 2.0 * a
tau_ff = alpha_SI * L_chord
p('n_e = %.5g cm^-3, T = %.5g K, h nu / kT = %.5g' % (n_e_cgs, T_K, x))
p('alpha_ff = %.5g cm^-1 = %.5g m^-1' % (alpha_cgs, alpha_SI))
p('沿弦长 L = 2a = %.5g m 的光学深度 tau_ff = %.5g' % (L_chord, tau_ff))
p('tau_ff << 1：轫致辐射光学薄，几乎全部逃逸')
lam_mfp = 1.0 / alpha_SI
p('平均自由程 1/alpha_ff = %.5g m（比装置尺寸大 %.5g 倍）' % (lam_mfp, lam_mfp / L_chord))

head('5. 线辐射与辐射崩塌阈值（钨，L_W 取 10 keV 量级）')
n_W = f_W * n_e
P_line = n_e * n_W * L_W
P_LINE = P_line * V
p('L_W = %.4g W m^3（声明：钨在 10 keV 的量级值）' % L_W)
p('n_W/n_e = %.4g -> n_W = %.5g m^-3' % (f_W, n_W))
p('P_line = n_e n_W L_W = %.5g W/m^3 = %.5g kW/m^3' % (P_line, P_line / 1.0e3))
p('全等离子体 P_LINE = %.5g W = %.5g MW' % (P_LINE, P_LINE / 1.0e6))
p('P_LINE / P_BR = %.5g%%' % (100.0 * P_LINE / P_BR))
f_W_star = P_alpha / (V * n_e * L_W) / n_e
p('辐射崩塌阈值（P_LINE = P_alpha）：n_W* = %.5g m^-3，f_W* = %.5g%%'
  % (f_W_star * n_e, 100.0 * f_W_star))

head('6. 辐射在功率账里的位置')
P_rad = P_BR + P_LINE
p('P_rad = P_BR + P_LINE = %.5g MW' % (P_rad / 1.0e6))
p('P_rad / P_fus = %.5g%%' % (100.0 * P_rad / P_fus))
p('自加热余量 P_alpha / P_rad = %.5g' % (P_alpha / P_rad))
tau_rad = W_th / P_rad
p('辐射冷却时间 tau_rad = W_th / P_rad = %.5g s' % tau_rad)
p('对照：附录 V 的电流猝灭时间 tau_CQ = 76.9074 ms')
flux_wall = P_BR / A_wall
p('轫致辐射摊到第一壁：P_BR / A = %.5g W/m^2 = %.5g kW/m^2'
  % (flux_wall, flux_wall / 1.0e3))

head('8. 练习解答量与对照比例（正文引用的每一个派生值都在这里）')
p_Z15 = C_BR * 1.5 * n_e ** 2 * math.sqrt(T_keV) * V
p('练习1  Z_eff=1.5 -> P_BR = %.6g MW' % (p_Z15 / 1.0e6))
T2 = 20.0
P_br2 = C_BR * Z_eff * n_e ** 2 * math.sqrt(T2)
p('练习2  T_e=20 keV -> P_br = %.6g kW/m^3（比值 sqrt(20/9.194) = %.5g，涨幅 %.5g%%）'
  % (P_br2 / 1.0e3, math.sqrt(T2 / T_keV), 100.0 * (math.sqrt(T2 / T_keV) - 1.0)))
n_W_br = (P_BR / V) / (n_e * L_W)
p('练习3  P_line = P_BR 所需 n_W = %.6g m^-3，n_W/n_e = %.6g%%'
  % (n_W_br, 100.0 * n_W_br / n_e))
T_bb = T_keV * (P_BR / 1.0e6) / (P_RJ / 1.0e6)
p('练习4  黑体上限达到 %.6g MW 所需亮度温度 = %.6g keV' % (P_BR / 1.0e6, T_bb))
p('练习5  tau_ff = 1 所需弦长 = %.6g m = %.6g km（装置尺度的 %.5g 倍）'
  % (lam_mfp, lam_mfp / 1.0e3, lam_mfp / L_chord))
p('练习6  屏蔽层吸收 %.6g MW / %.6g m^2 = %.6g kW/m^2 = %.6g MW/m^2'
  % (P_BR / 1.0e6, A_wall, flux_wall / 1.0e3, flux_wall / 1.0e6))
p('对照   P_LINE / P_alpha = %.6g%%' % (100.0 * P_LINE / P_alpha))
p('对照   P_rad / P_alpha   = %.6g%%' % (100.0 * P_rad / P_alpha))
p('对照   tau_rad / tau_CQ  = %.6g（tau_CQ = 0.0769074 s）' % (tau_rad / 0.0769074))
p('对照   P_CYC(光学薄) / P_fus = %.6g -> %.6g%%' % (P_CYC_thin / P_fus, 100.0 * P_CYC_thin / P_fus))

head('7. 印刷精度自检（正文与门禁只引用这些值）')
vals = [
    ('[1]  P_br', P_br, 'W/m^3'),
    ('[2]  P_BR', P_BR / 1.0e6, 'MW'),
    ('[3]  P_BR/P_alpha', 100.0 * P_BR / P_alpha, '%'),
    ('[4]  P_BR/P_fus', 100.0 * P_BR / P_fus, '%'),
    ('[5]  T_e/(m_e c^2)', T / (m_e * c * c), ''),
    ('[6]  相对论修正因子', rel, ''),
    ('[7]  修正后 P_BR', P_BR_rel / 1.0e6, 'MW'),
    ('[8]  平均光子波长', lam_ph * 1.0e9, 'nm'),
    ('[9]  f_ce', nu_ce / 1.0e9, 'GHz'),
    ('[10] u_B', u_B, 'J/m^3'),
    ('[11] P_1', P1, 'W'),
    ('[12] P_cyc(光学薄)', P_cyc_thin / 1.0e6, 'MW/m^3'),
    ('[13] P_CYC(光学薄)', P_CYC_thin / 1.0e6, 'MW'),
    ('[14] P_CYC(薄)/P_BR', P_CYC_thin / P_BR, ''),
    ('[15] P_CYC(薄)/P_fus', P_CYC_thin / P_fus, ''),
    ('[16] A_wall', A_wall, 'm^2'),
    ('[17] F_RJ', F_RJ, 'W/m^2'),
    ('[18] P_CYC(黑体)', P_RJ / 1.0e6, 'MW'),
    ('[19] 再吸收压低倍数', P_CYC_thin / P_RJ, ''),
    ('[20] P_CYC(黑体)/P_BR', 100.0 * P_RJ / P_BR, '%'),
    ('[21] alpha_ff', alpha_SI, 'm^-1'),
    ('[22] tau_ff', tau_ff, ''),
    ('[23] 平均自由程', lam_mfp, 'm'),
    ('[24] P_line', P_line / 1.0e3, 'kW/m^3'),
    ('[25] P_LINE', P_LINE / 1.0e6, 'MW'),
    ('[26] P_LINE/P_BR', 100.0 * P_LINE / P_BR, '%'),
    ('[27] f_W*', 100.0 * f_W_star, '%'),
    ('[28] P_rad', P_rad / 1.0e6, 'MW'),
    ('[29] P_rad/P_fus', 100.0 * P_rad / P_fus, '%'),
    ('[30] P_alpha/P_rad', P_alpha / P_rad, ''),
    ('[31] tau_rad', tau_rad, 's'),
    ('[32] 壁面轫致热流', flux_wall / 1.0e3, 'kW/m^2'),
    ('[33] 练习1 P_BR(Z=1.5)', p_Z15 / 1.0e6, 'MW'),
    ('[34] 练习2 P_br(20keV)', P_br2 / 1.0e3, 'kW/m^3'),
    ('[35] 练习2 涨幅', 100.0 * (math.sqrt(T2 / T_keV) - 1.0), '%'),
    ('[36] 练习3 n_W/n_e', 100.0 * n_W_br / n_e, '%'),
    ('[37] 练习4 亮度温度', T_bb, 'keV'),
    ('[38] 练习5 弦长', lam_mfp / 1.0e3, 'km'),
    ('[39] 练习6 壁面沉积', flux_wall / 1.0e6, 'MW/m^2'),
    ('[40] P_LINE/P_alpha', 100.0 * P_LINE / P_alpha, '%'),
    ('[41] P_rad/P_alpha', 100.0 * P_rad / P_alpha, '%'),
    ('[42] tau_rad/tau_CQ', tau_rad / 0.0769074, ''),
]
for name, v, u in vals:
    p('%-24s = %.6g %s' % (name, v, u))

dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'appX_calc.txt')
io.open(dst, 'w', encoding='utf-8', newline='').write('\n'.join(out) + '\n')
print('\n[写出] %s（%d 行）' % (dst, len(out)))
