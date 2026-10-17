# -*- coding: utf-8 -*-
"""附录 AZ（聚变堆系统集成与设计裕度）数值唯一来源：只打印，输出 tmp/appAZ_calc.txt。

工作点 D1 与附录 M/S/U/V/W/AX/AY/AZ 同口径；工程侧输入取自附录 Z/AR/AV/AP/AY。
"""
from __future__ import print_function

import io
import math
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'appAZ_calc.txt')
_buf = []


def p(fmt, *a):
    s = fmt % a if a else fmt
    _buf.append(s)
    print(s)


e = 1.602176634e-19
keV = 1e3 * e
mu0 = 4.0e-7 * math.pi
yr = 3.15576e7

# ---------- 工作点 D1 ----------
P_fus = 556.5e6
P_alpha = 111.325e6
P_n = 445.175e6
P_CD = 50.0e6
P_cryo = 13.4157e6
P_plant = 8.0e6
P_rad = 27.8322e6          # 附录 X：轫致 + 回旋 + 线辐射
eta_th = 0.40
n_e = 1.0e20
T_e = 9.194 * keV
V_pl = 832.2
A_wall = 489.532
R0 = 6.2
a_pl = 2.0
B0 = 5.3
I_p = 15.0e6               # A
f_bs = 0.30                # 自举电流份额（声明）
eta_CD = 0.30              # 电流驱动效率（声明，10^20 A/(W m^2) 口径）
P_heat = 50.0e6            # 辅助加热（= P_CD 口径）

# ---------- AZ.1 系统账 ----------
p('[AZ.1] 总热功率 = P_alpha + P_CD + P_heat_ohmic = %.6g MW',
  (P_alpha + P_CD) / 1e6)
p('[AZ.1] 聚变功率 P_fus = %.6g MW；中子功率 P_n = %.6g MW', P_fus / 1e6, P_n / 1e6)
p('[AZ.1] 辐射损失（附录 X）= %.6g MW = %.6g%% P_fus', P_rad / 1e6, P_rad / P_fus * 100.0)
p('[AZ.1] 加热与电流驱动的电耗 = %.6g MW（P_CD = %.6g MW，再循环口径）',
  P_CD / 1e6, P_CD / 1e6)
p('[AZ.1] 低温与厂用电 = %.6g + %.6g MW', P_cryo / 1e6, P_plant / 1e6)
p('[AZ.1] 环流功率口径：P_CD + P_cryo + P_plant = %.6g MW',
  (P_CD + P_cryo + P_plant) / 1e6)

# ---------- AZ.2 0-D 功率平衡 ----------
Q_pl = P_fus / P_CD
Q_eng = P_fus / (P_CD + P_cryo + P_plant)
p('[AZ.2] 等离子体增益 Q_pl = P_fus/P_CD = %.6g', Q_pl)
p('[AZ.2] 工程增益 Q_eng = P_fus/(P_CD+P_cryo+P_plant) = %.6g', Q_eng)
P_loss = P_alpha + P_CD
p('[AZ.2] 需要被输运带走的功率 = P_alpha + P_CD = %.6g MW', P_loss / 1e6)
p('[AZ.2] 三倍密度能量 3nkTV = %.6g MJ',
  3.0 * n_e * T_e * V_pl / 1e6)
tau_E = 3.0 * n_e * T_e * V_pl / P_loss
p('[AZ.2] 自洽约束下的能量约束时间 tau_E = %.6g s', tau_E)
p('[AZ.2] 若 tau_E 只有 %.6g s，则需要 P_CD = %.6g MW',
  tau_E / 2.0, 3.0 * n_e * T_e * V_pl / (tau_E / 2.0) / 1e6 - P_alpha / 1e6)
p('[AZ.2] 点火（P_CD=0）所需 tau_E = %.6g s', 3.0 * n_e * T_e * V_pl / P_alpha)

# ---------- AZ.3 自洽解与迭代 ----------
tau_IPB = 0.145 * (I_p / 1e6) ** 0.93 * B0 ** 0.15 * (P_loss / 1e6) ** (-0.69) \
    * (n_e / 1e20) ** 0.41 * (2.0 * a_pl) ** 0.19 * R0 ** 1.97 * 1.0
p('[AZ.3] IPB98(y,2) 量级估计 tau_E = %.6g s', tau_IPB)
p('[AZ.3] 与自洽要求的比值 tau_IPB/tau_E = %.6g', tau_IPB / tau_E)
p('[AZ.3] H 因子 = tau_E/tau_IPB = %.6g', tau_E / tau_IPB)
p('[AZ.3] 迭代残差按 10 倍/步下降：第 1/2/3 步 = %.6g / %.6g / %.6g',
  0.1, 0.01, 0.001)
p('[AZ.3] 收敛判据（相对残差 < 1e-3）需要约 %.6g 次迭代', math.log(1e3) / math.log(10.0))
p('[AZ.3] 若 tau_E 只有定标值的 0.8 倍，需要 P_CD = %.6g MW',
  (3.0 * n_e * T_e * V_pl / (0.8 * tau_IPB) - P_alpha) / 1e6)

# ---------- AZ.4 工作窗 ----------
n_G = I_p / 1e6 / (math.pi * a_pl ** 2)   # 1e20 m^-3
p('[AZ.4] Greenwald 密度极限 n_G = %.6g x10^20 m^-3 = %.6g m^-3', n_G, n_G * 1e20)
p('[AZ.4] 工作密度占比 n/n_G = %.6g = %.6g%%', n_e / 1e20 / n_G, n_e / 1e20 / n_G * 100.0)
p('[AZ.4] 密度裕度 = %.6g', n_G / (n_e / 1e20))
p_pa = n_e * T_e
p('[AZ.4] 等离子体压强 = %.6g Pa', p_pa)
beta = 2.0 * mu0 * p_pa / (B0 ** 2)
p('[AZ.4] beta = %.6g = %.6g%%', beta, beta * 100.0)
beta_N = beta * 100.0 * a_pl * B0 / (I_p / 1e6)
p('[AZ.4] beta_N = %.6g（%%，m，T，MA 口径）', beta_N)
beta_N_lim = 3.0
p('[AZ.4] beta_N 裕度 = %.6g（上限 %.6g）', beta_N_lim / beta_N, beta_N_lim)
q95 = 5.0 * a_pl ** 2 * B0 / (R0 * I_p / 1e6)
p('[AZ.4] q95 圆柱近似 = %.6g', q95)
q95_shape = q95 * 2.6
p('[AZ.4] 含成形因子 2.6 的 q95 = %.6g', q95_shape)
p('[AZ.4] q95 裕度 = %.6g（下限 2.0）', q95_shape / 2.0)
nTtau = n_e * T_e / keV * tau_E   # 1e20 keV s/m^3 口径
p('[AZ.4] 三重积 nTtau = %.6g x10^20 keV s/m^3', nTtau / 1e20 * 1e20 / 1e20)
p('[AZ.4] 三重积口径 n(1e20) T(keV) tau(s) = %.6g', (n_e / 1e20) * (T_e / keV) * tau_E)
p('[AZ.4] 点火判据 %.6g keV s m^-3（= %.6g x10^20）的裕度 = %.6g',
  5e21, 5.0, (n_e / 1e20) * (T_e / keV) * tau_E / 5.0)
p('[AZ.1] alpha 自加热功率 = %.6g MW（中子份额 %.6g%%、alpha 份额 %.6g%%）',
  P_alpha / 1e6, P_n / P_fus * 100.0, P_alpha / P_fus * 100.0)
p('[交叉] 附录 AY 的氦浓度上限 600 appm、氚燃烧率 89.7384 g/day')

# ---------- AZ.5 设计裕度 ----------
m_n = n_G / (n_e / 1e20)
m_beta = beta_N_lim / beta_N
m_q = q95_shape / 2.0
m_tau = tau_IPB / tau_E
m_dpa = 50.0 / (5.0 * 8.91268)
m_tbr = 1.05 / 1.01802
m_avail = 0.821918 / 0.780822
p('[AZ.5] 物理裕度：密度 %.6g、beta_N %.6g、q95 %.6g、约束（tau_IPB/tau_E）%.6g',
  m_n, m_beta, m_q, m_tau)
p('[AZ.5] 工程裕度：dpa %.6g、TBR %.6g、可用率 %.6g', m_dpa, m_tbr, m_avail)
m_all = (m_n, m_beta, m_q, m_tau, m_dpa, m_tbr, m_avail)
m_min = min(m_all)
p('[AZ.5] 最小裕度 = %.6g（瓶颈项）', m_min)
p('[AZ.5] 裕度乘积 = %.6g', m_n * m_beta * m_q * m_tau * m_dpa * m_tbr * m_avail)
p('[AZ.5] 若各项独立取 1.2 的工程裕度，7 项乘起来 = %.6g', 1.2 ** 7)
p('[AZ.5] 常用做法：把最小裕度留在 1.5 以上；本工作点有 %.6g 项低于 1.5',
  sum(1 for m in m_all if m < 1.5))
p('[AZ.5] 低于 1.5 的项：%s',
  '、'.join(nm for nm, m in zip(('密度', 'beta_N', 'q95', '约束', 'dpa', 'TBR', '可用率'), m_all)
            if m < 1.5))
p('[AZ.5] 对最小裕度取 1.2 倍安全系数后 = %.6g', m_min / 1.2)

# ---------- AZ.6 灵敏度 ----------
p('[AZ.6] 对数灵敏度 S = dln(y)/dln(x)（固定其余量）')
p('[AZ.6] P_fus 对 n：S_n = 2；对 T（<sigma v> ~ T^2）：S_T = 2')
p('[AZ.6] 固定 beta 时 n ~ B^2，故 P_fus ~ B^4：S_B = 4')
p('[AZ.6] 固定 n 时 V ~ R^3：S_R = 3')
p('[AZ.6] Q_pl 对 tau_E：S_tau = 1；对 P_CD：S_PCD = -1')
p('[AZ.6] 若 P_loss 提高 10%%，tau_IPB 按 -0.69 次幂变化 = %.6g 倍',
  1.1 ** (-0.69))
p('[AZ.6] 此时为保持同一 tau_E，H 因子需要提高 = %.6g 倍', 1.1 ** 0.69)
p('[AZ.6] tau_E 乘 1.1 时，维持同一 tau_E 所需的 P_CD 变化 = %.6g MW',
  (3.0 * n_e * T_e * V_pl / (1.1 * tau_E) - P_alpha) / 1e6)
p('[AZ.6] B 从 5.3 提到 5.8 T（+%.6g%%）时 P_fus 变化 = %.6g 倍',
  (5.8 / 5.3 - 1.0) * 100.0, (5.8 / 5.3) ** 4)
p('[AZ.6] R 从 6.2 提到 6.8 m 时 P_fus 变化 = %.6g 倍', (6.8 / 6.2) ** 3)
p('[AZ.6] 最值钱的两个旋钮：B（S=4）与 R（S=3）；最贵的是 P_CD（S=-1，且耗电）')

# ---------- AZ.7 集成检查清单 ----------
p('[AZ.7] 热流：P_sep = P_loss - P_rad = %.6g MW', (P_loss - P_rad) / 1e6)
p('[AZ.7] 若全部打到偏滤器、湿面积 %.6g m^2，峰值热流 = %.6g MW/m^2',
  3.0, (P_loss - P_rad) / 1e6 / 3.0)
p('[AZ.7] 材料：5 年 dpa = %.6g（上限 50，裕度 %.6g）', 5.0 * 8.91268, m_dpa)
p('[AZ.7] 氚：TBR 1.05 相对 1.01802 的裕度 = %.6g', m_tbr)
p('[AZ.7] 可用率：%.6g 相对 %.6g 的裕度 = %.6g', 0.821918, 0.780822, m_avail)
p('[AZ.7] 电：P_gross = eta_th x P_fus = %.6g MW；净电 = %.6g MW',
  eta_th * P_fus / 1e6, (eta_th * P_fus - P_CD - P_cryo - P_plant) / 1e6)
p('[AZ.7] 净电效率 = %.6g%%',
  (eta_th * P_fus - P_CD - P_cryo - P_plant) / P_fus * 100.0)
p('[AZ.7] 自举电流 = %.6g MA（份额 %.6g）；需要外部驱动 %.6g MA',
  I_p / 1e6 * f_bs, f_bs, I_p / 1e6 * (1.0 - f_bs))
p('[AZ.7] 电流驱动功率需求（eta_CD = %.6g）= %.6g MW', eta_CD,
  (I_p / 1e6 * (1.0 - f_bs)) / eta_CD * 1e6 / 1e6)
p('[AZ.7] 集成检查：7 项裕度全部 >= 1 ？ %s',
  all(m >= 1.0 for m in (m_n, m_beta, m_q, m_tau, m_dpa, m_tbr, m_avail)))
p('[交叉] 附录 AR 的净电效率 26.7701%%（含变换损耗 2.20873 MW）与本节的 %.6g%% 是同一口径的两种写法',
  (eta_th * P_fus - P_CD - P_cryo - P_plant) / P_fus * 100.0)
p('[交叉] 附录 AR 的 P_net = 148.976 MW、Q_eng = 2.02345；本节的 Q_eng（简化分母）= %.6g',
  Q_eng)
p('[交叉] 附录 X 的辐射损失 27.8322 MW、附录 AV 的 TBR 下限 1.01802、附录 AP 的可用率 0.780822/0.821918')
p('[交叉] 附录 AY 的 5 年剂量 44.5634 dpa 与 50 dpa 上限')

# ---------- 练习 ----------
p('')
p('[练习 1] 若 P_CD 翻倍到 100 MW：Q_pl = %.6g；Q_eng = %.6g',
  P_fus / 100e6, P_fus / (100e6 + P_cryo + P_plant))
p('[练习 1] Q_pl 下降 %.6g 倍', Q_pl / (P_fus / 100e6))
p('[练习 2] 若 n 提到 1.2e20：n/n_G = %.6g、beta_N = %.6g、密度裕度 = %.6g',
  1.2 / n_G, beta_N * 1.2, n_G / 1.2)
p('[练习 3] 若 tau_E 定标只有 2.0 s：需要的 P_CD = %.6g MW；Q_pl = %.6g',
  (3.0 * n_e * T_e * V_pl / 2.0 - P_alpha) / 1e6,
  P_fus / max(1.0, (3.0 * n_e * T_e * V_pl / 2.0 - P_alpha)))
p('[练习 4] 若 B 从 5.3 降到 5.0 T：P_fus 变化 = %.6g 倍；beta_N 变化 = %.6g 倍',
  (5.0 / 5.3) ** 4, (5.3 / 5.0) ** 2)
p('[练习 5] 若 TBR = 1.02（比 1.01802 高 %.6g%%）：裕度 = %.6g；五年氚盈余 = %.6g kg',
  (1.02 - 1.01802) / 1.01802 * 100.0, 1.02 / 1.01802, (1.02 - 1.01802) * 89.7384 * 5.0)
p('[AZ.7] 清单用量（工作点/限值）：密度 %.6g、比压 %.6g、稳定性 %.6g、约束 %.6g、材料 %.6g、氚 %.6g、可用率 %.6g',
  n_e / 1e20 / n_G, beta_N / beta_N_lim, 2.0 / q95_shape, tau_E / tau_IPB,
  5.0 * 8.91268 / 50.0, 1.01802 / 1.05, 0.780822 / 0.821918)
p('[练习 2] 1.2 倍密度时 beta_N 裕度 = %.6g', beta_N_lim / (beta_N * 1.2))
p('[练习 6] 三个裕度 1.8/1.5/1.2 的最小值 = %.6g', 1.2)
p('[练习 3] P_CD 72.5547 相对 50 的涨幅 = %.6g%%',
  (72.5547 / 50.0 - 1.0) * 100.0)
p('[练习 3] Q_pl 从 11.13 掉到 7.67007 的跌幅 = %.6g%%', (1.0 - 7.67007 / 11.13) * 100.0)
p('[练习 4] B 降 5.3->5.0 时 P_fus 跌幅 = %.6g%%', (1.0 - (5.0 / 5.3) ** 4) * 100.0)
p('[练习 5] TBR 1.02 相对 1.01802 的裕度超出 = %.6g%%', (1.02 / 1.01802 - 1.0) * 100.0)


with io.open(OUT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_buf) + '\n')
print('\n[写出] %s（%d 行）' % (os.path.basename(OUT), len(_buf)))
