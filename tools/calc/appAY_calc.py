# -*- coding: utf-8 -*-
"""附录 AY（聚变材料与辐照损伤）数值唯一来源：只打印，不写正文。

工作点 D1 与附录 M/S/U/V/W/AK/AP/AV/AX 同口径；额外声明材料与损伤输入。
输出：tmp/appAY_calc.txt
"""
from __future__ import print_function

import io
import math
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'appAY_calc.txt')
_buf = []


def p(fmt, *a):
    s = fmt % a if a else fmt
    _buf.append(s)
    print(s)


# ---------- 常量与换算 ----------
e = 1.602176634e-19
MeV = 1e6 * e
keV = 1e3 * e
u = 1.66053906660e-27
N_A = 6.02214076e23
yr = 3.15576e7          # s
day = 86400.0
hour = 3600.0

# ---------- 工作点 D1 ----------
P_fus = 556.5e6          # W
P_n = 445.175e6          # W（附录 AV 的 79.99545% 份额）
E_n = 14.068 * MeV       # J
A_wall = 489.532         # m^2
V_pl = 832.2             # m^3

# ---------- 声明输入（材料与损伤） ----------
sigma_d = 7.0e-25        # m^2 = 7000 b，14 MeV 中子在铁基材料里的 NRT 位移截面量级
E_d = 40.0               # eV，Fe 的离位阈能（NRT 口径）
T_op_C = 500.0           # 摄氏度，结构材料运行温度
He_per_dpa = 12.0        # appm/dpa，聚变中子谱的氦产生率量级
swell_inc = 10.0         # dpa，肿胀孵化期
swell_rate = 0.2         # %/dpa，RAFM 钢的肿胀速率
dpa_limit = 50.0         # dpa，设计剂量上限
dpa_limit_hi = 70.0      # dpa，上限口径
dsig_sat = 400.0e6       # Pa，辐照硬化饱和增量
d0_hard = 0.1            # dpa，硬化饱和尺度
dbtt_coef = 0.4          # K/MPa，DBTT 位移系数
DBTT_unirr_C = -60.0     # 摄氏度，未辐照 DBTT
Y_sput = 1.0e-3          # 溅射产额（W，量级）
Gamma_D = 1.0e23         # m^-2 s^-1，偏滤器氘离子通量
f_redep = 0.95           # 再沉积份额
t_W_mm = 5.0             # mm，钨装甲厚度
rho_W = 19300.0          # kg/m^3
M_W = 0.18384            # kg/mol
E_mod = 200.0e9          # Pa
alpha_th = 12.0e-6       # 1/K
nu_pois = 0.3
sigma_y = 550.0e6        # Pa
dT_thin = 100.0          # K，薄壁温差（工况 A）
dT_thick = 500.0         # K，工况 B
eps_f = 0.2              # Coffin-Manson 延性系数
c_CM = -0.5
n_pulse_yr = 3000.0      # 次/年
life_creep_h = 1.0e5     # h，100 MPa / 550 摄氏度下的蠕变寿命
rho_steel = 7850.0       # kg/m^3
M_steel = 0.055845       # kg/mol（Fe 基）
gamma_surf = 1.0         # J/m^2
r_bubble = 1.0e-9        # m

# ---------- AY.2 位移损伤 ----------
NWL = 0.909388e6         # W/m^2，中子壁负载（附录 AV/Z 口径）
Gamma_n = NWL / E_n
p('[AY.2] 中子壁负载 = %.6g MW/m^2 = %.6g W/m^2（附录 AV/Z 口径）', NWL / 1e6, NWL)
p('[AY.2] 由 P_n/A_wall 复算 = %.6g MW/m^2（四舍五入差 1 个单位）',
  P_n / A_wall / 1e6)
p('[AY.2] 14.068 MeV 中子能量 = %.6g J', E_n)
p('[AY.2] 中子通量 = %.6g n/(m^2 s)', Gamma_n)
dpa_s = sigma_d * Gamma_n
p('[AY.2] dpa 率 = sigma_d x Gamma = %.6g dpa/s', dpa_s)
dpa_fpy = dpa_s * yr
p('[AY.2] dpa/fpy = %.6g（1 fpy = %.6g s）', dpa_fpy, yr)
p('[AY.2] 单位壁负载的损伤 = %.6g dpa/fpy per MW/m^2', dpa_fpy / (NWL / 1e6))
p('[AY.2] NRT 位移截面 = %.6g m^2 = %.6g b（离位阈能 %.6g eV）',
  sigma_d, sigma_d * 1e28, E_d)
p('[AY.2] 到 %.6g dpa 的寿命 = %.6g fpy = %.6g 年 = %.6g s',
  dpa_limit, dpa_limit / dpa_fpy, dpa_limit / dpa_fpy, dpa_limit / dpa_s)
p('[AY.2] 到 %.6g dpa 的寿命 = %.6g fpy', dpa_limit_hi, dpa_limit_hi / dpa_fpy)
p('[AY.2] 5 年运行累积 = %.6g dpa（上限的 %.6g%%）',
  5.0 * dpa_fpy, 5.0 * dpa_fpy / dpa_limit * 100.0)
p('[AY.2] 每 1 dpa 对应通量 = %.6g n/m^2', 1.0 / sigma_d)

# ---------- AY.3 氦与氢 ----------
n_steel = rho_steel / M_steel * N_A
p('[AY.3] 钢的原子密度 = %.6g m^-3', n_steel)
he_appm50 = He_per_dpa * dpa_limit
p('[AY.3] %.6g dpa 后氦浓度 = %.6g appm = %.6g at.%%', dpa_limit, he_appm50, he_appm50 * 1e-4)
n_He = he_appm50 * 1e-6 * n_steel
p('[AY.3] 氦原子密度 = %.6g m^-3', n_He)
V_STP_per_m3 = n_He / N_A * 22.41397e-3
p('[AY.3] 每立方米钢里的氦（STP）= %.6g m^3 = %.6g L', V_STP_per_m3, V_STP_per_m3 * 1000.0)
p_Laplace = 2.0 * gamma_surf / r_bubble
p('[AY.3] r = %.6g nm 气泡的内压 = %.6g Pa = %.6g GPa',
  r_bubble * 1e9, p_Laplace, p_Laplace / 1e9)
p('[AY.3] 氦产生率 = %.6g appm/dpa；%.6g dpa 时 %.6g appm',
  He_per_dpa, 5.0 * dpa_fpy, He_per_dpa * 5.0 * dpa_fpy)

# ---------- AY.4 肿胀 ----------
swell50 = (dpa_limit - swell_inc) * swell_rate
p('[AY.4] 孵化 %.6g dpa、速率 %.6g %%/dpa -> %.6g dpa 时肿胀 = %.6g%%',
  swell_inc, swell_rate, dpa_limit, swell50)
p('[AY.4] 线性尺寸变化 = %.6g%%', ((1.0 + swell50 / 100.0) ** (1.0 / 3.0) - 1.0) * 100.0)
swell70 = (dpa_limit_hi - swell_inc) * swell_rate
p('[AY.4] %.6g dpa 时肿胀 = %.6g%%', dpa_limit_hi, swell70)
dpa_inc_equiv = swell_inc / dpa_fpy
p('[AY.4] 孵化期折合运行时间 = %.6g fpy', dpa_inc_equiv)
p('[AY.4] 肿胀 1%% 所需剂量 = %.6g dpa', swell_inc + 1.0 / swell_rate)
swell5y = max(0.0, 5.0 * dpa_fpy - swell_inc) * swell_rate
p('[AY.4] 5 年（%.6g dpa）时肿胀 = %.6g%%', 5.0 * dpa_fpy, swell5y)
p('[AY.4] 线性变化：5 年 = %.6g%%、70 dpa = %.6g%%、1%% 肿胀 = %.6g%%',
  ((1.0 + swell5y / 100.0) ** (1.0 / 3.0) - 1.0) * 100.0,
  ((1.0 + swell70 / 100.0) ** (1.0 / 3.0) - 1.0) * 100.0,
  ((1.01) ** (1.0 / 3.0) - 1.0) * 100.0)
p('[AY.4] 1%% 肿胀报警线：剂量 %.6g dpa = %.6g fpy', swell_inc + 1.0 / swell_rate,
  (swell_inc + 1.0 / swell_rate) / dpa_fpy)

# ---------- AY.5 硬化与脆化 ----------
dsig = dsig_sat * (1.0 - math.exp(-dpa_limit / d0_hard))
p('[AY.5] 辐照硬化增量 = %.6g Pa = %.6g MPa（饱和 %.6g MPa）',
  dsig, dsig / 1e6, dsig_sat / 1e6)
p('[AY.5] 1 dpa 时的硬化 = %.6g MPa', dsig_sat * (1.0 - math.exp(-1.0 / d0_hard)) / 1e6)
p('[AY.5] 0.1 dpa 时的硬化 = %.6g MPa',
  dsig_sat * (1.0 - math.exp(-0.1 / d0_hard)) / 1e6)
dDBTT = dbtt_coef * dsig / 1e6
p('[AY.5] DBTT 位移 = %.6g K（系数 %.6g K/MPa）', dDBTT, dbtt_coef)
DBTT_irr_C = DBTT_unirr_C + dDBTT
p('[AY.5] 辐照后 DBTT = %.6g 摄氏度', DBTT_irr_C)
p('[AY.5] 运行温度下限（DBTT + 100 K）= %.6g 摄氏度', DBTT_irr_C + 100.0)
p('[AY.5] 硬化后的屈服强度 = %.6g MPa', (sigma_y + dsig) / 1e6)
p('[AY.5] 屈服强度提高 = %.6g%%', dsig / sigma_y * 100.0)

# ---------- AY.6 溅射侵蚀与再沉积 ----------
R_gross = Y_sput * Gamma_D * M_W / (rho_W * N_A)
p('[AY.6] 总侵蚀速率 = %.6g m/s = %.6g mm/fpy = %.6g mm/年',
  R_gross, R_gross * yr * 1000.0, R_gross * yr * 1000.0)
R_net = (1.0 - f_redep) * R_gross
p('[AY.6] 净侵蚀速率 = %.6g mm/fpy（再沉积份额 %.6g）', R_net * yr * 1000.0, f_redep)
p('[AY.6] 再沉积速率 = %.6g mm/fpy', f_redep * R_gross * yr * 1000.0)
life_W = t_W_mm / (R_net * yr * 1000.0)
p('[AY.6] %.6g mm 钨装甲的净侵蚀寿命 = %.6g fpy', t_W_mm, life_W)
p('[AY.6] 每 1e26 ions/m^2 的侵蚀深度 = %.6g mm',
  Y_sput * 1e26 * M_W / (rho_W * N_A) * 1000.0)
p('[AY.6] 净/总 比 = %.6g；若再沉积 0.99 则寿命 = %.6g fpy',
  1.0 - f_redep, t_W_mm / (0.01 * R_gross * yr * 1000.0))

# ---------- AY.7 热疲劳与蠕变 ----------
sig_th = E_mod * alpha_th * dT_thin / (1.0 - nu_pois)
p('[AY.7] 薄壁温差 %.6g K 的热应力 = %.6g Pa = %.6g MPa',
  dT_thin, sig_th, sig_th / 1e6)
p('[AY.7] 与屈服强度之比 = %.6g', sig_th / sigma_y)
eps_th = alpha_th * dT_thin
eps_shakedown = 2.0 * sigma_y / E_mod
p('[AY.7] 热应变（100 K）= %.6g，弹性安定界限 = %.6g', eps_th, eps_shakedown)
p('[AY.7] 是否弹性安定 = %s', eps_th < eps_shakedown)
sig_th_b = E_mod * alpha_th * dT_thick / (1.0 - nu_pois)
eps_th_b = alpha_th * dT_thick
eps_pl_b = eps_th_b - eps_shakedown
N_f = 0.5 * (eps_pl_b / (2.0 * eps_f)) ** (1.0 / c_CM)
p('[AY.7] 工况 B（%.6g K）：热应力 %.6g MPa、热应变 %.6g、塑性应变 %.6g',
  dT_thick, sig_th_b / 1e6, eps_th_b, eps_pl_b)
p('[AY.7] 工况 B 的 Coffin-Manson 寿命 = %.6g 次循环', N_f)
n_5y = n_pulse_yr * 5.0
p('[AY.7] 5 年脉冲数 = %.6g 次（余量 %.6g 倍）', n_5y, N_f / n_5y)
t_mission_h = 5.0 * yr / hour
p('[AY.7] 5 年 = %.6g h；蠕变寿命 %.6g h；余量 = %.6g 倍',
  t_mission_h, life_creep_h, life_creep_h / t_mission_h)
p('[AY.7] 运行温度窗口 = %.6g--%.6g 摄氏度（宽 %.6g K），运行点 %.6g 摄氏度',
  DBTT_irr_C + 100.0, 550.0, 550.0 - (DBTT_irr_C + 100.0), T_op_C)
p('[AY.7] 运行点距下限 = %.6g K、距上限 = %.6g K',
  T_op_C - (DBTT_irr_C + 100.0), 550.0 - T_op_C)

# ---------- 练习 ----------
p('')
p('[练习 1] NWL = 2.0 MW/m^2 -> %.6g dpa/fpy；到 50 dpa 需 %.6g fpy',
  2.0 * dpa_fpy / (P_n / A_wall / 1e6), dpa_limit / (2.0 * dpa_fpy / (P_n / A_wall / 1e6)))
p('[练习 1] 与基准之比 = %.6g', 2.0 / (P_n / A_wall / 1e6))
he500 = 10.0 * 50.0
p('[练习 2] 10 appm/dpa x 50 dpa = %.6g appm；每立方米 %.6g L（STP）',
  he500, he500 * 1e-6 * n_steel / N_A * 22.41397e-3 * 1000.0)
p('[练习 2] 每立方米氦原子数 = %.6g m^-3', he500 * 1e-6 * n_steel)
sw60 = (60.0 - 15.0) * 0.3
p('[练习 3] 孵化 15 dpa、速率 0.3 %%/dpa、60 dpa -> 肿胀 %.6g%%', sw60)
p('[练习 3] 线性尺寸变化 = %.6g%%', ((1.0 + sw60 / 100.0) ** (1.0 / 3.0) - 1.0) * 100.0)
R_gross2 = 2.0e-3 * Gamma_D * M_W / (rho_W * N_A)
p('[练习 4] Y = 2e-3 -> 总侵蚀 %.6g mm/fpy；再沉积 0.97 时净 %.6g mm/fpy',
  R_gross2 * yr * 1000.0, 0.03 * R_gross2 * yr * 1000.0)
p('[练习 4] 3 mm 装甲寿命 = %.6g fpy', 3.0 / (0.03 * R_gross2 * yr * 1000.0))
sig_th2 = E_mod * alpha_th * 200.0 / (1.0 - nu_pois)
p('[练习 5] 200 K：热应力 %.6g MPa、热应变 %.6g、弹性安定界限 %.6g',
  sig_th2 / 1e6, alpha_th * 200.0, eps_shakedown)
p('[练习 5] 是否超过屈服 = %s；是否弹性安定 = %s',
  sig_th2 > sigma_y, alpha_th * 200.0 < eps_shakedown)
dpa_5 = 5.0e-25 * Gamma_n * yr
p('[练习 6] sigma_d = 5e-25 m^2 -> %.6g dpa/fpy（为基准的 %.6g 倍）',
  dpa_5, dpa_5 / dpa_fpy)
p('[练习 6] 到 50 dpa 需 %.6g fpy', 50.0 / dpa_5)

# ---------- 与相邻附录的交叉口径 ----------
p('')
# ---------- 声明输入一览（供交叉核对） ----------
p('')
p('[输入] sigma_d = %.6g m^2、E_d = %.6g eV、He 产生率 = %.6g appm/dpa',
  sigma_d, E_d, He_per_dpa)
p('[输入] 肿胀孵化 = %.6g dpa、速率 = %.6g %%/dpa、上限 = %.6g / %.6g dpa',
  swell_inc, swell_rate, dpa_limit, dpa_limit_hi)
p('[输入] 硬化饱和 = %.6g MPa、饱和尺度 = %.6g dpa、DBTT 系数 = %.6g K/MPa、未辐照 DBTT = %.6g 摄氏度',
  dsig_sat / 1e6, d0_hard, dbtt_coef, DBTT_unirr_C)
p('[输入] 溅射产额 = %.6g、氘通量 = %.6g m^-2 s^-1、再沉积 = %.6g、钨装甲 = %.6g mm',
  Y_sput, Gamma_D, f_redep, t_W_mm)
p('[输入] rho_W = %.6g kg/m^3、M_W = %.6g kg/mol、rho_steel = %.6g kg/m^3、M_steel = %.6g kg/mol',
  rho_W, M_W, rho_steel, M_steel)
p('[输入] E = %.6g Pa、alpha = %.6g /K、nu = %.6g、sigma_y = %.6g Pa',
  E_mod, alpha_th, nu_pois, sigma_y)
p('[输入] 工况 A 温差 = %.6g K、工况 B 温差 = %.6g K、eps_f = %.6g、c = %.6g',
  dT_thin, dT_thick, eps_f, c_CM)
p('[输入] 热循环 = %.6g 次/年、蠕变寿命 = %.6g h、上限温度 = %.6g 摄氏度、运行点 = %.6g 摄氏度',
  n_pulse_yr, life_creep_h, 550.0, T_op_C)
p('[输入] 摩尔气体体积 = %.6g m^3/mol、N_A = %.6g /mol、gamma_surf = %.6g J/m^2、r_bubble = %.6g m',
  22.41397e-3, N_A, gamma_surf, r_bubble)
p('[AY.8] 5 年 = %.6g h；工况 B 的 N_f = %.6g 次', t_mission_h, N_f)
p('[AY.5] 0.1 dpa 达到饱和值的百分比 = %.6g%%',
  (1.0 - math.exp(-0.1 / d0_hard)) * 100.0)
p('[AY.7] 100 K 热应力占屈服的百分比 = %.6g%%', sig_th / sigma_y * 100.0)
p('[练习 5] 685.714 MPa / 550 MPa = %.6g', sig_th2 / sigma_y)
p('[交叉] 附录 AP 的 50.001 dpa 口径与本附录的 %.6g dpa 上限一致', dpa_limit)
p('[交叉] 附录 AV 的中子壁负载 %.6g MW/m^2、P_n = %.6g MW', P_n / A_wall / 1e6, P_n / 1e6)
p('[交叉] 附录 Z 的包层体积 %.6g m^3 与壁面积 %.6g m^2', 195.813, A_wall)
p('[交叉] 1 fpy = %.6g s；5 fpy = %.6g s', yr, 5.0 * yr)

with io.open(OUT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_buf) + '\n')
print('\n[写出] %s（%d 行）' % (os.path.basename(OUT), len(_buf)))
