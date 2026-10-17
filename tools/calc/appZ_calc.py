# -*- coding: utf-8 -*-
"""附录 Z（聚变堆工程：从等离子体到电网）的唯一数值来源。

口径：
  * 聚变侧沿用 D1 工作点（与附录 M/N/S/U/V/W/X/Y 同一台装置）：
    P_fus = 556.5 MW、P_alpha = 111.3 MW、P_n = 445.2 MW（本书统一取 1/5 与 4/5，
    精确的 D-T 能量分配是 3.5/17.6 = 0.1989、14.1/17.6 = 0.8011）。
  * 工程侧的新参数全部在 DECLARED 里声明（包层倍增、热效率、再循环份额、TBR、
    屏蔽衰减长度、dpa 率、可用率、冷却剂参数），正文照抄并注明是量级值。
  * 所有打印值即为正文与门禁 P31 可以引用的值。
"""
import math

# ---------------- 常数 ----------------
MeV = 1.602176634e-13           # J
u = 1.66053906660e-27           # kg
m_T = 3.0160492 * u             # 氚原子质量 kg
NA = 6.02214076e23
yr = 365.25 * 24 * 3600.0       # s
g0 = 9.80665

# ---------------- DECLARED：D1 聚变侧 ----------------
P_fus = 556.5                   # MW
f_alpha = 0.2                   # 本书口径（精确 0.1989）
f_n = 0.8                       # 本书口径（精确 0.8011）
P_alpha = P_fus * f_alpha       # 111.3 MW
P_n = P_fus * f_n               # 445.2 MW
W_th = 371.4                    # MJ 热储能（附录 U/V/X 口径）
A_wall = 489.532                # m^2 第一壁面积（附录 X 口径）
q_n_limit = 2.0                 # MW/m^2 中子壁负载工程上限（声明）
E_fus = 17.6                    # MeV 每次 D-T 反应释放

# ---------------- DECLARED：工程侧 ----------------
M_blanket = 1.15                # 包层能量倍增（声明）
eta_th = 0.40                   # 热效率（声明）
f_recirc = 0.15                 # 再循环份额（声明：低温、泵、加热、控制）
P_aux = 50.0                    # MW 辅助加热功率（声明，与附录 U 口径同量级）
TBR = 1.10                      # 氚增殖比（声明）
TBR_min_eng = 1.05              # 工程下限（声明）
T_half_T = 12.32                # 年，氚半衰期（声明）
inv_T = 2.0                     # kg，氚库存（声明）
tau_proc = 4.0                  # h，氚处理循环时间（声明）
lambda_sh = 0.15                # m，屏蔽中子衰减长度（声明）
dpa_per_MWm2 = 5.0              # dpa/年 per (MW/m^2)（声明）
dpa_limit = 50.0                # dpa，包层结构材料寿命上限（声明）
avail = 0.75                    # 可用率（声明）
cp_cool = 5.2                   # kJ/(kg K) 氦冷却剂比热（声明）
dT_cool = 200.0                 # K 冷却剂温升（声明）
f_pump = 0.01                   # 泵功率占热功率份额（声明）
COP_cryo = 300.0                # 4 K 制冷的电功/冷功比（声明）

DAO = '《等离子体物理自学教材》附录 Z 数值脚本'
print('=' * 72)
print(DAO, '—— 全部数字的唯一来源')
print('=' * 72)

# ---------------- Z.1 三本账接起来 ----------------
P_th = M_blanket * P_fus                 # MW
P_gross = eta_th * P_th                  # MW
P_recirc = f_recirc * P_gross            # MW
P_net = P_gross - P_recirc               # MW
Q = P_fus / P_aux
Q_eng = P_net / P_aux
eta_net = P_net / P_th
print('\n[Z.1] 从聚变功率到电功率')
print('  P_fus = %.1f MW, f_alpha = %.2f, f_n = %.2f' % (P_fus, f_alpha, f_n))
print('  P_alpha = %.4g MW, P_n = %.4g MW' % (P_alpha, P_n))
print('  包层倍增 M = %.2f -> P_th = %.5g MW' % (M_blanket, P_th))
print('  热效率 %.2f -> P_gross = %.5g MW' % (eta_th, P_gross))
print('  再循环份额 %.2f -> P_recirc = %.5g MW' % (f_recirc, P_recirc))
print('  净电功率 P_net = %.5g MW' % P_net)
print('  净效率 eta_net = P_net/P_th = %.5g' % eta_net)
print('  P_aux = %.1f MW -> Q = %.5g, Q_eng = P_net/P_aux = %.5g' % (P_aux, Q, Q_eng))
print('  Q/Q_eng = %.5g' % (Q / Q_eng))
_e_a = 3.5 / E_fus
_e_n = 14.1 / E_fus
print('  精确分配（3.5/17.6 与 14.1/17.6）= %.6g 与 %.6g' % (_e_a, _e_n))
print('  本书 1/5 与 4/5 的相对偏差 = %.5g %% 与 %.5g %%'
      % ((f_alpha / _e_a - 1) * 100, (f_n / _e_n - 1) * 100))
print('  每 1 MW 聚变功率对应的净电功率 = %.5g MW' % (P_net / P_fus))

# ---------------- Z.2 中子：壁负载与损伤 ----------------
q_n = P_n / A_wall                       # MW/m^2
margin_q = q_n_limit / q_n
dpa_rate = dpa_per_MWm2 * q_n            # dpa/年
life_y = dpa_limit / dpa_rate            # 年
E_n = P_n * 1e6 * yr                     # J/年
print('\n[Z.2] 中子壁负载与材料损伤')
print('  第一壁面积 A_wall = %.6g m^2' % A_wall)
print('  中子壁负载 q_n = P_n/A_wall = %.5g MW/m^2' % q_n)
print('  工程上限 %.1f MW/m^2 -> 余量 %.5g 倍' % (q_n_limit, margin_q))
print('  dpa 率 = %.5g * %.5g = %.5g dpa/年' % (dpa_per_MWm2, q_n, dpa_rate))
print('  寿命上限 %.1f dpa -> %.5g 年' % (dpa_limit, life_y))
print('  中子功率年输运能量 = %.6g J = %.5g GJ' % (E_n, E_n / 1e9))

# ---------------- Z.3 氚：燃耗、自持与库存 ----------------
n_T_per_J = 1.0 / (E_fus * MeV)          # 每焦耳聚变能消耗的氚原子数
mdot_T = n_T_per_J * (P_fus * 1e6) * m_T # kg/s
burn_day = mdot_T * 86400 * 1e3          # g/天
burn_year = mdot_T * yr                  # kg/年
lam_T = math.log(2) / T_half_T           # /年
decay_year = lam_T * inv_T               # kg/年
decay_day = decay_year * 1e3 / 365.25    # g/天
TBR_break = 1 + decay_year / burn_year
excess = TBR * burn_year - burn_year - decay_year
throughput = inv_T / tau_proc            # kg/h
throughput_day = throughput * 24         # kg/天
# 每 GW 聚变功率的年燃耗（与文献口径对照）
burn_per_GW = burn_year / (P_fus / 1e3)
print('\n[Z.3] 氚循环')
print('  每次 D-T 反应 %.1f MeV -> 每焦耳 %.6g 个氚原子' % (E_fus, n_T_per_J))
print('  氚原子质量 = %.6g kg' % m_T)
print('  燃耗率 = %.6g kg/s = %.5g g/天 = %.5g kg/年' % (mdot_T, burn_day, burn_year))
print('  折算到每 GW 聚变功率 = %.5g kg/年' % burn_per_GW)
print('  半衰期 %.2f 年 -> 衰变常数 %.6g /年' % (T_half_T, lam_T))
print('  库存 %.1f kg 的自衰变 = %.5g kg/年 = %.5g g/天' % (inv_T, decay_year, decay_day))
print('  自持所需 TBR 下限 = 1 + 衰变/燃耗 = %.5g' % TBR_break)
print('  声明 TBR = %.2f -> 净增 = %.5g kg/年' % (TBR, excess))
print('  处理流量 = 库存/循环时间 = %.5g kg/h = %.5g kg/天' % (throughput, throughput_day))
print('  P_alpha/P_n = %.5g' % (P_alpha / P_n))

# ---------------- Z.4 屏蔽与超导磁体 ----------------
atten_target = 1e-6
x_shield = lambda_sh * math.log(1.0 / atten_target)   # m
P_n_leak = P_n * atten_target                          # MW
P_cryo = P_n_leak * 1e3 * COP_cryo                     # kW（MW->kW 后乘 COP）
print('\n[Z.4] 屏蔽与磁体热负荷')
print('  衰减长度 lambda = %.2f m' % lambda_sh)
print('  要衰减到 %.0e 需要厚度 x = %.5g m' % (atten_target, x_shield))
print('  漏到磁体的中子功率 = %.5g MW = %.5g W' % (P_n_leak, P_n_leak * 1e6))
print('  4 K 制冷电功（COP=%.0f）= %.5g kW' % (COP_cryo, P_cryo))
print('  制冷电功占净电功率 = %.5g %%' % (P_cryo / 1e3 / P_net * 100))

# ---------------- Z.5 热工与冷却剂 ----------------
mdot_cool = P_th * 1e6 / (cp_cool * 1e3 * dT_cool)     # kg/s
P_pump = f_pump * P_th                                  # MW
print('\n[Z.5] 热工与冷却剂')
print('  冷却剂比热 %.2f kJ/(kg K)、温升 %.1f K' % (cp_cool, dT_cool))
print('  质量流量 = %.5g kg/s' % mdot_cool)
print('  泵功率（%.0f%% 热功率）= %.5g MW' % (f_pump * 100, P_pump))
print('  泵功率占毛电功率 = %.5g %%' % (P_pump / P_gross * 100))

# ---------------- Z.6 可用率与年发电量 ----------------
hours = 8760.0
E_year = P_net * 1e6 * hours * 3600 * avail            # J
GWh = E_year / 3.6e12
E_year_85 = P_net * 1e6 * hours * 3600 * 0.85 / 3.6e12
downtime = hours * (1 - avail)
gain85 = E_year_85 / GWh - 1
coal = GWh * 1e6 * 0.34 / 1e9                          # 声明：0.34 kg 标煤/kWh -> 百万吨
print('\n[Z.6] 可用率与年发电量')
print('  可用率 %.2f -> 年发电量 = %.6g J = %.5g GWh' % (avail, E_year, GWh))
print('  年停机 = %.1f h' % downtime)
print('  可用率提到 0.85 -> %.5g GWh（增幅 %.5g %%）' % (E_year_85, gain85 * 100))
print('  折合标煤（0.34 kg/kWh）= %.5g 百万吨' % coal)
print('  容量因子 = 可用率 = %.2f' % avail)

# ---------------- DECLARED 一览（正文照抄的声明值全在这里） ----------------
print('\n[DECLARED] 本附录声明的工程参数一览')
print('  M = %.2f, eta_th = %.2f, f_recirc = %.2f, P_aux = %.1f MW'
      % (M_blanket, eta_th, f_recirc, P_aux))
print('  q_n_limit = %.1f MW/m^2, dpa 系数 = %.1f, dpa 上限 = %.1f'
      % (q_n_limit, dpa_per_MWm2, dpa_limit))
print('  TBR = %.2f, TBR 工程下限 = %.2f, 半衰期 = %.2f 年, 库存 = %.1f kg, 处理循环 = %.1f h'
      % (TBR, TBR_min_eng, T_half_T, inv_T, tau_proc))
print('  TBR 工程下限的增量 = %.2f（= %.2f - 1）' % (TBR_min_eng - 1, TBR_min_eng))
print('  lambda_sh = %.2f m, 可用率 = %.2f（对照 %.2f）, 比热 = %.1f kJ/(kg K), 温升 = %.1f K'
      % (lambda_sh, avail, 0.85, cp_cool, dT_cool))
print('  泵份额 = %.2f, 制冷 COP = %.1f, 年小时数 = %.1f h, 标煤系数 = %.2f kg/kWh'
      % (f_pump, COP_cryo, hours, 0.34))

# ---------------- Z.7 练习解答 ----------------
print('\n[Z.7] 练习解答')
# 1) M = 1.30
P_th1 = 1.30 * P_fus
Pg1 = eta_th * P_th1
Pa1 = f_recirc * Pg1
Pn1 = Pg1 - Pa1
print('  1) M=1.30: P_th = %.5g MW, P_gross = %.5g MW, P_net = %.5g MW（涨 %.5g 倍）'
      % (P_th1, Pg1, Pn1, Pn1 / P_net))
# 2) eta_th = 0.45
Pg2 = 0.45 * P_th
Pn2 = Pg2 * (1 - f_recirc)
print('  2) eta_th=0.45: P_gross = %.5g MW, P_net = %.5g MW（涨 %.5g %%）'
      % (Pg2, Pn2, (Pn2 / P_net - 1) * 100))
# 3) 壁面积加倍
q_n3 = P_n / (2 * A_wall)
print('  3) A_wall 加倍: q_n = %.5g MW/m^2, 余量 %.5g 倍' % (q_n3, q_n_limit / q_n3))
# 4) 衰减到 1e-8
x8 = lambda_sh * math.log(1e8)
print('  4) 衰减到 1e-8: x = %.5g m（比 1e-6 厚 %.5g m）' % (x8, x8 - x_shield))
# 5) TBR = 1.05
exc5 = 1.05 * burn_year - burn_year - decay_year
print('  5) TBR=1.05: 净增 = %.5g kg/年（TBR=1.10 时 %.5g kg/年）' % (exc5, excess))
# 6) 可用率 0.85 的年发电量已在 Z.6
print('  6) A=0.85: %.5g GWh，增幅 %.5g %%，折合标煤 %.5g 百万吨'
      % (E_year_85, gain85 * 100, E_year_85 * 1e6 * 0.34 / 1e9))

# ---------------- 对照与比例 ----------------
print('\n[对照]')
print('  P_net/P_fus = %.5g' % (P_net / P_fus))
print('  eta_net = %.5g' % eta_net)
print('  q_n/q_n_limit = %.5g' % (q_n / q_n_limit))
print('  burn_year/P_fus[GW] = %.5g' % burn_per_GW)
print('  TBR/TBR_break = %.5g' % (TBR / TBR_break))
print('  年发电量/年燃耗 = %.5g GWh per kg' % (GWh / burn_year))
print('\n完成：以上即正文可以引用的全部数值。')
