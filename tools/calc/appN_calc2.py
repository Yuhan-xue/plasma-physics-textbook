# -*- coding: utf-8 -*-
"""附录 N 数值复算 v2：修两处口径（年侵蚀按等离子体秒数；两点模型给出自洽解表）。"""
import math

m_u = 1.66053906660e-27
m_D = 2.013553212745 * m_u
E_EV = 1.602176634e-19


def out(label, val, unit=''):
    print('%-48s %-16s %s' % (label, val, unit))


def c_s(T_eV, m=m_D):
    return math.sqrt(T_eV * E_EV / m)


print('== A. 功率账与靶板几何（N.1）==')
P_fus, P_ext, f_rad_core = 300e6, 40e6, 0.30
P_alpha = P_fus / 5
P_sol = (P_alpha + P_ext) * (1 - f_rad_core)
R0, lam_q, f_x, share = 6.2, 5e-3, 5.0, 0.5
A_wet = 2 * math.pi * R0 * lam_q * f_x
q_par = P_sol * share / A_wet
out('P_alpha (MW)', P_alpha / 1e6)
out('P_SOL (MW)', P_sol / 1e6)
out('A_wet 单靶 (m^2)', A_wet)
out('q_par (MW/m^2)', q_par / 1e6)
for deg in (90, 10, 5, 3, 1):
    out('  q_perp @ %2d deg (MW/m^2)' % deg, q_par * math.sin(math.radians(deg)) / 1e6)
out('杠杆 1/sin3deg', 1 / math.sin(math.radians(3)))
wall = 4 * math.pi ** 2 * R0 * 2.0 * 2.0
out('A_wall (m^2)', wall)
out('70 MW 均匀落全壁 (kW/m^2)', P_sol / wall / 1e3)

print()
print('== B. 两点模型（N.2）：q_par 固定，扫靶板温度 ==')
kappa0, L_conn, gamma_sh = 2000.0, 40.0, 7.0
T_t_list = [5, 10, 15, 30, 50]
rows = []
for T_t in T_t_list:
    cs = c_s(T_t)
    n_t = q_par / (gamma_sh * (T_t * E_EV) * cs)
    T_u = (T_t ** 3.5 + 3.5 * q_par * L_conn / kappa0) ** (2 / 7.0)
    n_u = n_t * T_t / T_u
    p_u = n_u * T_u
    rows.append((T_t, cs, n_t, T_u, n_u, p_u))
    out('T_t=%2d eV' % T_t, 'c_s=%.0f n_t=%.3g T_u=%.1f n_u=%.3g n_uT_u=%.3g' % (cs, n_t, T_u, n_u, p_u))
out('低再循环对照：n_t=2e19, T_t=15 eV 的 q_par (MW/m^2)',
    gamma_sh * 2e19 * (15 * E_EV) * c_s(15) / 1e6)
out('Gamma_i = n_t c_s @ n_t=2e19,T=15 (m^-2 s^-1)', 2e19 * c_s(15))

print()
print('== C. 脱靶与杂质辐射（N.3）==')
P_target = P_sol * share
q_allow = 5e6
f_rad_need = 1 - q_allow / q_par
P_rad_div = f_rad_need * P_target
n_e_div, V_div = 1e20, 20.0
out('单靶入射功率 (MW)', P_target / 1e6)
out('需辐射掉的比例（靶面允许 5 MW/m^2）', f_rad_need)
out('需辐射功率 (MW)', P_rad_div / 1e6)
for Lz in (3e-32, 3e-33):
    n_imp = P_rad_div / (n_e_div * Lz * V_div)
    out('  Lz=%.0e' % Lz, 'n_imp=%.3g, 比例=%.3g%%' % (n_imp, 100 * n_imp / n_e_div))

print()
print('== D. 溅射与侵蚀（N.4）==')
n_W = 6.3e28
Gamma_i = 2e19 * c_s(15)
E_i = 5 * 15
out('离子轰击能量 E_i=5kT_t (eV)', E_i)
Y_steady, Y_elm = 1e-5, 1e-3
depth_rate = Y_steady * Gamma_i / n_W
out('稳态侵蚀速率 (nm/s)', depth_rate * 1e9)
plasma_s_per_year = 1e4 * 10.0
out('等离子体秒数/年 (1e4 次 x 10 s)', plasma_s_per_year)
out('稳态：um/yr（无再沉积）', depth_rate * plasma_s_per_year * 1e6)
out('稳态：um/yr（再沉积 90%）', depth_rate * plasma_s_per_year * 1e6 * 0.1)
out('稳态：5 mm 靶板寿命 (yr)', 5e-3 / (depth_rate * plasma_s_per_year * 0.1))
N_ion_elm = 1e6 / (300 * E_EV)
d_elm = Y_elm * N_ion_elm / n_W
out('单次 ELM 离子通量 (m^-2)', N_ion_elm)
out('单次 ELM 侵蚀深度 (nm)', d_elm * 1e9)
out('单次放电（1000 ELM）(um)', d_elm * 1000 * 1e6)
out('ELM 侵蚀（um/yr，无再沉积）', d_elm * 1000 * 1e4 * 1e6)
out('ELM 侵蚀（mm/yr，再沉积 90%）', d_elm * 1000 * 1e4 * 1e3 * 0.1)
out('ELM 主导下的 5 mm 靶板寿命 (yr)', 5e-3 / (d_elm * 1000 * 1e4 * 0.1))
out('ELM 稳态侵蚀之比', (d_elm * 1000 * 1e4) / (depth_rate * plasma_s_per_year))

print()
print('== E. ELM 瞬态热负荷（N.5）==')
W_th = 3 * 1e20 * (10e3 * E_EV) * 840.0
W_elm = 0.10 * W_th
q_elm = W_elm * 0.5 / 20.0 / 0.3e-3
out('W_th (MJ)', W_th / 1e6)
out('DeltaW_ELM=10% (MJ)', W_elm / 1e6)
out('峰值瞬态热流 (GW/m^2)', q_elm / 1e9)
out('相对稳态 q_perp(3deg)=1.88 (倍)', q_elm / (q_par * math.sin(math.radians(3))))
out('ELM 频率 100 Hz 时平均功率份额（40 MJ/次 / 400 MJ）', 40 / 400)

print()
print('== F. 氚（N.6）==')
E_fus = 17.6e6 * E_EV
for P in (500e6, 300e6):
    rate = P / E_fus
    g_day = rate * 3.0160492 / 6.02214076e23 * 86400
    out('P_fus=%.0f MW: 消耗 (1/s)' % (P / 1e6), rate)
    out('  g/day', g_day)
    out('  滞留 1%% -> g/day', g_day * 0.01)
    out('  1 kg 上限 / 1%% 滞留 (yr)', 1000 / (g_day * 0.01) / 365)
    out('  10 年压 1 kg 需滞留率', 1000 / (10 * 365) / g_day)
lam = math.log(2) / (12.32 * 365.25)
out('1 kg 每天衰变 (g/day)', 1000 * lam)
out('衰变/(消耗 76.7 g/day)', lam * 1000 / 76.727)
