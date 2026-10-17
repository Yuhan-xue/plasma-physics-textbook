# -*- coding: utf-8 -*-
"""附录 AO 数值唯一来源：真空与壁处理（只 print，不写文件）。

用法：python tmp\\appAO_calc.py | Out-File -Encoding utf8 tmp\\appAO_calc.txt
"""
from __future__ import print_function

import math

# ---------- DECLARED ----------
kB = 1.380649e-23
e = 1.602176634e-19
NA = 6.02214076e23
u = 1.66053906660e-27
sigma_SB = 5.670374419e-8

V_vessel = 832.2          # m^3，与附录 M/S/U/V/W/X/Y/Z/AD/AE/AI/AJ/AK/AL/AM/AN 同点
A_wall = 489.532          # m^2，附录 X 口径
T_gas = 300.0             # K
P_div = 5.0               # Pa，附录 AJ 口径
P_vessel = 1.0e-6         # Pa，真空室本底目标
S_pump = 100.0            # m^3/s，附录 AJ 口径
d_duct, L_duct = 1.0, 1.0  # m，分子流导管（偏滤器短管）
q_unbaked = 1.0e-5        # Pa m^3/(s m^2)，未烘烤不锈钢放气率（声明）
q_baked = 1.0e-8          # Pa m^3/(s m^2)，350 °C 烘烤后（声明）
T_bake = 623.15           # K = 350 °C
eps_wall = 0.1            # 壁发射率（声明）
Q_leak_meas = 3.0e-9      # Pa m^3/s，实测总漏率（声明）
Q_leak_spec = 1.0e-8      # Pa m^3/s，允许上限（声明）
sens_leak = 1.0e-11       # Pa m^3/s，氦检漏灵敏度（声明）
monolayer = 1.0e19        # m^-2，单层水分子面密度（声明）
m_H2O = 18.015            # u
film_B = 100e-9           # m，硼化膜厚（声明）
rho_B = 2.34e3            # kg/m^3，硼密度
M_B = 10.811              # u
R_rec = 0.95              # 再循环系数（附录 AJ 口径）
tau_p = 2.5               # s，粒子约束时间（附录 AJ 口径）
A_dep = 100.0             # m^2，共沉积区面积（声明）
t_dep = 50e-6             # m，共沉积层厚（声明）
rho_dep = 1.8e3           # kg/m^3，共沉积层密度（声明）
DC_ratio = 0.4            # D/C 原子比（声明）
M_C = 12.011              # u
M_D = 2.014              # u
m_fuel_year = 31.191      # kg/年，氚燃耗（附录 Z/AJ/AK 口径）
m_fuel_flow = 162.72      # kg/年，加料氚流（附录 AJ 口径）
n_core = 1.0e20           # m^-3，芯部电子密度（全书工作点）
retention = 0.01          # 滞留份额（附录 AJ/AK 口径）
inv_limit = 2.0           # kg，库存上限（附录 Z/AJ/AK 口径）
S_He_ratio = 0.1          # 低温泵对氦的抽速比（声明）
N_fus = 1.97352e20        # /s，聚变反应率（附录 AJ 口径）
S_ext = 2.0591e21         # /s，外部加料率（附录 AJ 口径）


def p(tag, label, val, unit=''):
    print('%-16s %-46s %.6g %s' % (tag, label, val, unit))


print('=== [AO.1] 真空是什么单位：粒子数、通量与吞吐量 ===')
n_div = P_div / (kB * T_gas)
p('[AO.1]', '偏滤器压强 5 Pa 对应的粒子数密度', n_div, 'm^-3')
N_vessel = n_div * V_vessel
p('[AO.1]', '真空室内总粒子数', N_vessel, '')
Q_throughput = P_div * S_pump
p('[AO.1]', '吞吐量 P·S', Q_throughput, 'Pa m^3/s')
p('[AO.1]', '吞吐量换算 mbar L/s', Q_throughput * 10.0, 'mbar L/s')
N_pump = Q_throughput / (kB * T_gas)
p('[AO.1]', '每秒抽走的粒子数', N_pump, 's^-1')
sccm = 1e-6 / 60.0 * 2.6867811e25
p('[AO.1]', '1 sccm 对应的粒子流', sccm, 's^-1')
p('[AO.1]', '外部加料 2.0591e21 /s 相当于', S_ext / sccm, 'sccm')
p('[AO.1]', '本底 1e-6 Pa 对应密度', P_vessel / (kB * T_gas), 'm^-3')
p('[AO.1]', '偏滤器与本底压强之比', P_div / P_vessel, '')

print('=== [AO.2] 抽气：有效抽速由导管决定 ===')
C_duct = 12.1 * d_duct ** 3 / L_duct
p('[AO.2]', '分子流圆导管流导 C（直径 1 m、长 1 m）', C_duct, 'm^3/s')
S_eff = S_pump * C_duct / (S_pump + C_duct)
p('[AO.2]', '有效抽速 S_eff', S_eff, 'm^3/s')
p('[AO.2]', 'S_eff 相对泵抽速的份额', S_eff / S_pump, '')
p('[AO.2]', '损失倍数 S_pump/S_eff', S_pump / S_eff, '')
C_req = 9.0 * S_pump
p('[AO.2]', '损失不超过一成所需的流导', C_req, 'm^3/s')
p('[AO.2]', '对应单管直径（长 1 m）', (C_req * L_duct / 12.1) ** (1.0 / 3.0), 'm')
p('[AO.2]', '若只用 1 m 短管，需要并联根数', C_req / C_duct, '')
p('[AO.2]', '低温泵对氦的抽速（声明比 0.1）', S_pump * S_He_ratio, 'm^3/s')
p('[AO.2]', '氦分压与本底压强之比', 1.0 / S_He_ratio, '')

print('=== [AO.3] 漏率与放气：本底压强是「放气/抽速」 ===')
Q_out_unbaked = q_unbaked * A_wall
Q_out_baked = q_baked * A_wall
p('[AO.3]', '未烘烤放气量', Q_out_unbaked, 'Pa m^3/s')
p('[AO.3]', '350 度烘烤后放气量', Q_out_baked, 'Pa m^3/s')
P_unbaked = Q_out_unbaked / S_pump
P_baked = Q_out_baked / S_pump
p('[AO.3]', '未烘烤平衡压强', P_unbaked, 'Pa')
p('[AO.3]', '烘烤后平衡压强', P_baked, 'Pa')
p('[AO.3]', '烘烤后相对本底目标的余量', P_vessel / P_baked, '')
p('[AO.3]', '不烘烤时要把压强压到 1e-6 Pa 所需抽速', Q_out_unbaked / P_vessel, 'm^3/s')
p('[AO.3]', '实测漏率占放气的份额（烘烤后）', Q_leak_meas / Q_out_baked, '')
p('[AO.3]', '实测漏率贡献的压强', Q_leak_meas / S_pump, 'Pa')
p('[AO.3]', '检漏灵敏度相对实测漏率的余量', Q_leak_meas / sens_leak, '')
p('[AO.3]', '允许漏率上限相对实测值的余量', Q_leak_spec / Q_leak_meas, '')
N_water = monolayer * A_wall
p('[AO.3]', '单层水的分子数', N_water, '')
p('[AO.3]', '单层水的质量', N_water / NA * m_H2O, 'g')

print('=== [AO.4] 烘烤与壁处理：加热功率、硼化、放电清洗 ===')
P_rad = eps_wall * sigma_SB * A_wall * (T_bake ** 4 - T_gas ** 4)
p('[AO.4]', '烘烤辐射加热功率（黑体近似）', P_rad, 'W')
p('[AO.4]', '烘烤加热功率（含对流余量 2 倍）', 2 * P_rad, 'W')
p('[AO.4]', '烘烤时长 100 h 的能量', 2 * P_rad * 100 * 3600 / 3.6e9, 'GJ')
m_B = A_wall * film_B * rho_B
p('[AO.4]', '硼化膜所需硼质量', m_B * 1e3, 'g')
p('[AO.4]', '硼化膜原子面密度', film_B * rho_B / M_B / u, 'm^-2')
p('[AO.4]', '硼化膜相当于单原子层的层数', film_B * rho_B / M_B / u / 1e19, '')
p('[AO.4]', '放电清洗 1 h 抽走的粒子数（1e22 /s 通量）', 1e22 * 3600, '')
p('[AO.4]', '放电清洗 1 h 相当于多少层水膜', 1e22 * 3600 / N_water, '')

print('=== [AO.5] 壁条件与再循环：净加料只是损失的百分之几 ===')
N_core = n_core * V_vessel
p('[AO.5]', '芯部粒子数 N', N_core, '')
N_wall = N_core / tau_p
p('[AO.5]', '芯部粒子损失率 N/tau_p', N_wall, 's^-1')
Phi_rec = R_rec * N_wall
p('[AO.5]', '再循环返回壁的粒子流', Phi_rec, 's^-1')
Phi_net = (1 - R_rec) * N_wall
p('[AO.5]', '净流出（未再循环）粒子流', Phi_net, 's^-1')
p('[AO.5]', '外部加料与净流出的比值', S_ext / Phi_net, '')
p('[AO.5]', '再循环流是外部加料的倍数', Phi_rec / S_ext, '')
p('[AO.5]', '壁面粒子通量', N_wall / A_wall, 'm^-2 s^-1')
p('[AO.5]', '溅射产额 1e-4 时的溅射率', N_wall / A_wall * 1e-4, 'm^-2 s^-1')
p('[AO.5]', '壁面单层原子被打掉一遍所需时间', 1e19 / (N_wall / A_wall * 1e-4), 's')

print('=== [AO.6] 共沉积与氚滞留：库存账 ===')
V_dep = A_dep * t_dep
m_dep = V_dep * rho_dep
p('[AO.6]', '共沉积层体积', V_dep, 'm^3')
p('[AO.6]', '共沉积层质量', m_dep * 1e3, 'g')
n_C = m_dep / (M_C * u)
p('[AO.6]', '共沉积层碳原子数', n_C, '')
n_D = DC_ratio * n_C
p('[AO.6]', '共沉积层氘氚原子数', n_D, '')
m_D = n_D * M_D * u
p('[AO.6]', '共沉积层氘氚质量（按氘计）', m_D * 1e3, 'g')
p('[AO.6]', '占 2 kg 库存上限的份额', m_D / inv_limit, '')
p('[AO.6]', '年滞留量（加料流的 1%）', m_fuel_flow * retention, 'kg/年')
p('[AO.6]', '达到 2 kg 上限所需时间', inv_limit / (m_fuel_flow * retention), '年')
p('[AO.6]', '共沉积层清扫周期（若每周期清 100 g）', m_D * 1e3 / 100.0, '周期/年')

print('=== [AO.7] 练习解答 ===')
p('[AO.7-1]', 'S=200 m^3/s 时烘烤后压强', Q_out_baked / 200.0, 'Pa')
p('[AO.7-1]', '相对 1e-6 Pa 的余量', P_vessel / (Q_out_baked / 200.0), '')
p('[AO.7-2]', '导管直径 0.5 m 时流导', 12.1 * 0.5 ** 3 / L_duct, 'm^3/s')
p('[AO.7-2]', '对应有效抽速', S_pump * (12.1 * 0.125 / L_duct) / (S_pump + 12.1 * 0.125 / L_duct), 'm^3/s')
p('[AO.7-3]', '漏率 1e-8 Pa m^3/s 贡献的压强', 1e-8 / S_pump, 'Pa')
p('[AO.7-3]', '相对 1e-6 Pa 目标的份额', (1e-8 / S_pump) / P_vessel, '')
p('[AO.7-4]', '硼化膜厚 200 nm 时的硼质量', A_wall * 200e-9 * rho_B * 1e3, 'g')
p('[AO.7-5]', 'R=0.9 时的净流出', (1 - 0.9) * N_wall, 's^-1')
p('[AO.7-5]', 'R=0.9 时外部加料与净流出之比', S_ext / ((1 - 0.9) * N_wall), '')
p('[AO.7-6]', '共沉积面积 200 m^2 时的氘氚质量', n_C * 2 * DC_ratio * M_D * u * 1e3, 'g')
p('[AO.7-6]', '占库存上限份额', n_C * 2 * DC_ratio * M_D * u / inv_limit, '')

print('=== [AO.8] 派生比值 ===')
p('[AO.8]', '抽走的粒子流与聚变反应率之比', N_pump / N_fus, '')
p('[AO.8]', '氦产率（= 反应率）与抽出粒子流之比', N_fus / N_pump, '')
p('[AO.8]', '偏滤器压强下 1 m^3 的气体质量（按氘）', n_div * M_D * u / 1e3, 'kg')
p('[AO.8]', '整室 5 Pa 氘气质量', n_div * V_vessel * M_D * u * 1e3, 'g')
p('[AO.8]', '烘烤前后放气率之比', q_unbaked / q_baked, '')
p('[AO.8]', '烘烤前后平衡压强之比', P_unbaked / P_baked, '')


# ---- [AO.8] 派生比值（正文里写出的倍数与百分比）----
print("=== [AO.8b] 派生比值（供交叉核对） ===")
p_unbaked = P_unbaked
print("%-46s %10.4f" % ("未烘烤压强 / 本底目标", p_unbaked / 1e-6))
print("%-46s %10.2f %%" % ("共沉积氘氚 / 2 kg 上限", 100.0 * m_D / 2.0))
print("%-46s %10.2f %%" % ("面积加倍后 / 2 kg 上限", 100.0 * 2 * m_D / 2.0))
