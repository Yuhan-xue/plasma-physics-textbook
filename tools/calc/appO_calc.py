# -*- coding: utf-8 -*-
"""附录 O（电推进）全部数值的独立复算脚本。正文与门禁都以本脚本输出为准。"""
from __future__ import print_function
import math

g0 = 9.80665
e = 1.602176634e-19
u = 1.66053906660e-27
mXe = 131.293 * u
mAr = 39.948 * u
k_B = 1.380649e-23
m_e = 9.1093837015e-31
eps0 = 8.8541878128e-12


def show(name, val, unit='', fmt='%.6g'):
    print('%-46s %s %s' % (name, fmt % val, unit))


print('== 常数 ==')
show('mXe [kg]', mXe, '', '%.6e')
show('mAr [kg]', mAr, '', '%.6e')

print('\n== O.1 比冲与火箭方程 ==')
for isp in (450.0, 1800.0, 3000.0, 5000.0):
    show('Isp=%.0f s -> v_ex [m/s]' % isp, isp * g0, 'm/s', '%.6g')
dv = 10000.0
for isp in (450.0, 1800.0, 3000.0, 5000.0):
    ve = isp * g0
    mr = math.exp(dv / ve)
    show('dv=10 km/s, Isp=%.0f -> m0/mf' % isp, mr, '', '%.4f')
    show('  推进剂质量分数 [%%]', 100 * (1 - 1 / mr), '%', '%.3f')
# 1000 kg 卫星，Δv=3 km/s
m0 = 1000.0
for isp in (300.0, 3000.0):
    ve = isp * g0
    mp = m0 * (1 - math.exp(-3000.0 / ve))
    show('m0=1000 kg, dv=3 km/s, Isp=%.0f -> 推进剂 [kg]' % isp, mp, 'kg', '%.4g')

print('\n== O.2 推力与功率账 ==')
for T, isp, eta in ((0.1, 3000.0, 0.60), (10.0, 3000.0, 0.60), (0.092, 3100.0, 0.61)):
    ve = isp * g0
    mdot = T / ve
    pjet = 0.5 * T * ve
    pin = pjet / eta
    show('T=%.3g N, Isp=%.0f s -> mdot [kg/s]' % (T, isp), mdot, 'kg/s', '%.6g')
    show('   mdot [mg/s]', mdot * 1e6, 'mg/s', '%.4g')
    show('   P_jet [W]', pjet, 'W', '%.6g')
    show('   P_in (eta=%.2f) [W]' % eta, pin, 'W', '%.6g')
    show('   T/P [mN/kW]', 1e3 * T / (pin / 1e3), 'mN/kW', '%.4g')
T, isp = 10.0, 3000.0
mdot = T / (isp * g0)
show('T=10 N 年耗推进剂 [kg/yr]', mdot * 3.15576e7, 'kg', '%.4g')

print('\n== O.3 静电加速（氙离子） ==')
for V in (300.0, 1000.0, 1500.0):
    v = math.sqrt(2 * e * V / mXe)
    show('V=%.0f V -> v_i [m/s]' % V, v, 'm/s', '%.6g')
    show('   Isp [s]', v / g0, 's', '%.4g')
Ib, Vb = 2.0, 1000.0
v = math.sqrt(2 * e * Vb / mXe)
Tb = Ib * math.sqrt(2 * mXe * Vb / e)
show('I_b=2 A, V_b=1000 V -> T [N]', Tb, 'N', '%.6g')
show('   P_jet = I_b V_b [W]', Ib * Vb, 'W', '%.6g')
show('   P_in (eta=0.65) [W]', Ib * Vb / 0.65, 'W', '%.6g')
show('   mdot = I_b m/q [kg/s]', Ib * mXe / e, 'kg/s', '%.6g')
show('   mdot [mg/s]', Ib * mXe / e * 1e6, 'mg/s', '%.4g')
# Child-Langmuir 空间电荷限制（单电荷氙离子）
for V, d in ((1000.0, 1.5e-3), (1500.0, 1.5e-3)):
    j = 4 * eps0 / 9 * math.sqrt(2 * e / mXe) * V ** 1.5 / d ** 2
    show('CL: V=%.0f V, d=%.1f mm -> j [A/m^2]' % (V, d * 1e3), j, 'A/m^2', '%.6g')
    show('   面积以通过 3 A [cm^2]', 3.0 / j * 1e4, 'cm^2', '%.4g')

print('\n== O.4 霍尔推力器 ==')
Vd, L, B = 300.0, 25e-3, 0.02
E = Vd / L
show('E = V_d/L [V/m]', E, 'V/m', '%.6g')
vE = E / B
show('v_E = E/B [m/s]', vE, 'm/s', '%.6g')
vi = math.sqrt(2 * e * Vd / mXe)
show('v_i (300 V) [m/s]', vi, 'm/s', '%.6g')
show('Isp [s]', vi / g0, 's', '%.4g')
wce = e * B / m_e
wci = e * B / mXe
show('omega_ce [rad/s]', wce, 'rad/s', '%.6g')
show('omega_ci [rad/s]', wci, 'rad/s', '%.6g')
Te = 20.0
vth = math.sqrt(2 * e * Te / m_e)
show('v_th(e, 20 eV) [m/s]', vth, 'm/s', '%.6g')
nn, sig = 1e19, 2e-19
nu = nn * sig * vth
show('nu_en = n_n sigma v [1/s]', nu, '1/s', '%.6g')
show('omega_ce tau_e = omega_ce/nu', wce / nu, '', '%.4g')
tau_i = L / vi
show('tau_i = L/v_i [s]', tau_i, 's', '%.4g')
show('omega_ci tau_i', wci * tau_i, '', '%.4g')
for Ib2, Vd2, eta2 in ((5.0, 300.0, 0.50), (10.0, 300.0, 0.50)):
    T2 = Ib2 * math.sqrt(2 * mXe * Vd2 / e)
    pj = Ib2 * Vd2
    pin = pj / eta2
    show('I_b=%.0f A -> T [N]' % Ib2, T2, 'N', '%.6g')
    show('   P_jet [W] / P_in [W]', pj, 'W', '%.6g')
    show('   P_in [W]', pin, 'W', '%.6g')
    show('   T/P [mN/kW]', 1e3 * T2 / (pin / 1e3), 'mN/kW', '%.4g')
show('霍尔电流/束流 = v_E/v_i', vE / vi, '', '%.4g')

print('\n== O.5 磁喷嘴与 VASIMR ==')
for P, isp, eta in ((200e3, 5000.0, 0.70),):
    ve = isp * g0
    pj = eta * P
    T3 = 2 * pj / ve
    show('P=%.0f kW, Isp=%.0f s, eta=%.2f -> P_jet [kW]' % (P / 1e3, isp, eta), pj / 1e3, 'kW', '%.6g')
    show('   T [N]', T3, 'N', '%.6g')
    show('   mdot [mg/s]', T3 / ve * 1e6, 'mg/s', '%.4g')
show('P=200 kW 在 Isp=2000 s 时的推力 [N]', 2 * 0.70 * 200e3 / (2000 * g0), 'N', '%.4g')
show('   对应 mdot [mg/s]', (2 * 0.70 * 200e3 / (2000 * g0)) / (2000 * g0) * 1e6, 'mg/s', '%.4g')
show('P=200 kW 在 Isp=5000 s 时的推力 [N]', 2 * 0.70 * 200e3 / (5000 * g0), 'N', '%.4g')
show('   对应 mdot [mg/s]', (2 * 0.70 * 200e3 / (5000 * g0)) / (5000 * g0) * 1e6, 'mg/s', '%.4g')
show('P=200 kW 在 Isp=10000 s 时的推力 [N]', 2 * 0.70 * 200e3 / (10000 * g0), 'N', '%.4g')
show('   对应 mdot [mg/s]', (2 * 0.70 * 200e3 / (10000 * g0)) / (10000 * g0) * 1e6, 'mg/s', '%.4g')
show('氙电离能 12.13 eV 占 300 V 束能 [%%]', 100 * 12.13 / 300, '%', '%.4g')
show('氙电离能 12.13 eV 占 1000 V 束能 [%%]', 100 * 12.13 / 1000, '%', '%.4g')
show('1 N 所需功率（40.79 mN/kW）[kW]', 1000 / 40.79, 'kW', '%.4g')

print('\n== O.6 中和、CEX 与羽流 ==')
Ie, dV = 2.0, 15.0
show('中和器功耗 I_e dV [W]', Ie * dV, 'W', '%.4g')
show('   占束功率 2 kW 的比例 [%%]', 100 * Ie * dV / (Ie * 1000.0), '%', '%.4g')
Lp, sig_cex = 0.30, 5e-19
for nbg in (1e17, 1e18):
    P_cex = nbg * sig_cex * Lp
    show('n_n=%.0e -> CEX 概率' % nbg, P_cex, '', '%.4g')
    show('   束离子损失 [%%]', 100 * P_cex, '%', '%.3g')
show('羽流发散 15 deg 的推力损失 [%%]', 100 * (1 - math.cos(math.radians(15))), '%', '%.4g')
show('羽流发散 25 deg 的推力损失 [%%]', 100 * (1 - math.cos(math.radians(25))), '%', '%.4g')

print('\n== O.7 练习题的答案 ==')
show('练1: Isp=2500 s -> v_ex [m/s]', 2500 * g0, 'm/s', '%.6g')
show('练1: dv=6 km/s 的质量比', math.exp(6000 / (2500 * g0)), '', '%.4f')
show('练2: T=0.2 N, Isp=2500 s -> mdot [mg/s]', 0.2 / (2500 * g0) * 1e6, 'mg/s', '%.4g')
show('练2: P_jet [W]', 0.5 * 0.2 * 2500 * g0, 'W', '%.6g')
show('练2: eta=0.55 -> P_in [W]', 0.5 * 0.2 * 2500 * g0 / 0.55, 'W', '%.6g')
show('练3: V=1200 V Xe+ -> v [m/s]', math.sqrt(2 * e * 1200 / mXe), 'm/s', '%.6g')
show('练3: I_b=4 A -> T [N]', 4 * math.sqrt(2 * mXe * 1200 / e), 'N', '%.6g')
show('练3: P_jet [W]', 4 * 1200.0, 'W', '%.6g')
j = 4 * eps0 / 9 * math.sqrt(2 * e / mXe) * 1200.0 ** 1.5 / (1.0e-3) ** 2
show('练4: CL V=1200 V d=1 mm -> j [A/m^2]', j, 'A/m^2', '%.6g')
show('练4: 3 A 所需面积 [cm^2]', 3.0 / j * 1e4, 'cm^2', '%.4g')
vth2 = math.sqrt(2 * e * 25.0 / m_e)
nu2 = 5e18 * 3e-19 * vth2
show('练5: n_n=5e18 s=3e-19 Te=25 eV -> v_th', vth2, 'm/s', '%.6g')
show('练5: nu_en [1/s]', nu2, '1/s', '%.6g')
show('练5: omega_ce tau (B=0.015 T)', e * 0.015 / m_e / nu2, '', '%.4g')
show('练5: v_i(250 V) [m/s]', math.sqrt(2 * e * 250 / mXe), 'm/s', '%.6g')
show('练5: Isp [s]', math.sqrt(2 * e * 250 / mXe) / g0, 's', '%.4g')
show('练6: 300 kW, Isp=4000 s, eta=0.65 -> T [N]', 2 * 0.65 * 300e3 / (4000 * g0), 'N', '%.6g')
show('练6: 推进剂 1 t 可用时间 [天]', 1000.0 / (2 * 0.65 * 300e3 / (4000 * g0) / (4000 * g0)) / 86400.0, 'day', '%.4g')
print('\n完成')
