# -*- coding: utf-8 -*-
"""附录 S（能量约束与输运）与门禁 P19 的数值唯一来源。

用法：python tmp/appS_calc.py > tmp/appS_calc.txt

口径（正文必须与之一致）：
  * 装置输入只有 (I_p, B, n19, P_tot, R, a, kappa, M)；tau_E 由 IPB98(y,2) 在 H=1 给出，
    体积平均温度 T_avg 由功率平衡反解（使 W/P = tau_E 恰好成立），因此每个装置都是自洽工作点。
  * V = 2*pi^2*R*a^2*kappa（椭圆截面近似；强三角变位形会高估 20%~30%）
  * W = 3 n T V；tau_E = W/P_tot
  * 输运系数约定：柱位形、抛物型温度剖面、均匀加热 -> tau_E = 3 a^2/(8 chi)，故 chi_eff = 3 a^2/(8 tau_E)
  * 玻姆：chi_B = T_e[eV]/(16 B[T])  [m^2/s]
  * 回旋玻姆：chi_gB = chi_B * (rho_i/a)，rho_i = sqrt(2 m_i T_i)/(e B)（即 v_thi = sqrt(2 T_i/m_i)）
  * 书中旧写法（v1.15.1 及以前，错）：0.0562 H I^0.93 B^0.15 n19^0.41 P^-0.58 a^1.97 R^1.39 M^0.19
"""
import math

e = 1.602176634e-19
u = 1.66053906660 * 1e-27
keV = 1000.0 * e
m_D = 2.013553212745 * u

DEV = [
    # 名称, I[MA], B[T], n19, P_tot[MW], R[m], a[m], kappa, M
    ('D1 ITER 量级', 15.0, 5.3, 10.1, 150.0, 6.20, 2.00, 1.7, 2.5),
    ('D2 JET 量级', 3.0, 3.0, 4.0, 15.0, 2.96, 1.25, 1.4, 2.0),
    ('D3 DIII-D 量级', 1.5, 2.1, 5.0, 10.0, 1.67, 0.67, 1.8, 2.0),
]


def ipb98y2(I, B, n19, P, R, a, kap, M, H=1.0):
    eps = a / R
    return (0.0562 * H * I ** 0.93 * B ** 0.15 * n19 ** 0.41 * P ** -0.69
            * R ** 1.97 * eps ** 0.58 * kap ** 0.78 * M ** 0.19)


def book_old(I, B, n19, P, R, a, kap, M, H=1.0):
    return (0.0562 * H * I ** 0.93 * B ** 0.15 * n19 ** 0.41 * P ** -0.58
            * a ** 1.97 * R ** 1.39 * M ** 0.19)


def volume(R, a, kap):
    return 2.0 * math.pi ** 2 * R * a ** 2 * kap


print('=== 常数 ===')
print('e = %.9e C, m_D = %.9e kg, 1 keV = %.6e J' % (e, m_D, keV))

print()
print('=== S.4 三个装置：IPB98(y,2)（H=1）与自洽工作点 ===')
print('%-16s %5s %5s %6s %7s %6s %6s %5s %5s %8s %8s %9s' % (
    '装置', 'I/MA', 'B/T', 'n19', 'P/MW', 'R/m', 'a/m', 'kap', 'M', 'tau_E/s', 'V/m^3', 'T_avg/keV'))
rows = []
for name, I, B, n19, P, R, a, kap, M in DEV:
    tau = ipb98y2(I, B, n19, P, R, a, kap, M)
    V = volume(R, a, kap)
    n = n19 * 1e19
    # W = tau*P = 3 n T V  ->  T
    T = tau * P * 1e6 / (3.0 * n * e * V) / 1000.0
    rows.append((name, I, B, n19, P, R, a, kap, M, tau, V, T))
    print('%-16s %5.1f %5.1f %6.1f %7.1f %6.2f %6.2f %5.1f %5.1f %8.4f %8.1f %9.3f' % (
        name, I, B, n19, P, R, a, kap, M, tau, V, T))

print()
print('=== S.4 旧写法（v1.15.1 及以前的公式卡）偏差 ===')
print('%-16s %8s %8s %8s' % ('装置', 'IPB98/s', '旧写法/s', '旧写法/IPB98'))
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    t_old = book_old(I, B, n19, P, R, a, kap, M)
    print('%-16s %8.4f %8.4f %8.3f' % (name, tau, t_old, t_old / tau))

print()
print('=== S.4 P 口径的影响（D1：只改 P） ===')
for P in (87.0, 150.0):
    print('P = %5.1f MW -> tau_E = %.4f s' % (P, ipb98y2(15.0, 5.3, 10.1, P, 6.2, 2.0, 1.7, 2.5)))

print()
print('=== S.5 储能、功率平衡与 H 因子 ===')
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    W = 3.0 * n19 * 1e19 * T * keV * V
    print('%-16s W = %7.1f MJ, tau = W/P = %6.4f s, H = tau/tau_IPB98 = %5.3f' % (
        name, W / 1e6, W / (P * 1e6), W / (P * 1e6) / tau))

print()
print('=== S.5 设计杠杆：tau_E 提高 30% 时各参数需要的倍数 ===')
for k, v in (('I_p', 0.93), ('B', 0.15), ('n19', 0.41), ('P', -0.69),
             ('R', 1.97), ('a', 0.58), ('kappa', 0.78), ('M', 0.19)):
    print('%-6s 幂次 %+5.2f -> x %.3f' % (k, v, 1.3 ** (1.0 / v)))

print()
print('=== S.2/S.3 有效输运系数与玻姆/回旋玻姆 ===')
print('%-16s %8s %9s %9s %9s %9s %9s' % (
    '装置', 'chi_eff', 'chi_B', 'chi_gB', 'rho_i/mm', 'chi/chi_B', 'chi/chi_gB'))
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    chi = 3.0 * a ** 2 / (8.0 * tau)
    chi_B = T * 1000.0 / (16.0 * B)
    rho = math.sqrt(2.0 * m_D * T * keV) / (e * B)
    chi_gB = chi_B * (rho / a)
    print('%-16s %8.3f %9.2f %9.4f %9.3f %9.4f %9.3f' % (
        name, chi, chi_B, chi_gB, rho * 1e3, chi / chi_B, chi / chi_gB))

print()
print('=== S.3 回旋玻姆与玻姆的比值 chi_B/chi_gB = a/rho_i ===')
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    rho = math.sqrt(2.0 * m_D * T * keV) / (e * B)
    print('%-16s rho_i/a = %.3e -> chi_B/chi_gB = %7.1f' % (name, rho / a, a / rho))

print()
print('=== S.6 台基、ELM 与偏滤器热负荷（D1） ===')
name, I, B, n19, P, R, a, kap, M, tau, V, T = rows[0]
W = 3.0 * n19 * 1e19 * T * keV * V
dW = 0.10 * W
A_wet, tau_elm = 10.0, 0.5e-3
print('W = %.1f MJ, dW/W = 10%% -> dW = %.1f MJ' % (W / 1e6, dW / 1e6))
print('湿面积 %.0f m^2、沉积 %.1f ms -> 峰值 %.2f GW/m^2' % (A_wet, tau_elm * 1e3, dW / (A_wet * tau_elm) / 1e9))
print('稳态 P_SOL/A = %.0f MW / %.0f m^2 = %.1f MW/m^2；峰值是稳态的 %.0f 倍' % (
    P, A_wet, P / A_wet, (dW / (A_wet * tau_elm)) / (P * 1e6 / A_wet)))
print('功率平衡要求 ELM 平均频率 f > P/dW = %.2f Hz（不能更慢，否则台基储能无界增长）' % (P * 1e6 / dW))
print('若把单次 dW 压到 W 的 0.5%%（%.1f MJ），频率需升到 f > %.1f Hz' % (
    0.005 * W / 1e6, P * 1e6 / (0.005 * W)))
n_ped, T_ped, f_ped = 7.0e19, 3.5, 0.20
W_ped = 3.0 * n_ped * T_ped * keV * f_ped * V
print('台基 n=%.1e m^-3、T=%.1f keV、占体积 %.0f%% -> W_ped = %.1f MJ（%.0f%% 的 W）' % (
    n_ped, T_ped, f_ped * 100, W_ped / 1e6, 100 * W_ped / W))

print()
print('=== S.7 Greenwald 密度极限 ===')
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    nG = I / (math.pi * a ** 2)
    print('%-16s n_G = %.3f x10^20 m^-3, n/n_G = %.3f' % (name, nG, n19 / 10.0 / nG))

print()
print('=== S.8 测地声模 f_GAM = c_s/(2 pi R)，c_s = sqrt(T_e/m_i) ===')
for name, I, B, n19, P, R, a, kap, M, tau, V, T in rows:
    for lab, Tk in (('芯部 %.2f keV' % T, T), ('边缘 1.00 keV', 1.0)):
        cs = math.sqrt(Tk * keV / m_D)
        f = cs / R / (2 * math.pi)
        print('%-16s %-14s c_s = %7.3e m/s, f_GAM = %6.2f kHz' % (name, lab, cs, f / 1e3))

print()
print('=== S.9 无量纲参数与三种外推（都以 D3 为起点、D1 为目标） ===')
D = [dict(name=r[0], I=r[1], B=r[2], n19=r[3], P=r[4], R=r[5], a=r[6],
          kap=r[7], M=r[8], tau=r[9], V=r[10], T=r[11]) for r in rows]
for d in D:
    rho = math.sqrt(2.0 * m_D * d['T'] * keV) / (e * d['B'])
    print('%-16s rho* = rho_i/a = %.3e, T_avg = %5.2f keV' % (d['name'], rho / d['a'], d['T']))
d1, d3 = D[0], D[2]
# chi_B ~ T/B  ->  tau_B ~ a^2 B / T
# chi_gB ~ T^1.5/(B^2 a)  ->  tau_gB ~ a^3 B^2 / T^1.5
t_bohm = d3['tau'] * (d1['a'] ** 2 * d1['B'] / d1['T']) / (d3['a'] ** 2 * d3['B'] / d3['T'])
t_gb = d3['tau'] * (d1['a'] ** 3 * d1['B'] ** 2 / d1['T'] ** 1.5) / (d3['a'] ** 3 * d3['B'] ** 2 / d3['T'] ** 1.5)
print('从 %s 外推到 %s（假定 beta、nu*、q 不变）：' % (d3['name'], d1['name']))
print('  玻姆型 tau ~ a^2 B/T     -> %.3f s，比 IPB98 低 %.2f 倍' % (t_bohm, d1['tau'] / t_bohm))
print('  回旋玻姆型 tau ~ a^3 B^2/T^1.5 -> %.3f s，比 IPB98 低 %.2f 倍' % (t_gb, d1['tau'] / t_gb))
print('  IPB98(y,2) 本身          -> %.3f s' % d1['tau'])

print()
print('=== 六道练习答案 ===')
name, I, B, n19, P, R, a, kap, M, tau, V, T = rows[0]
print('练1: P=150 MW -> tau_E = %.3f s；P=87 MW -> %.3f s；相差 %.2f 倍' % (
    tau, ipb98y2(I, B, n19, 87.0, R, a, kap, M), ipb98y2(I, B, n19, 87.0, R, a, kap, M) / tau))
chi = 3.0 * a ** 2 / (8.0 * tau)
chi_B = T * 1000.0 / (16.0 * B)
rho = math.sqrt(2.0 * m_D * T * keV) / (e * B)
chi_gB = chi_B * rho / a
print('练2: chi_eff = %.3f m^2/s；chi_B = %.1f、chi_gB = %.4f；chi/chi_B = %.4f、chi/chi_gB = %.2f' % (
    chi, chi_B, chi_gB, chi / chi_B, chi / chi_gB))
W = 3.0 * n19 * 1e19 * T * keV * V
dW = 0.10 * W
print('练3: W = %.1f MJ、dW = %.1f MJ、峰值 %.2f GW/m^2、f_max = %.2f Hz' % (
    W / 1e6, dW / 1e6, dW / (10.0 * 0.5e-3) / 1e9, P * 1e6 / dW))
nG3 = rows[2][1] / (math.pi * rows[2][6] ** 2)
print('练4: D3 的 n_G = %.3f x10^20 m^-3，n19=5 即 %.3f n_G' % (nG3, 5.0 / 10.0 / nG3))
print('练5: I x %.3f、B x %.3f、R x %.3f、a x %.3f' % (
    1.3 ** (1 / 0.93), 1.3 ** (1 / 0.15), 1.3 ** (1 / 1.97), 1.3 ** (1 / 0.58)))
print('练6: D1 芯部 f_GAM = %.2f kHz、边缘 f_GAM = %.2f kHz' % (
    math.sqrt(T * keV / m_D) / R / (2 * math.pi) / 1e3,
    math.sqrt(1.0 * keV / m_D) / R / (2 * math.pi) / 1e3))
print()
print('以上数字是正文与门禁的唯一来源')
