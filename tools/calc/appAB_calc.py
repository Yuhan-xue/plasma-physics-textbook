# -*- coding: utf-8 -*-
"""附录 AB（空间与天体等离子体）的唯一数值来源。

工作点（DECLARED，1 AU 典型太阳风 + 磁层对照）：
  太阳风：n = 5 cm^-3、T = 10 eV（电子与离子同温）、B = 5 nT、v = 400 km/s
  磁层：赤道偶极场 B0 = 3.1e-5 T（地表赤道）、R_E = 6371 km；尾区 L = 1e7 m
  辐射带：1 MeV 电子在 L=4 处的赤道磁场 B = 100 nT
  电离层：以阿尔芬电导为对照
所有打印值即正文与门禁可以引用的值。
"""
import math

# ---------------- 常数 ----------------
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6
RE = 6.371e6                     # 地球半径 m

# ---------------- DECLARED ----------------
n_sw = 5e6                       # m^-3（5 cm^-3）
T_sw = 10.0 * e                  # J（10 eV）
B_sw = 5e-9                      # T（5 nT）
v_sw = 4.0e5                     # m/s（400 km/s）
B0_eq = 3.1e-5                   # T（地表赤道偶极场）
L_tail = 1.0e7                   # m（磁尾特征长度）
lnL = 20.0                       # 库仑对数（声明）
B_rb = 100e-9                    # T（L=4 赤道场）
E_rb = 1.0e6 * e                 # J（1 MeV 电子）

DAO = '《等离子体物理自学教材》附录 AB 数值脚本'
print('=' * 72)
print(DAO, '—— 全部数字的唯一来源')
print('=' * 72)

# ---------------- AB.2 太阳风：为什么是超声速的 ----------------
wpe = math.sqrt(n_sw * e ** 2 / (eps0 * me))
wpi = math.sqrt(n_sw * e ** 2 / (eps0 * mp))
c_s = math.sqrt(2 * T_sw / mp)
rho_sw = n_sw * mp
vA = B_sw / math.sqrt(mu0 * rho_sw)
M_s = v_sw / c_s
M_A = v_sw / vA
p_sw = n_sw * T_sw
beta = 2 * mu0 * p_sw / B_sw ** 2
print('\n[AB.2] 太阳风')
print('  omega_pe = %.6g rad/s = %.4g Hz' % (wpe, wpe / (2 * math.pi)))
print('  omega_pi = %.6g rad/s = %.4g Hz' % (wpi, wpi / (2 * math.pi)))
print('  c_s = %.6g m/s = %.5g km/s' % (c_s, c_s / 1e3))
print('  v_A = %.6g m/s = %.5g km/s' % (vA, vA / 1e3))
print('  声速马赫数 M_s = %.5g' % M_s)
print('  阿尔芬马赫数 M_A = %.5g' % M_A)
print('  等离子体 beta = %.5g' % beta)
print('  质量密度 = %.6g kg/m^3' % rho_sw)

# ---------------- AB.3 弓激波与磁层顶 ----------------
P_ram = rho_sw * v_sw ** 2
P_dyn = 0.5 * P_ram
P_mag0 = B0_eq ** 2 / (2 * mu0)
r_mp = RE * (P_mag0 / P_ram) ** (1.0 / 6.0)
A_eff = math.pi * RE ** 2
# 注意：动压 rho v^2 是压力（Pa），要乘速度才是能流密度（W/m^2）；
# 迎风面能流用 1/2 rho v^3。
F_sw = 0.5 * rho_sw * v_sw ** 3
P_sw_total = F_sw * A_eff
print('\n[AB.3] 弓激波与磁层顶')
print('  太阳风动压 rho v^2 = %.6g Pa（动压系数 1/2 时为 %.6g Pa）' % (P_ram, P_dyn))
print('  地表赤道磁压 B^2/(2mu0) = %.6g Pa' % P_mag0)
print('  磁层顶平衡距离 r_mp = %.5g R_E' % (r_mp / RE))
print('  迎风面积 pi R_E^2 = %.6g m^2' % A_eff)
print('  太阳风对磁层的功率量级 = %.6g W = %.5g GW' % (P_sw_total, P_sw_total / 1e9))

# ---------------- AB.5 磁重联：Sweet--Parker 与 Petschek ----------------
eta_sp = 5.2e-5 * lnL / (10.0 ** 1.5)        # Spitzer 电阻率（Z=1, T=10 eV）
S_lund = mu0 * L_tail * vA / eta_sp
rate_sp = 1.0 / math.sqrt(S_lund)
rate_pk = math.pi / (8.0 * math.log(S_lund))
print('\n[AB.5] 磁重联')
print('  Spitzer 电阻率 eta = %.6g ohm m（lnL = %.1f, T = 10 eV）' % (eta_sp, lnL))
print('  Lundquist 数 S = mu0 L v_A/eta = %.6g' % S_lund)
print('  Sweet--Parker 重联率 S^-1/2 = %.6g' % rate_sp)
print('  Petschek 重联率 pi/(8 lnS) = %.6g' % rate_pk)
print('  两者之比 = %.5g' % (rate_pk / rate_sp))
print('  S-P 重联时间 = L/(rate_sp*v_A) = %.6g s = %.5g 天'
      % (L_tail / (rate_sp * vA), L_tail / (rate_sp * vA) / 86400))
print('  Petschek 重联时间 = %.6g s = %.5g s（对照）'
      % (L_tail / (rate_pk * vA), L_tail / (rate_pk * vA)))

# ---------------- AB.6 辐射带与阿尔芬电导 ----------------
p_rb = math.sqrt(E_rb ** 2 + 2 * E_rb * me * c ** 2) / c
rg_rb = p_rb / (e * B_rb)
Sigma_A = 1.0 / (mu0 * vA)
tau_gyro = 2 * math.pi * me / (e * B_sw)
print('\n[AB.6] 辐射带与电导')
print('  1 MeV 电子动量 = %.6g kg m/s' % p_rb)
print('  回旋半径 r_g = %.6g m = %.5g km' % (rg_rb, rg_rb / 1e3))
print('  太阳风电子的回旋周期 = %.6g s' % tau_gyro)
print('  阿尔芬电导 Sigma_A = 1/(mu0 v_A) = %.5g S' % Sigma_A)

# ---------------- 对照 ----------------
print('\n[对照]')
print('  v_sw/c = %.6g' % (v_sw / c))
print('  c_s/v_A = %.5g' % (c_s / vA))
print('  r_mp/RE = %.5g，观测约 10 R_E' % (r_mp / RE))
print('  S 的平方根 = %.6g' % math.sqrt(S_lund))
print('\n完成：以上即正文可以引用的全部数值。')
