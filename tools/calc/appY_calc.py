# -*- coding: utf-8 -*-
"""附录 Y（强场与激光等离子体）的唯一数值来源。

口径：
  * 工作点 D1（与附录 M/S/U/V/W/X 同一台装置的不同侧面）只在需要对照时出现；
    本附录的主角是一台「超强超短激光」，其参数全部在下面 DECLARED 里声明。
  * 所有打印值即为正文与门禁 P30 可以引用的值（4-6 位有效数字）。
  * 公式一律给出出处或「本书自拟」标注；文献值放 LIT 段，不进正文数值断言。
"""
import math

# ---------------- 常数（CODATA 2018） ----------------
c = 2.99792458e8          # m/s
e = 1.602176634e-19       # C
me = 9.1093837015e-31     # kg
mp = 1.67262192369e-27    # kg
eps0 = 8.8541878128e-12   # F/m
mu0 = 1.25663706212e-6    # H/m
hbar = 1.054571817e-34    # J s
keV = 1e3 * e
MeV = 1e6 * e
mec2 = me * c * c
meV = mec2 / e            # 电子的静止能量，eV

# ---------------- DECLARED：本附录声明的激光与靶参数 ----------------
lam = 0.8e-6              # 激光波长 0.8 um（Ti:sapphire）
I22 = 1e22                # W/cm^2 峰值强度
I = I22 * 1e4             # W/m^2
w0 = 5e-6                 # 焦斑半径 5 um
tau = 30e-15              # 脉冲宽度 30 fs
ne_gas = 1e18 * 1e6       # 气体射流/毛细管电子密度 1e18 cm^-3 -> m^-3
ne_solid = 1e23 * 1e6     # 固体密度全电离 1e23 cm^-3
rhoR_icf = 3.0            # g/cm^2 高增益面密度（声明）
rho_icf = 200.0           # g/cm^3 压缩后 DT 密度（声明）
Edrv = 2.0e6              # J 驱动器能量（NIF 量级，声明）
HB = 6.0                  # g/cm^2 燃耗常数（Atzeni & Meyer-ter-Vehn 式）
E_dt = 3.396e14           # J/kg DT 聚变能（17.6 MeV / 5 u）

DAO = '《等离子体物理自学教材》附录 Y 数值脚本'
print('=' * 72)
print(DAO, '—— 全部数字的唯一来源')
print('=' * 72)

# ---------------- Y.1 归一化振幅与有质动力势 ----------------
omega = 2 * math.pi * c / lam
k0 = omega / c
# a0 = sqrt(I lambda^2 / 1.37e18)  (I: W/cm^2, lambda: um, 线偏振)
a0 = math.sqrt(I22 * (lam * 1e6) ** 2 / 1.37e18)
E_L = a0 * me * omega * c / e                 # V/m
Up_eV = meV * (math.sqrt(1 + a0 ** 2 / 2) - 1)
gam_osc = math.sqrt(1 + a0 ** 2 / 2)
beta_osc = math.sqrt(1 - 1 / gam_osc ** 2)
E_S = me ** 2 * c ** 3 / (e * hbar)           # Schwinger 场 V/m
I_S = c * eps0 * E_S ** 2 / 2                 # W/m^2（口径：含 1/2）
print('\n[Y.1] 归一化振幅与有质动力势')
print('  激光波长 lam = %.4g m = %.1f um' % (lam, lam * 1e6))
print('  角频率 omega = %.6g rad/s = %.4g s^-1' % (omega, omega))
print('  波数 k0 = %.6g m^-1' % k0)
print('  归一化振幅 a0 = %.4f' % a0)
print('  激光电场 E_L = %.6g V/m = %.4g V/m' % (E_L, E_L))
print('  电子抖动 gamma = %.4f' % gam_osc)
print('  抖动速度 beta = v/c = %.6f' % beta_osc)
print('  有质动力势 U_p = %.4f MeV = %.5g eV' % (Up_eV / 1e6, Up_eV))
print('  相对论因子 sqrt(1+a0^2/2) = %.5f' % gam_osc)
print('  Schwinger 场 E_S = %.6g V/m' % E_S)
print('  Schwinger 强度 I_S = %.6g W/m^2 = %.4g W/cm^2' % (I_S, I_S / 1e4))
print('  E_L/E_S = %.6g' % (E_L / E_S))
print('  I/I_S = %.6g' % (I / I_S))
print('  E_L/E_a（原子场 5.14e11 V/m）= %.4f' % (E_L / 5.14e11))
print('  a0^2 = %.5g' % a0 ** 2)

# ---------------- Y.2 等离子体波：频率、波长、波破场 ----------------
wpe = math.sqrt(ne_gas * e ** 2 / (eps0 * me))
lam_p = 2 * math.pi * c / wpe
nc = eps0 * me * omega ** 2 / e ** 2
E_WB = me * c * wpe / e
print('\n[Y.2] 等离子体波')
print('  电子密度 n_e = %.4g m^-3 = %.4g cm^-3' % (ne_gas, ne_gas / 1e6))
print('  等离子体频率 omega_pe = %.6g rad/s' % wpe)
print('  等离子体波长 lam_p = %.6g m = %.5g um' % (lam_p, lam_p * 1e6))
print('  临界密度 n_c = %.6g m^-3 = %.5g cm^-3' % (nc, nc / 1e6))
print('  n_e/n_c = %.6g' % (ne_gas / nc))
print('  n_c/n_e = %.6g' % (nc / ne_gas))
print('  波破场 E_WB = %.6g V/m = %.4f GV/m' % (E_WB, E_WB / 1e9))
print('  E_WB / 1 GV/m = %.4f' % (E_WB / 1e9))
print('  E_WB/E_L = %.6g' % (E_WB / E_L))
print('  lam_p/lam = %.4f' % (lam_p / lam))

# ---------------- Y.3 尾场加速：三个长度与能量账 ----------------
vg_c = math.sqrt(1 - ne_gas / nc)
dv_c = ne_gas / (2 * nc)
L_d = lam_p * nc / ne_gas
dW = e * E_WB * L_d
P_las = I * math.pi * w0 ** 2
E_pulse = P_las * tau
Pc_GW = 17.4 * (nc / ne_gas)
P_over_Pc = P_las / (Pc_GW * 1e9)
A_spot = math.pi * w0 ** 2
print('\n[Y.3] 尾场加速')
print('  v_g/c = %.8f' % vg_c)
print('  (c-v_g)/c = n_e/(2n_c) = %.6g' % dv_c)
print('  失相长度 L_d = lam_p(n_c/n_e) = %.6g m = %.5g cm' % (L_d, L_d * 100))
print('  理想能量增益 dW = e E_WB L_d = %.6g J = %.5g eV = %.4f GeV'
      % (dW, dW / e, dW / e / 1e9))
print('  焦斑面积 A = %.6g m^2' % A_spot)
print('  激光功率 P = %.6g W = %.5g TW' % (P_las, P_las / 1e12))
print('  脉冲能量 E = P tau = %.6g J' % E_pulse)
print('  自聚焦临界功率 P_c = 17.4(n_c/n_e) GW = %.5g GW = %.5g TW'
      % (Pc_GW, Pc_GW / 1e3))
print('  P/P_c = %.5f' % P_over_Pc)
print('  平均加速梯度（理想）E_WB = %.4f GeV/m' % (E_WB / 1e9))
print('  1 GeV 需要的长度（E_WB 梯度）= %.5g m = %.4g cm'
      % (1e9 * e / (e * E_WB), 1e9 * e / (e * E_WB) * 100))

# ---------------- Y.4 相对论透明与自聚焦 ----------------
ratio_solid = ne_solid / nc
a0_tr = ratio_solid
I_tr_cm2 = a0_tr ** 2 * 1.37e18 / (lam * 1e6) ** 2
print('\n[Y.4] 相对论透明')
print('  固体密度 n_e = %.4g m^-3 = %.4g cm^-3' % (ne_solid, ne_solid / 1e6))
print('  n_e/n_c（固体）= %.5f' % ratio_solid)
print('  透明条件 a0 > n_e/n_c -> 需要的 a0 = %.5f' % a0_tr)
print('  本工作点 a0 = %.4f，余量倍数 = %.5f' % (a0, a0 / a0_tr))
print('  对应阈值强度 = %.6g W/cm^2 = %.4g W/cm^2' % (I_tr_cm2, I_tr_cm2))
print('  相对论临界密度倍数 gamma = %.4f（a0 处）' % gam_osc)
print('  n_c,rel = gamma n_c = %.6g m^-3' % (gam_osc * nc))

# ---------------- Y.5 辐射压强与强场 QED 参数 ----------------
P_rad = I / c
P_mag_tok = 1.11767e7     # Pa，附录 S/V 的磁压（对照用，取自附录 S 口径）
chi = gam_osc * E_L / E_S
print('\n[Y.5] 辐射压强与 QED 参数')
print('  辐射压强 I/c = %.6g Pa' % P_rad)
print('  = %.4g TPa = %.4g Mbar = %.4g Tbar' % (P_rad / 1e12, P_rad / 1e11, P_rad / 1e12))
print('  反射时 2I/c = %.6g Pa' % (2 * P_rad))
print('  托卡马克磁压（附录 S 口径）= %.6g Pa' % P_mag_tok)
print('  比值 (I/c)/P_mag = %.5g' % (P_rad / P_mag_tok))
print('  QED 参数 chi = gamma E/E_S = %.6g' % chi)
print('  1/chi = %.5g' % (1 / chi))
print('  chi/1 = %.4g %%' % (chi * 100))

# ---------------- Y.6 惯性约束：面密度与增益 ----------------
rhoR_SI = rhoR_icf * 10.0          # g/cm^2 -> kg/m^2
rho_SI = rho_icf * 1e3             # g/cm^3 -> kg/m^3
R_icf = rhoR_SI / rho_SI
V_icf = 4.0 / 3.0 * math.pi * R_icf ** 3
m_icf = V_icf * rho_SI
E_fus = m_icf * E_dt
f_burn = rhoR_icf / (rhoR_icf + HB)
E_burn = f_burn * E_fus
gain = E_burn / Edrv
print('\n[Y.6] 惯性约束')
print('  面密度 rhoR = %.4f g/cm^2 = %.5g kg/m^2' % (rhoR_icf, rhoR_SI))
print('  压缩密度 rho = %.4f g/cm^3 = %.5g kg/m^3' % (rho_icf, rho_SI))
print('  压缩半径 R = %.6g m = %.5g um' % (R_icf, R_icf * 1e6))
print('  燃料体积 V = %.6g m^3' % V_icf)
print('  燃料质量 m = %.6g kg = %.5g mg' % (m_icf, m_icf * 1e6))
print('  DT 比能 = %.6g J/kg' % E_dt)
print('  总聚变能 E_fus = %.6g J = %.5g MJ' % (E_fus, E_fus / 1e6))
print('  燃耗份额 f = rhoR/(rhoR+6) = %.6f = %.5g %%' % (f_burn, f_burn * 100))
print('  实际释放 E_burn = %.6g J = %.5g MJ' % (E_burn, E_burn / 1e6))
print('  驱动器能量 E_drv = %.6g J = %.4g MJ' % (Edrv, Edrv / 1e6))
print('  靶增益 G = E_burn/E_drv = %.5f' % gain)
print('  若驱动器效率 10%%，墙插增益 = %.4f' % (gain * 0.1))

# ---------------- Y.7 练习解答 ----------------
print('\n[Y.7] 练习解答')
# 练习 1：I = 1e20 W/cm^2 的 a0 与 U_p
a0_1 = math.sqrt(1e20 * (lam * 1e6) ** 2 / 1.37e18)
Up_1 = meV * (math.sqrt(1 + a0_1 ** 2 / 2) - 1) / 1e6
E_1 = a0_1 * me * omega * c / e
print('  1) a0 = %.5f, U_p = %.5f MeV, E_L = %.6g V/m = %.4g V/m'
      % (a0_1, Up_1, E_1, E_1))
print('     a0 比工作点小 %.4f 倍，U_p 小 %.4f 倍' % (a0 / a0_1, Up_eV / 1e6 / Up_1))
# 练习 2：n_e = 1e17 cm^-3 的 E_WB 与 lam_p
ne2 = 1e17 * 1e6
wpe2 = math.sqrt(ne2 * e ** 2 / (eps0 * me))
E_WB2 = me * c * wpe2 / e
lam_p2 = 2 * math.pi * c / wpe2
L_d2 = lam_p2 * (nc / ne2)
print('  2) n_e = 1e17 cm^-3: E_WB = %.5g V/m = %.4f GV/m, lam_p = %.5g m = %.4g um'
      % (E_WB2, E_WB2 / 1e9, lam_p2, lam_p2 * 1e6))
print('     L_d = %.5g m = %.4g cm' % (L_d2, L_d2 * 100))
print('     E_WB 随 sqrt(n_e) 变化：比 1e18 小 %.4f 倍' % (E_WB / E_WB2))
# 练习 3：CO2 激光（10.6 um）
lam3 = 10.6e-6
omega3 = 2 * math.pi * c / lam3
nc3 = eps0 * me * omega3 ** 2 / e ** 2
L_d3 = lam_p * (nc3 / ne_gas)
print('  3) CO2 激光 lam = 10.6 um: n_c = %.6g m^-3 = %.5g cm^-3' % (nc3, nc3 / 1e6))
print('     n_e/n_c = %.6g' % (ne_gas / nc3))
print('     L_d = %.6g m = %.5g mm（比 0.8 um 短 %.4f 倍）'
      % (L_d3, L_d3 * 1e3, L_d / L_d3))
# 练习 4：a0 = 1 需要的强度与脉冲能量
I_a1 = 1.37e18 / (lam * 1e6) ** 2
I_a1_SI = I_a1 * 1e4
P_a1 = I_a1_SI * math.pi * w0 ** 2
E_a1 = P_a1 * tau
print('  4) a0 = 1 需要 I = %.6g W/cm^2 = %.4g W/cm^2' % (I_a1, I_a1))
print('     P = %.6g W = %.5g TW, E = %.6g J = %.4g mJ'
      % (P_a1, P_a1 / 1e12, E_a1, E_a1 * 1e3))
print('     与工作点相比 I 小 %.4g 倍、能量小 %.5g 倍'
      % (I22 / I_a1, E_pulse / E_a1))
# 练习 5：相对论透明阈值强度（固体）
print('  5) 固体靶透明阈值 a0 = %.4f -> I = %.6g W/cm^2 = %.5g W/cm^2'
      % (a0_tr, I_tr_cm2, I_tr_cm2))
print('     比工作点低 %.5f 倍' % (I22 / I_tr_cm2))
print('     对应 P = %.6g W = %.5g TW'
      % (I_tr_cm2 * 1e4 * math.pi * w0 ** 2, I_tr_cm2 * 1e4 * math.pi * w0 ** 2 / 1e12))
# 练习 6：rhoR = 1 g/cm^2, rho = 100 g/cm^3
rhoR6, rho6 = 1.0, 100.0
R6 = rhoR6 * 10.0 / (rho6 * 1e3)
V6 = 4.0 / 3.0 * math.pi * R6 ** 3
m6 = V6 * rho6 * 1e3
E6 = m6 * E_dt
f6 = rhoR6 / (rhoR6 + HB)
print('  6) rhoR = 1 g/cm^2, rho = 100 g/cm^3: R = %.5g m = %.4g um' % (R6, R6 * 1e6))
print('     m = %.6g kg = %.5g mg, E_fus = %.6g J = %.5g MJ'
      % (m6, m6 * 1e6, E6, E6 / 1e6))
print('     f = %.6f = %.5g %%, E_burn = %.6g J = %.5g MJ'
      % (f6, f6 * 100, f6 * E6, f6 * E6 / 1e6))
print('     G = %.5f（E_drv = 2 MJ）' % (f6 * E6 / Edrv))

# ---------------- 对照与比例（正文引用） ----------------
print('\n[对照]')
print('  P_rad / P_mag(托卡马克) = %.5g' % (P_rad / P_mag_tok))
print('  E_WB / (1 GV/m) = %.4f' % (E_WB / 1e9))
print('  a0 / a0_tr = %.5f' % (a0 / a0_tr))
print('  L_d / lam_p = %.5g' % (L_d / lam_p))
print('  E_pulse / E_a1 = %.5g' % (E_pulse / E_a1))
print('  chi = %.6g（千分之 %.4g）' % (chi, chi * 1e3))
print('  I/I_S = %.6g' % (I / I_S))
print('  1/chi = %.5g' % (1 / chi))
# ---------------- 附加：单位换算与导出比值（正文引用） ----------------
print('\n[附加]')
print('  I/c = %.6g Pa = %.5g bar = %.5g Mbar'
      % (I / c, I / c / 1e5, I / c / 1e11))
print('  2I/c = %.6g Pa = %.5g Mbar' % (2 * I / c, 2 * I / c / 1e11))
print('  n_e 1e18 -> 1e17: E_WB 比 = %.5g, L_d 比 = %.5g, 增益 = %.5g GeV'
      % (math.sqrt(10.0), 10.0 ** 1.5, dW / e / 1e9 * 10.0))
print('\n完成：以上即正文可以引用的全部数值。')
