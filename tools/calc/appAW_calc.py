# -*- coding: utf-8 -*-
"""附录 AW《高能量密度物理：从温稠密物质到实验室天体物理》数值唯一来源（只 print，不写文件）。

口径：
  - 两个声明工况，都取 3 倍频 Nd 激光（lambda = 0.351 um，光子能量 3.53 eV）：
      工况 N（纳秒驱动，接附录 AI 的间接驱动）：E_L = 10.0 kJ、tau = 3.0 ns、光斑半径 500 um；
      工况 F（飞秒强场，接附录 Y 的强场）：E_L = 10.0 J、tau = 30 fs、光斑半径 5 um；
  - 温稠密物质（WDM）声明：n_e = 1e29 m^-3、Zbar = 3.5、T_e = 10 eV、T_i = 10 eV、平均离子质量 6.5 u；
  - 冲击声明：D = 50 km/s、rho_0 = 1000 kg/m^3、gamma = 5/3（强冲击压缩比 (gamma+1)/(gamma-1)）；
  - 烧蚀压强用经典标度律 P[Mbar] = 8.6 (I_14)^{2/3} lambda_um^{-2/3}（直接驱动标度，见来源节）；
    黑腔反照率 alpha = 0.8、激光-X 射线转换效率 eta_X = 0.70；
  - 实验室天体声明：sigma = 1e6 S/m、v = 100 km/s、L = 1 mm、B = 10 T、P = 1 Mbar；
  - 常数：m_e、e、c、hbar、k_B、eps0、mu0、sigma_SB 取 CODATA-2018。

运行：python tmp\\appAW_calc.py > tmp\\appAW_calc.txt（PowerShell 用 Out-File -Encoding utf8）
"""
from __future__ import print_function

import math

# ---------------- 常数 ----------------
E_CHG = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 2.99792458e8
HBAR = 1.054571817e-34
K_B = 1.380649e-23
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
SIGMA_SB = 5.670374419e-8
U_KG = 1.66053906660e-27
EV_K = 11604.51812          # 1 eV 对应的开尔文
MBAR = 1.0e11               # Pa

# ---------------- 声明输入 ----------------
LAM_UM = 0.351              # 激光波长（um）
E_PHOT_EV = 3.53            # 光子能量（eV）

E_L_N = 10.0e3              # 工况 N：10 kJ
TAU_N = 3.0e-9              # 3 ns
R_N = 500e-6                # 500 um

E_L_F = 10.0                # 工况 F：10 J
TAU_F = 30e-15              # 30 fs
R_F = 5e-6                  # 5 um

# WDM 声明
NE_WDM = 1.0e29             # m^-3
ZBAR = 3.5
TE_WDM_EV = 10.0
TI_WDM_EV = 10.0
M_ION_U = 6.5

# 冲击声明
D_SHOCK = 50e3              # m/s
RHO0 = 1000.0               # kg/m^3
GAMMA = 5.0 / 3.0

# 辐射声明
ETA_X = 0.70
ALPHA_WALL = 0.80

# 实验室天体声明
SIGMA_COND = 1.0e6          # S/m
V_ASTRO = 100e3             # m/s
L_ASTRO = 1.0e-3            # m
B_ASTRO = 10.0              # T
P_ASTRO = 1.0 * MBAR        # Pa

# 强场声明（势垒抑制用氢的 13.6 eV）
E_ION_H_EV = 13.6
E_ION_H_J = E_ION_H_EV * E_CHG


def ev(joule):
    return joule / E_CHG


def kelvin(ev_):
    return ev_ * EV_K


def main():
    print('=' * 72)
    print('附录 AW《高能量密度物理：从温稠密物质到实验室天体物理》数值基线')
    print('激光 lambda = %.6g um，光子能量 %.6g eV' % (LAM_UM, E_PHOT_EV))
    print('=' * 72)

    # ================= AW.1 三个刻度 =================
    print('\n[AW.1] 三个刻度：功率密度、压强、温度')
    area_n = math.pi * R_N ** 2
    area_f = math.pi * R_F ** 2
    p_las_n = E_L_N / TAU_N
    i_n = p_las_n / area_n                     # W/m^2
    i_f = (E_L_F / TAU_F) / area_f
    i_n_cm2 = i_n / 1.0e4
    i_f_cm2 = i_f / 1.0e4
    print('工况 N：光斑面积 = %.6g m^2' % area_n)
    print('工况 N：激光功率 = %.6g W' % p_las_n)
    print('工况 N：强度 I = %.6g W/m^2 = %.6g W/cm^2' % (i_n, i_n_cm2))
    print('工况 F：光斑面积 = %.6g m^2' % area_f)
    print('工况 F：激光功率 = %.6g W' % (E_L_F / TAU_F))
    print('工况 F：强度 I = %.6g W/m^2 = %.6g W/cm^2' % (i_f, i_f_cm2))
    p_rad_n = i_n / C_LIGHT                    # 辐射压（工况 N）
    p_rad_f = i_f / C_LIGHT
    print('工况 N：辐射压 I/c = %.6g Pa = %.6g GPa' % (p_rad_n, p_rad_n / 1e9))
    print('工况 F：辐射压 I/c = %.6g Pa = %.6g GPa' % (p_rad_f, p_rad_f / 1e9))
    u_las_n = i_n / C_LIGHT                    # 能量密度（同量纲）
    print('工况 N：激光能量密度 I/c = %.6g J/m^3' % u_las_n)
    print('工况 F：激光能量密度 I/c = %.6g J/m^3' % (i_f / C_LIGHT))
    print('温度换算：1 eV = %.6g K' % EV_K)
    print('压强换算：1 Mbar = %.6g Pa = %.6g atm' % (MBAR, MBAR / 101325.0))

    # ================= AW.2 温稠密物质 =================
    print('\n[AW.2] 温稠密物质：简并与耦合同时在场')
    ni = NE_WDM / ZBAR
    e_f = (HBAR ** 2 / (2 * M_E)) * (3 * math.pi ** 2 * NE_WDM) ** (2.0 / 3.0)
    t_f_k = e_f / K_B
    t_f_ev = ev(e_f)
    a_i = (3.0 / (4 * math.pi * ni)) ** (1.0 / 3.0)
    gamma_ii = (ZBAR ** 2 * E_CHG ** 2) / (4 * math.pi * EPS0 * a_i * (TI_WDM_EV * E_CHG))
    theta = TE_WDM_EV / t_f_ev
    lam_de = math.sqrt(EPS0 * (TE_WDM_EV * E_CHG) / (NE_WDM * E_CHG ** 2))
    n_d = (4.0 / 3.0) * math.pi * NE_WDM * lam_de ** 3
    print('离子密度 n_i = n_e/Zbar = %.6g m^-3' % ni)
    print('费米能 E_F = %.6g J = %.6g eV' % (e_f, t_f_ev))
    print('费米温度 T_F = %.6g K' % t_f_k)
    print('离子球半径 a_i = %.6g m' % a_i)
    print('离子耦合参数 Gamma_ii = %.6g' % gamma_ii)
    print('简并参数 Theta = T_e/T_F = %.6g' % theta)
    print('德拜长度 lambda_De = %.6g m' % lam_de)
    print('等离子体参数 N_D = %.6g' % n_d)

    # ================= AW.3 状态方程 =================
    print('\n[AW.3] 状态方程：理想气体、简并压与声速')
    p_ideal = NE_WDM * (TE_WDM_EV * E_CHG)
    p_deg = 0.4 * NE_WDM * e_f
    p_ion = ni * (TI_WDM_EV * E_CHG)
    p_tot = p_ideal + p_deg + p_ion
    m_ion = M_ION_U * U_KG
    c_s = math.sqrt(ZBAR * (TE_WDM_EV * E_CHG) / m_ion)
    print('理想电子热压 n_e k T_e = %.6g Pa = %.6g Mbar' % (p_ideal, p_ideal / MBAR))
    print('简并压 (2/5) n_e E_F = %.6g Pa = %.6g Mbar' % (p_deg, p_deg / MBAR))
    print('离子热压 n_i k T_i = %.6g Pa = %.6g Mbar' % (p_ion, p_ion / MBAR))
    print('总压 = %.6g Pa = %.6g Mbar' % (p_tot, p_tot / MBAR))
    print('简并压/热压 = %.6g' % (p_deg / p_ideal))
    print('离子声速 c_s = %.6g m/s' % c_s)
    print('离子质量 m_i = %.6g kg' % m_ion)

    # ================= AW.4 强场电离与场效应 =================
    print('\n[AW.4] 强场电离与场效应（工况 F）')
    a0_f = math.sqrt(i_f_cm2 * LAM_UM ** 2 / 1.37e18)
    a0_n = math.sqrt(i_n_cm2 * LAM_UM ** 2 / 1.37e18)
    up_f_ev = (M_E * C_LIGHT ** 2 / E_CHG) * a0_f ** 2 / 4.0
    up_n_ev = (M_E * C_LIGHT ** 2 / E_CHG) * a0_n ** 2 / 4.0
    v_osc_f = a0_f / math.sqrt(1 + a0_f ** 2)
    gam_rel_f = math.sqrt(1 + a0_f ** 2)
    gamma_k_f = math.sqrt(E_ION_H_EV / (2 * up_f_ev))
    f_field_f = math.sqrt(2 * i_f / (EPS0 * C_LIGHT))
    f_bsi = (E_ION_H_J ** 2 * 4 * math.pi * EPS0) / (4 * E_CHG ** 3)
    i_bsi = 0.5 * EPS0 * C_LIGHT * f_bsi ** 2
    lam_m = LAM_UM * 1e-6
    omega_l = 2 * math.pi * C_LIGHT / lam_m
    n_c = EPS0 * M_E * omega_l ** 2 / E_CHG ** 2
    ratio_nc = NE_WDM / n_c
    a0_trans = ratio_nc                       # 透明条件 gamma = sqrt(1+a_0^2) > n_e/n_c
    i_trans = 1.37e18 * a0_trans ** 2 / LAM_UM ** 2
    print('工况 F：归一化矢量势 a_0 = %.6g' % a0_f)
    print('工况 N：归一化矢量势 a_0 = %.6g' % a0_n)
    print('工况 N：抖动速度 v_osc/c = %.6g' % (a0_n / math.sqrt(1 + a0_n ** 2)))
    print('工况 N：相对论因子 gamma = %.6g' % math.sqrt(1 + a0_n ** 2))
    print('工况 N：激光电场 F = %.6g V/m' % math.sqrt(2 * i_n / (EPS0 * C_LIGHT)))
    print('工况 F：有质动力势 U_p = %.6g eV = %.6g keV' % (up_f_ev, up_f_ev / 1e3))
    print('工况 N：有质动力势 U_p = %.6g eV' % up_n_ev)
    print('工况 F：电子抖动速度 v_osc/c = %.6g' % v_osc_f)
    print('工况 F：电子抖动速度百分数 = %.6g%%' % (v_osc_f * 100))
    print('工况 F：相对论因子 gamma = %.6g' % gam_rel_f)
    print('工况 F：Keldysh 参数 gamma_K = %.6g' % gamma_k_f)
    print('工况 F：激光电场 F = %.6g V/m' % f_field_f)
    print('势垒抑制场 F_BSI = %.6g V/m' % f_bsi)
    print('势垒抑制强度 I_BSI = %.6g W/m^2 = %.6g W/cm^2' % (i_bsi, i_bsi / 1e4))
    print('临界密度 n_c = %.6g m^-3' % n_c)
    print('固体密度/临界密度 = %.6g' % ratio_nc)
    print('相对论透明所需 a_0 = %.6g' % a0_trans)
    print('相对论透明所需强度 = %.6g W/cm^2 = %.6g W/m^2' % (i_trans, i_trans * 1e4))

    # ================= AW.5 冲击与烧蚀 =================
    print('\n[AW.5] 冲击与烧蚀：Rankine-Hugoniot 与 Mbar 量级')
    r_comp = (GAMMA + 1) / (GAMMA - 1)
    u_p = D_SHOCK / r_comp
    p_shock = RHO0 * D_SHOCK * u_p
    kt_shock = 0.1875 * m_ion * u_p ** 2
    t_shock_ev = ev(kt_shock)
    p_abl_mbar = 8.6 * (i_n_cm2 / 1e14) ** (2.0 / 3.0) * LAM_UM ** (-2.0 / 3.0)
    p_abl = p_abl_mbar * MBAR
    print('强冲击压缩比 r = (gamma+1)/(gamma-1) = %.6g' % r_comp)
    print('粒子速度 u_p = D/r = %.6g m/s' % u_p)
    print('冲击压强 P = rho_0 D u_p = %.6g Pa = %.6g Mbar' % (p_shock, p_shock / MBAR))
    print('冲击离子温度 k T_i = %.6g J = %.6g eV' % (kt_shock, t_shock_ev))
    print('烧蚀压强（标度律） = %.6g Mbar = %.6g Pa = %.6g GPa' % (p_abl_mbar, p_abl, p_abl / 1e9))
    print('烧蚀压强/辐射压 = %.6g' % (p_abl / p_rad_n))

    # ================= AW.6 辐射温度与黑腔 =================
    print('\n[AW.6] 辐射温度与黑腔：四次方根')
    f_x = ETA_X * i_n
    t_rad_k = (f_x / SIGMA_SB) ** 0.25
    t_rad_ev = t_rad_k / EV_K
    u_rad = 4 * SIGMA_SB * t_rad_k ** 4 / C_LIGHT
    p_rad_bb = u_rad / 3.0
    t_wall_k = ((1 - ALPHA_WALL) * f_x / SIGMA_SB) ** 0.25
    t_wall_ev = t_wall_k / EV_K
    e_x = ETA_X * E_L_N
    print('X 射线通量 eta_X I = %.6g W/m^2' % f_x)
    print('辐射温度 T_rad = %.6g K = %.6g eV' % (t_rad_k, t_rad_ev))
    print('黑体能量密度 u_rad = %.6g J/m^3' % u_rad)
    print('黑体辐射压 u_rad/3 = %.6g Pa = %.6g Mbar' % (p_rad_bb, p_rad_bb / MBAR))
    print('壁吸收后有效温度 = %.6g K = %.6g eV' % (t_wall_k, t_wall_ev))
    print('X 射线能量 eta_X E_L = %.6g J' % e_x)
    print('四次方根灵敏度：通量 +10%% 时温度涨 %.6g%%' % ((1.1 ** 0.25 - 1) * 100))

    # ================= AW.7 实验室天体物理 =================
    print('\n[AW.7] 实验室天体物理：无量纲相似')
    rm = MU0 * SIGMA_COND * V_ASTRO * L_ASTRO
    tau_eta = MU0 * SIGMA_COND * L_ASTRO ** 2
    v_a = B_ASTRO / math.sqrt(MU0 * RHO0)
    beta_pl = 2 * MU0 * P_ASTRO / B_ASTRO ** 2
    rm_want = 1.0e4
    print('磁雷诺数 Rm = mu0 sigma v L = %.6g' % rm)
    print('磁扩散时间 tau_eta = mu0 sigma L^2 = %.6g s = %.6g us' % (tau_eta, tau_eta * 1e6))
    print('阿尔芬速度 v_A = B/sqrt(mu0 rho) = %.6g m/s' % v_a)
    print('等离子体 beta = 2 mu0 P/B^2 = %.6g' % beta_pl)
    print('要达到天体量级 Rm = %.6g 需要提高 %.6g 倍' % (rm_want, rm_want / rm))

    # ================= AW.8 练习解答 =================
    print('\n[AW.8] 练习解答')
    i_ex1 = 1.0e15              # W/cm^2
    a0_ex1 = math.sqrt(i_ex1 * LAM_UM ** 2 / 1.37e18)
    up_ex1 = (M_E * C_LIGHT ** 2 / E_CHG) * a0_ex1 ** 2 / 4.0
    print('练习1：a_0 = %.6g，U_p = %.6g eV = %.6g keV，v_osc/c = %.6g'
          % (a0_ex1, up_ex1, up_ex1 / 1e3, a0_ex1 / math.sqrt(1 + a0_ex1 ** 2)))
    ne_ex2 = 1.0e30
    ef_ex2 = (HBAR ** 2 / (2 * M_E)) * (3 * math.pi ** 2 * ne_ex2) ** (2.0 / 3.0)
    ef_ex2_ev = ev(ef_ex2)
    theta_ex2 = 5.0 / ef_ex2_ev
    p_deg_ex2 = 0.4 * ne_ex2 * ef_ex2
    print('练习2：E_F = %.6g eV，Theta(5 eV) = %.6g，简并压 = %.6g Pa = %.6g Mbar'
          % (ef_ex2_ev, theta_ex2, p_deg_ex2, p_deg_ex2 / MBAR))
    d_ex3 = 100e3
    up_ex3 = d_ex3 / r_comp
    p_ex3 = RHO0 * d_ex3 * up_ex3
    t_ex3 = ev(0.1875 * m_ion * up_ex3 ** 2)
    print('练习3：u_p = %.6g m/s，P = %.6g Mbar，k T_i = %.6g eV' % (up_ex3, p_ex3 / MBAR, t_ex3))
    f_ex4 = 0.5 * 1e15 * 1e4
    t_ex4 = (f_ex4 / SIGMA_SB) ** 0.25 / EV_K
    print('练习4：X 射线通量 = %.6g W/m^2' % f_ex4)
    print('练习4：T_rad = %.6g eV' % t_ex4)
    lam_ex5 = 0.8
    a0_ex5 = math.sqrt(1e19 * lam_ex5 ** 2 / 1.37e18)
    up_ex5 = (M_E * C_LIGHT ** 2 / E_CHG) * a0_ex5 ** 2 / 4.0
    gk_ex5 = math.sqrt(E_ION_H_EV / (2 * up_ex5))
    print('练习5：a_0 = %.6g，U_p = %.6g keV = %.6g eV，gamma_K = %.6g'
          % (a0_ex5, up_ex5 / 1e3, up_ex5, gk_ex5))
    print('练习6：Rm 需提高 %.6g 倍' % (rm_want / rm))


if __name__ == '__main__':
    main()
