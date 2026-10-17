# -*- coding: utf-8 -*-
r"""附录 AM 磁体与超导：数值唯一来源（只 print，不写文件）。

用法：python tmp\appAM_calc.py | Out-File -Encoding utf8 tmp\appAM_calc.txt

工作点与附录 AG/AL 同一台装置：R0 = 6.2 m、a = 2.0 m、B_TF = 5.3 T、I_p = 11.0844 MA。
本脚本只做量级估算，所有「声明量」都打印在 [DECLARED] 段，供正文与门禁引用。
"""
from __future__ import print_function

import math

# ---------------------------------------------------------------- 常数
mu0 = 4.0e-7 * math.pi          # 1.2566370614359173e-06
e = 1.602176634e-19
kB = 1.380649e-23
u = 1.66053906660e-27

# ---------------------------------------------------------------- 声明量
DECLARED = [
    ('R0', 6.2, 'm', '大半径（与附录 AG/AL 同点）'),
    ('a', 2.0, 'm', '小半径'),
    ('B_TF', 5.3, 'T', '环向场（在 R0 处）'),
    ('I_p', 11.0844e6, 'A', '等离子体电流（附录 AG 的 q95=3 工作点）'),
    ('I_op', 68.0e3, 'A', '导体工作电流（ITER 量级）'),
    ('n_TF', 18, '-', '环向场线圈数（ITER 量级）'),
    ('J_op', 40.0e6, 'A/m^2', '电缆平均电流密度 40 A/mm^2'),
    ('f_noncu', 0.40, '-', '电缆中非铜（超导）面积份额'),
    ('f_pack', 0.50, '-', '绕组包中电缆的填充因子'),
    ('R_in', 4.2, 'm', '内腿半径'),
    ('h_inner', 8.0, 'm', '内腿高度'),
    ('a_B', 3.5, 'm', '场区等效半径（储能估算用）'),
    ('Jc0', 6000.0e6, 'A/m^2', 'Nb3Sn 参数化 Jc0 = 6000 A/mm^2'),
    ('B_c2', 26.0, 'T', 'Nb3Sn 上临界场参数化值'),
    ('beta', 1.5, '-', 'Jc(B) 的幂指数'),
    ('T_c0', 18.0, 'K', 'Nb3Sn 零场临界温度'),
    ('T_op', 4.5, 'K', '运行温度（超流氦）'),
    ('T_bath', 4.5, 'K', '冷媒温度'),
    ('T_warm', 300.0, 'K', '室温端温度'),
    ('eta_carnot', 0.25, '-', '制冷机相对卡诺效率'),
    ('tau_dump', 10.0, 's', '泄能时间常数'),
    ('tau_detect', 1.0, 's', '失超检测延迟'),
    ('tau_fail', 1000.0, 's', '保护失效时的等效时间常数'),
    ('C_cu', 1.5e6, 'J/(m^3 K)', '铜等效体积热容（4--150 K 平均）'),
    ('rho_cu', 5.0e-9, 'Ohm m', '铜等效电阻率（4--150 K 平均，RRR~100）'),
    ('rho_cu_warm', 1.7e-8, 'Ohm m', '室温铜电阻率'),
    ('J_cu_water', 5.0e6, 'A/m^2', '水冷铜允许电流密度 5 A/mm^2'),
    ('Q_cold', 5.0e3, 'W', '4.5 K 冷负荷（含交流损耗与支撑漏热）'),
    ('n_tau', 0.05, 's', '耦合损耗时间常数'),
    ('dBdt', 0.5, 'T/s', '磁场变化率（脉冲段）'),
    ('V_wp_total', 8.2, 'm^3', '全部 TF 绕组包体积（由本模型算得，此处只作校验）'),
    ('sigma_yield', 900.0e6, 'Pa', '高强度钢屈服强度量级'),
]


def show_declared():
    print('[DECLARED] 声明量（写进正文的输入，不是算出来的）')
    for name, val, unit, why in DECLARED:
        print('  %-14s = %-14.6g %-10s  # %s' % (name, val, unit, why))
    print('')


def main():
    show_declared()

    # ------------------------------------------------------------ AM.1
    # 安培环路：NI = 2 pi R B / mu0
    NI = 2.0 * math.pi * 6.2 * 5.3 / mu0
    N_per_coil = NI / 68.0e3 / 18.0
    L_turn = 2.0 * math.pi * 8.0            # 平均匝长（取 R_avg = 8 m）
    L_per_coil = N_per_coil * L_turn
    L_all = L_per_coil * 18.0
    A_cab = 68.0e3 / 40.0e6
    A_wp_per_coil = A_cab * N_per_coil / 0.50
    A_wp_all = A_wp_per_coil * 18.0
    V_wp = A_wp_all * L_turn                # 每匝绕组包体积 × 匝数
    print('[AM.1] 为什么必须超导：同一套安培匝数的两条路')
    print('  NI                     = %.6g A-turns' % NI)
    print('  N_per_coil             = %.6g 匝' % N_per_coil)
    print('  L_turn                 = %.6g m' % L_turn)
    print('  L_per_coil             = %.6g m' % L_per_coil)
    print('  L_all                  = %.6g m = %.6g km' % (L_all, L_all / 1e3))
    print('  A_cab                  = %.6g m^2 = %.6g mm^2' % (A_cab, A_cab * 1e6))
    print('  A_wp_per_coil          = %.6g m^2' % A_wp_per_coil)
    print('  V_wp                   = %.6g m^3' % V_wp)
    # 铜方案：同样电流、水冷铜 J = 5 A/mm^2
    A_cu = 68.0e3 / 5.0e6
    R_cu = 1.7e-8 * L_all / A_cu
    P_cu = (68.0e3) ** 2 * R_cu
    print('  A_cu (5 A/mm^2)        = %.6g m^2 = %.6g mm^2' % (A_cu, A_cu * 1e6))
    print('  R_cu                   = %.6g Ohm' % R_cu)
    print('  P_cu                   = %.6g W = %.6g MW = %.6g GW'
          % (P_cu, P_cu / 1e6, P_cu / 1e9))
    print('  P_cu / P_fus(556.5MW)  = %.6g 倍' % (P_cu / 556.5e6))
    # 超导方案：只有制冷电耗
    T_ratio = 300.0 / 4.5
    W_carnot = 5.0e3 * (T_ratio - 1.0)
    P_cryo = W_carnot / 0.25
    print('  T_ratio (300/4.5)      = %.6g' % T_ratio)
    print('  W_carnot (5 kW)        = %.6g W = %.6g kW' % (W_carnot, W_carnot / 1e3))
    print('  P_cryo (25%% Carnot)    = %.6g W = %.6g kW' % (P_cryo, P_cryo / 1e3))
    print('  P_cu / P_cryo          = %.6g 倍' % (P_cu / P_cryo))
    print('')

    # ------------------------------------------------------------ AM.2
    # 临界面：Jc(B)、Ic、工作点裕度；Tc(B)、温度裕度
    B = 5.3
    Jc_53 = 6000.0e6 * (1.0 - B / 26.0) ** 1.5
    Jc_12 = 6000.0e6 * (1.0 - 12.0 / 26.0) ** 1.5
    Jc_20 = 6000.0e6 * (1.0 - 20.0 / 26.0) ** 1.5
    Ic_cab = Jc_53 * A_cab * 0.40
    Ic_12 = Jc_12 * A_cab * 0.40
    Ic_20 = Jc_20 * A_cab * 0.40
    Tc_53 = 18.0 * math.sqrt(1.0 - 5.3 / 26.0)
    Tc_12 = 18.0 * math.sqrt(1.0 - 12.0 / 26.0)
    J_op_nc = 40.0e6 / 0.40
    T_cs = Tc_53 - (Tc_53 - 4.5) * J_op_nc / Jc_53
    print('[AM.2] 临界面：Jc(B, T) 与工作点裕度')
    print('  Jc(5.3 T)              = %.6g A/m^2 = %.6g A/mm^2' % (Jc_53, Jc_53 / 1e6))
    print('  Jc(12 T)               = %.6g A/m^2 = %.6g A/mm^2' % (Jc_12, Jc_12 / 1e6))
    print('  Jc(20 T)               = %.6g A/m^2 = %.6g A/mm^2' % (Jc_20, Jc_20 / 1e6))
    print('  Ic_cab(5.3 T)          = %.6g A = %.6g MA' % (Ic_cab, Ic_cab / 1e6))
    print('  Ic_cab(12 T)           = %.6g A = %.6g MA' % (Ic_12, Ic_12 / 1e6))
    print('  Ic_cab(20 T)           = %.6g A = %.6g MA' % (Ic_20, Ic_20 / 1e6))
    print('  I_op / Ic(5.3 T)       = %.6g' % (68.0e3 / Ic_cab))
    print('  I_op / Ic(12 T)        = %.6g' % (68.0e3 / Ic_12))
    print('  I_op / Ic(20 T)        = %.6g' % (68.0e3 / Ic_20))
    print('  Tc(5.3 T)              = %.6g K' % Tc_53)
    print('  Tc(12 T)               = %.6g K' % Tc_12)
    print('  J_op(非铜)             = %.6g A/m^2 = %.6g A/mm^2' % (J_op_nc, J_op_nc / 1e6))
    print('  T_cs (分流温度)        = %.6g K' % T_cs)
    print('  温度裕度 T_cs - T_op   = %.6g K' % (T_cs - 4.5))
    print('  n 值经验式裕度：I/Ic    = %.6g -> 电压 %.6g 倍基准'
          % (68.0e3 / Ic_cab, (68.0e3 / Ic_cab) ** 20))
    print('')

    # ------------------------------------------------------------ AM.2b
    # CICC 与稳定性判据（Stekly）
    A_cu = A_cab * 0.60
    rho_cu_cold = 1.0e-9
    R_prime = rho_cu_cold / A_cu
    h_kap = 1000.0
    P_wet = 0.20
    dT_stekly = 1.0
    Q_stekly = h_kap * P_wet * dT_stekly
    I_stekly = math.sqrt(Q_stekly / R_prime)
    alpha_stekly = (68.0e3) ** 2 * R_prime / Q_stekly
    print('[AM.2b] CICC 与稳定性：Stekly 判据够不够用')
    print('  A_cu（电缆中铜）       = %.6g m^2 = %.6g mm^2' % (A_cu, A_cu * 1e6))
    print('  R\' = rho/A_cu         = %.6g Ohm/m' % R_prime)
    print('  Q\' = h P dT           = %.6g W/m' % Q_stekly)
    print('  I_stekly               = %.6g A = %.6g kA' % (I_stekly, I_stekly / 1e3))
    print('  alpha (I=68 kA)        = %.6g' % alpha_stekly)
    print('  alpha=1 允许电流       = %.6g kA' % (68.0e3 / math.sqrt(alpha_stekly) / 1e3))
    print('')


    # 储能、电感、泄能、热点
    p_B = B ** 2 / (2.0 * mu0)
    V_B = 2.0 * math.pi * 6.2 * math.pi * 3.5 ** 2
    W = p_B * V_B
    Lamp = 2.0 * W / (68.0e3) ** 2
    R_dump = Lamp / 10.0
    V_peak = 68.0e3 * R_dump
    P_peak = (68.0e3) ** 2 * R_dump
    # 热点：铜里的积分 J^2 dt
    J_cu = 40.0e6 * 0.60 / 1.0            # 电缆中铜的电流密度（Cu:nonCu=1.5 -> 铜占 60%）
    J_cu = 68.0e3 / (A_cab * 0.60)
    miits_ok = J_cu ** 2 * (1.0 + 10.0 / 2.0)
    miits_fail = J_cu ** 2 * (1.0 + 1000.0 / 2.0)
    dT_ok = miits_ok * 5.0e-9 / 1.5e6
    dT_fail = miits_fail * 5.0e-9 / 1.5e6
    miits_150 = 1.5e6 * 150.0 / 5.0e-9
    print('[AM.3] 储能与失超：能量能去哪里')
    print('  p_B = B^2/(2 mu0)      = %.6g Pa = %.6g MPa = %.6g bar'
          % (p_B, p_B / 1e6, p_B / 1e5))
    print('  V_B (2 pi R pi a_B^2)  = %.6g m^3' % V_B)
    print('  W = p_B V_B            = %.6g J = %.6g MJ = %.6g GJ'
          % (W, W / 1e6, W / 1e9))
    print('  L (2W/I^2)             = %.6g H' % Lamp)
    print('  R_dump (tau=10 s)      = %.6g Ohm' % R_dump)
    print('  V_peak                 = %.6g V = %.6g kV' % (V_peak, V_peak / 1e3))
    print('  P_peak                 = %.6g W = %.6g GW' % (P_peak, P_peak / 1e9))
    print('  J_cu (铜中)            = %.6g A/m^2 = %.6g A/mm^2' % (J_cu, J_cu / 1e6))
    print('  MIITs 正常保护         = %.6g A^2 s/m^4' % miits_ok)
    print('  dT 正常保护            = %.6g K' % dT_ok)
    print('  热点温度（4.5 K 冷端） = %.6g K' % (4.5 + dT_ok))
    print('  MIITs 保护失效         = %.6g A^2 s/m^4' % miits_fail)
    print('  dT 保护失效            = %.6g K' % dT_fail)
    print('  MIITs 允许(150 K)      = %.6g A^2 s/m^4' % miits_150)
    print('  失效/允许              = %.6g 倍' % (miits_fail / miits_150))
    print('')

    # ------------------------------------------------------------ AM.4
    # 力与应力
    f_line = 68.0e3 * 5.3
    F_coil = f_line * 8.0 * 134.0
    F_all = F_coil * 18.0
    sigma = F_coil / A_wp_per_coil
    print('[AM.4] 洛伦兹力与应力：为什么线圈要装进钢盒')
    print('  f = I B                = %.6g N/m = %.6g kN/m' % (f_line, f_line / 1e3))
    print('  F 一个线圈（内腿）     = %.6g N = %.6g MN' % (F_coil, F_coil / 1e6))
    print('  F 全部 18 个线圈       = %.6g N = %.6g GN' % (F_all, F_all / 1e9))
    print('  sigma = F/A_wp         = %.6g Pa = %.6g MPa' % (sigma, sigma / 1e6))
    print('  sigma / sigma_yield    = %.6g' % (sigma / 900.0e6))
    print('  (F/A_wp) 若绕组包减半  = %.6g MPa' % (F_coil / (A_wp_per_coil / 2.0) / 1e6))
    print('')

    # ------------------------------------------------------------ AM.5
    # 制冷、交流损耗
    Q_couple = (0.05 / mu0) * 0.5 ** 2          # W/m^3，斜坡期间的耦合损耗功率密度
    V_wp_all = A_wp_per_coil * 18.0 * L_turn
    P_couple_inst = Q_couple * V_wp_all
    E_couple = P_couple_inst * 10.0             # 10 s 斜坡
    P_couple_avg = E_couple / 1000.0            # 1000 s 脉冲周期
    print('[AM.5] 制冷与交流损耗：每一瓦都要乘一个系数')
    print('  Q_couple 密度          = %.6g W/m^3（斜坡期间）' % Q_couple)
    print('  V_wp 全部              = %.6g m^3' % V_wp_all)
    print('  P_couple 瞬时          = %.6g W = %.6g MW' % (P_couple_inst, P_couple_inst / 1e6))
    print('  E_couple 一次脉冲      = %.6g J = %.6g MJ' % (E_couple, E_couple / 1e6))
    print('  平均功率(1000 s 周期)  = %.6g W = %.6g kW' % (P_couple_avg, P_couple_avg / 1e3))
    print('  冷负荷 5 kW 的卡诺功   = %.6g W = %.6g kW' % (W_carnot, W_carnot / 1e3))
    print('  制冷电耗(25%%)         = %.6g kW = %.6g MW' % (P_cryo / 1e3, P_cryo / 1e6))
    print('  交流损耗的制冷电耗     = %.6g W = %.6g MW' % (P_couple_avg * (300.0 / 4.5 - 1.0) / 0.25, P_couple_avg * (300.0 / 4.5 - 1.0) / 0.25 / 1e6))
    P_cryo_tot = P_cryo + P_couple_avg * (300.0 / 4.5 - 1.0) / 0.25
    print('  制冷总电耗             = %.6g W = %.6g MW' % (P_cryo_tot, P_cryo_tot / 1e6))
    print('  制冷总电耗 / D-T 输出  = %.6g' % (P_cryo_tot / 556.5e6))
    print('  制冷电耗 / D-T 输出    = %.6g' % (P_cryo / 556.5e6))
    print('')

    # ------------------------------------------------------------ AM.6
    # 练习解答
    print('[AM.6] 估算练习解答')
    # 1) B 翻倍：储能、力、成本
    W2 = (2 * 5.3) ** 2 / (2 * mu0) * V_B
    F2 = 68.0e3 * 10.6 * 8.0 * 134.0
    print('  1) B=10.6 T: W = %.6g GJ (%.6g 倍), F = %.6g MN (%.6g 倍)'
          % (W2 / 1e9, W2 / W, F2 / 1e6, F2 / F_coil))
    # 2) 铜方案的泵功耗：水冷铜压降 5 bar、流速
    P_cu2 = (68.0e3) ** 2 * (1.7e-8 * L_all / (68.0e3 / 5.0e6))
    print('  2) 铜方案 P = %.6g GW；若 10 A/mm^2 则 P = %.6g GW'
          % (P_cu2 / 1e9, ((68.0e3) ** 2 * (1.7e-8 * L_all / (68.0e3 / 1.0e7))) / 1e9))
    # 3) 泄能时间常数 5 s 与 20 s 的热点
    dT5 = J_cu ** 2 * (1.0 + 5.0 / 2.0) * 5.0e-9 / 1.5e6
    dT20 = J_cu ** 2 * (1.0 + 20.0 / 2.0) * 5.0e-9 / 1.5e6
    print('  3) tau=5 s: dT = %.6g K；tau=20 s: dT = %.6g K' % (dT5, dT20))
    # 4) 冷负荷 20 kW 的电耗
    print('  4) 20 kW 冷负荷：%.6g W = %.6g MW 电耗'
          % (20.0e3 * (300.0 / 4.5 - 1.0) / 0.25, 20.0e3 * (300.0 / 4.5 - 1.0) / 0.25 / 1e6))
    # 5) Ic 降 10% 对 n 值电压
    print('  5) Jc 降 10%% 后 I/Ic = %.6g，电压比 = %.6g'
          % (68.0e3 / (0.9 * Ic_cab), (68.0e3 / (0.9 * Ic_cab)) ** 20 / (68.0e3 / Ic_cab) ** 20))
    # 6) 若允用 40% 应力余量，绕组包要多大
    A_need = F_coil / (0.6 * 900.0e6)
    print('  6) 要把应力压到 540 MPa，A_wp 需 %.6g m^2（现 %.6g m^2）'
          % (A_need, A_wp_per_coil))
    print('')


if __name__ == '__main__':
    main()
