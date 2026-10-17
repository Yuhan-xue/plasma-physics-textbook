# -*- coding: utf-8 -*-
"""附录 AD（边界台基与 L-H 模）的数值唯一来源。

工作点 D1（与附录 M / S / U / V / W / X / Y / Z / AB / AC 同一套输入）：
    n_e = 1e20 m^-3, T_e = T_i = 9.194 keV, B = 5.3 T, R = 6.2 m, a = 2.0 m,
    V = 832.2 m^3, W_th = 371.4 MJ, P_fus = 556.5 MW, P_alpha = 111.3 MW,
    I_p = 15 MA, q_95 = 3.5, q_0 = 1.0, m_i = 2.5 m_p, Z_eff = 2,
    第一壁/分界面面积 S = 489.532 m^2（与附录 X 同一数值）

本脚本只打印，不写任何文件；正文里出现的每个数字都必须能在这里找到。
运行：python tmp\appAD_calc.py
"""
import math

# ---------------------------------------------------------------- 常数
e = 1.602176634e-19          # C
mp = 1.67262192369e-27       # kg
mu0 = 4.0e-7 * math.pi       # H/m
keV = 1.0e3 * e              # J
keVJ = 1.602176634e-16       # J per keV（便于手算核对）

print('=' * 78)
print('附录 AD 数值来源：边界台基与 L-H 模（工作点 D1）')
print('=' * 78)

# ---------------------------------------------------------------- [DECLARED]
n_e = 1.0e20            # m^-3
T_e = 9.194             # keV
T_i = 9.194             # keV
B = 5.3                 # T
R = 6.2                 # m
a = 2.0                 # m
V_pl = 832.2            # m^3
W_th = 371.4e6          # J
P_fus = 556.5e6         # W
P_alpha = 111.3e6       # W
P_aux = 20.0e6          # W
P_rad_core = 27.8322e6  # W（附录 X 的总辐射功率）
I_p = 15.0e6            # A
q_95 = 3.5
q_0 = 1.0
m_i = 2.5 * mp          # kg
Z_eff = 2.0
S_wall = 489.532        # m^2（= 4 pi^2 R a，与附录 X 一致）
# 边界/台基（声明值：台基顶取芯部密度 60%、温度 3 keV）
n_ped = 0.60 * n_e      # m^-3
T_ped = 3.0             # keV
# Martin 2008 L-H 阈值定标（声明系数）
C_th, ex_n, ex_B, ex_S = 0.0488, 0.717, 0.803, 0.941
# EPED 台基宽度定标（声明系数）：Delta = 0.076 sqrt(beta_p,ped) R
C_w, ex_w = 0.076, 0.5
# 理想气球模判据（声明）：alpha_crit = 2.5（无自举）/ 4.0（含自举）
alpha_crit_no_bs, alpha_crit_bs = 2.5, 4.0
# 自举电流的序系数（声明：Wesson 标准形式，序 1 系数取 1.5）
C_bs = 1.5
# ELM：每次损失 7% W_th，上升时间 0.3 ms，靶板极向长度 1 m
f_ELM_W, t_ELM, L_div = 0.07, 0.3e-3, 1.0
# 通量展宽与外侧偏滤器份额（声明）
f_exp, f_out = 5.0, 0.70
# 钨靶板能量密度阈值（声明：0.3 ms 量级脉冲）
E_th_W = 1.0e6          # J/m^2
# 稳态热流材料上限（声明）
q_limit = 10.0e6        # W/m^2
# 能量约束时间（声明，H 模）
tau_E = 3.7             # s
# H -> L 回滞比（声明）
hyst = 0.65
# 弹丸（声明：1e21 个原子，全部沉积在台基壳层）
N_pellet = 1.0e21
m_atom = 2.5 * mp

print('[DECLARED] 输入一览')
print('  n_e=%.3g m^-3  T_e=T_i=%.4g keV  B=%.3g T  R=%.3g m  a=%.3g m  V=%.5g m^3'
      % (n_e, T_e, B, R, a, V_pl))
print('  W_th=%.5g MJ  P_fus=%.5g MW  P_alpha=%.5g MW  P_aux=%.4g MW  P_rad_core=%.6g MW'
      % (W_th / 1e6, P_fus / 1e6, P_alpha / 1e6, P_aux / 1e6, P_rad_core / 1e6))
print('  I_p=%.3g MA  q_95=%.3g  q_0=%.3g  m_i=%.3g m_p  Z_eff=%.3g  S_wall=%.6g m^2'
      % (I_p / 1e6, q_95, q_0, m_i / mp, Z_eff, S_wall))
print('  n_ped=%.4g m^-3（0.6 n_e）  T_ped=%.3g keV' % (n_ped, T_ped))
print('  Martin: P_th=%.4g n20^%.3g B^%.3g S^%.3g MW' % (C_th, ex_n, ex_B, ex_S))
print('  EPED: Delta_ped=%.3g sqrt(beta_p,ped) R m；alpha_crit=%.3g/%.3g（无/含自举）'
      % (C_w, alpha_crit_no_bs, alpha_crit_bs))
print('  ELM: dW/W_th=%.3g  t_ELM=%.3g ms  L_div=%.3g m  E_th(W)=%.3g MJ/m^2'
      % (f_ELM_W, t_ELM * 1e3, L_div, E_th_W / 1e6))
print('  f_exp=%.3g  f_out=%.3g  q_limit=%.3g MW/m^2  tau_E=%.3g s  N_pellet=%.3g'
      % (f_exp, f_out, q_limit / 1e6, tau_E, N_pellet))

# ---------------------------------------------------------------- [AD.1]
print('')
print('[AD.1] 从 L 模到 H 模：阈值、密度极限与回滞')
P_th = C_th * (n_e / 1e20) ** ex_n * B ** ex_B * S_wall ** ex_S * 1e6   # W
P_SOL = P_alpha + P_aux - P_rad_core                                    # W
n_G = (I_p / 1e6) / (math.pi * a * a) * 1e20          # Greenwald: 10^20 m^-3 per MA/(pi a^2)
f_G = n_e / n_G
P_back = hyst * P_th
print('  L-H 阈值 P_th = %.6g MW' % (P_th / 1e6))
print('  P_alpha / P_th = %.6g' % (P_alpha / P_th))
print('  P_SOL = P_alpha + P_aux - P_rad_core = %.6g MW' % (P_SOL / 1e6))
print('  P_SOL / P_th = %.6g' % (P_SOL / P_th))
print('  余量 P_SOL - P_th = %.6g MW' % ((P_SOL - P_th) / 1e6))
print('  格林沃尔德密度 n_G = I_p/(pi a^2) = %.6g m^-3 = %.6g 10^19 m^-3'
      % (n_G, n_G / 1e19))
print('  密度分数 f_G = n_e / n_G = %.6g' % f_G)
print('  回滞：H->L 反向阈值 = %.3g P_th = %.6g MW' % (hyst, P_back / 1e6))
print('  阈值对密度的敏感度：0.5 n_e -> P_th = %.6g MW，2 n_e -> P_th = %.6g MW'
      % (P_th * 0.5 ** ex_n / 1e6, P_th * 2.0 ** ex_n / 1e6))
print('  阈值对磁场的敏感度：B=4 T -> P_th = %.6g MW，B=6 T -> P_th = %.6g MW'
      % (P_th * (4.0 / B) ** ex_B / 1e6, P_th * (6.0 / B) ** ex_B / 1e6))
c_s = math.sqrt(T_e * keV / m_i)
gamma_turb = c_s / a
print('  离子声速 c_s = sqrt(T_e/m_i) = %.6g m/s = %.6g km/s' % (c_s, c_s / 1e3))
print('  湍流速率（序估计）gamma ~ c_s/a = %.6g s^-1' % gamma_turb)
v_sh = gamma_turb * 0.119603          # 用稍后算出的 Delta_ped 量级，先按 12 cm 估
E_r_need = v_sh * 1.5                 # B_theta = 1.5 T
print('  抑制湍流所需极向速度 v_theta >= gamma Delta_ped(12 cm) = %.6g m/s' % v_sh)
print('  对应径向电场 E_r = v_theta B_theta = %.6g V/m = %.6g kV/m'
      % (E_r_need, E_r_need / 1e3))
E_r_typ = 30.0e3
print('  若实测 E_r = %.3g kV/m，则剪切率 omega = E_r/(B_theta Delta_ped) = %.6g s^-1 = %.6g gamma'
      % (E_r_typ / 1e3, E_r_typ / (1.5 * 0.119603), E_r_typ / (1.5 * 0.119603) / gamma_turb))

# ---------------------------------------------------------------- [AD.2]
print('')
print('[AD.2] 台基：宽度、梯度与自举电流')
B_theta = mu0 * I_p / (2.0 * math.pi * a)
rho_i = math.sqrt(2.0 * m_i * T_i * keV) / (e * B)
rho_theta = rho_i * B / B_theta
k_theta = 1.0 / rho_i
p_ped = 2.0 * n_ped * T_ped * keV
p_avg = 2.0 * W_th / (3.0 * V_pl)
beta_p_ped = 2.0 * mu0 * p_ped / B_theta ** 2
Delta_ped = C_w * beta_p_ped ** ex_w * R
grad_p = p_ped / Delta_ped
alpha_MHD = 2.0 * mu0 * R * q_95 ** 2 * grad_p / B ** 2
beta_ped = 2.0 * mu0 * p_ped / B ** 2
j_bs = C_bs * math.sqrt(a / R) / B_theta * grad_p
j_avg = I_p / (math.pi * a * a)
I_bs = j_bs * Delta_ped * 2.0 * math.pi * R
eps = a / R
r_ped = a - Delta_ped
print('  极向磁场（边界）B_theta = mu0 I_p/(2 pi a) = %.6g T' % B_theta)
print('  离子拉莫尔半径 rho_i(T_i,B) = %.6g m = %.6g mm' % (rho_i, rho_i * 1e3))
print('  极向拉莫尔半径 rho_theta = rho_i B/B_theta = %.6g m = %.6g cm'
      % (rho_theta, rho_theta * 1e2))
print('  湍流波数 k_theta = 1/rho_i = %.6g m^-1（k_theta rho_i = 1）' % k_theta)
print('  台基顶压强 p_ped = 2 n_ped T_ped = %.6g Pa = %.6g kPa' % (p_ped, p_ped / 1e3))
print('  体积平均压强 <p> = 2 W_th/(3 V) = %.6g Pa = %.6g kPa' % (p_avg, p_avg / 1e3))
print('  p_ped / <p> = %.6g' % (p_ped / p_avg))
print('  台基极向 beta：beta_p,ped = 2 mu0 p_ped/B_theta^2 = %.6g' % beta_p_ped)
print('  台基总 beta：beta_ped = 2 mu0 p_ped/B^2 = %.6g = %.4g %%'
      % (beta_ped, beta_ped * 100))
print('  EPED 宽度 Delta_ped = %.3g sqrt(beta_p,ped) R = %.6g m = %.6g cm'
      % (C_w, Delta_ped, Delta_ped * 1e2))
print('  Delta_ped / rho_theta = %.6g' % (Delta_ped / rho_theta))
print('  台基梯度 |dp/dr| = p_ped/Delta_ped = %.6g Pa/m' % grad_p)
print('  气球参数 alpha_MHD = 2 mu0 R q95^2 |dp/dr|/B^2 = %.6g' % alpha_MHD)
print('  反径比 eps = a/R = %.6g；台基顶半径 r_ped = a - Delta_ped = %.6g m'
      % (eps, r_ped))
print('  自举电流密度 j_bs = %.3g sqrt(eps)/B_theta |dp/dr| = %.6g A/m^2' % (C_bs, j_bs))
print('  平均电流密度 j_avg = I_p/(pi a^2) = %.6g A/m^2' % j_avg)
print('  j_bs / j_avg = %.6g' % (j_bs / j_avg))
print('  台基自举电流 I_bs = j_bs Delta_ped 2 pi R = %.6g A = %.6g MA'
      % (I_bs, I_bs / 1e6))
print('  I_bs / I_p = %.6g' % (I_bs / I_p))

# ---------------------------------------------------------------- [AD.3]
print('')
print('[AD.3] 台基压强上限：气球模与剥离模两条边界')
p_max_no_bs = alpha_crit_no_bs * B ** 2 / (2.0 * mu0 * R * q_95 ** 2) * Delta_ped
p_max_bs = alpha_crit_bs * B ** 2 / (2.0 * mu0 * R * q_95 ** 2) * Delta_ped
V_ped = S_wall * Delta_ped
W_ped = 1.5 * 0.5 * p_ped * V_ped
print('  气球模上限（无自举 alpha_crit=%.3g）：p_max = %.6g Pa = %.6g kPa'
      % (alpha_crit_no_bs, p_max_no_bs, p_max_no_bs / 1e3))
print('  气球模上限（含自举 alpha_crit=%.3g）：p_max = %.6g Pa = %.6g kPa'
      % (alpha_crit_bs, p_max_bs, p_max_bs / 1e3))
print('  声明的 p_ped = %.6g kPa -> alpha_MHD/alpha_crit(无自举) = %.6g，'
      'alpha_MHD/alpha_crit(含自举) = %.6g'
      % (p_ped / 1e3, alpha_MHD / alpha_crit_no_bs, alpha_MHD / alpha_crit_bs))
print('  剥离模：j_bs/j_avg = %.6g；I_bs/I_p = %.6g' % (j_bs / j_avg, I_bs / I_p))
print('  台基壳层体积 V_ped = S_wall Delta_ped = %.6g m^3 = %.6g %% V'
      % (V_ped, V_ped / V_pl * 100))
print('  壳层储能（线性降到 0）W_ped = (3/2)(p_ped/2)V_ped = %.6g J = %.6g MJ'
      % (W_ped, W_ped / 1e6))
print('  W_ped / W_th = %.6g' % (W_ped / W_th))
print('  与总 beta 的对照：beta_ped = %.4g %%，芯部 beta ≈ %.4g %%'
      % (beta_ped * 100, 2 * mu0 * (p_avg * 1.5) / B ** 2 * 100))

# ---------------------------------------------------------------- [AD.4]
print('')
print('[AD.4] ELM：储能坍塌、频率与靶板能量密度')
dW_ELM = f_ELM_W * W_th
A_div = 2.0 * math.pi * R * L_div
E_den = dW_ELM / A_div
q_ELM = E_den / t_ELM
f_ELM = P_SOL / dW_ELM
T_ELM = 1.0 / f_ELM
print('  每次损失 dW_ELM = %.3g W_th = %.6g J = %.6g MJ' % (f_ELM_W, dW_ELM, dW_ELM / 1e6))
print('  偏滤器面积 A_div = 2 pi R L_div = %.6g m^2' % A_div)
print('  能量密度 = dW_ELM/A_div = %.6g J/m^2 = %.6g MJ/m^2' % (E_den, E_den / 1e6))
print('  峰值热流（t_ELM=%.3g ms）= %.6g W/m^2 = %.6g GW/m^2'
      % (t_ELM * 1e3, q_ELM, q_ELM / 1e9))
print('  与钨阈值 E_th=%.3g MJ/m^2 之比 = %.6g（余量 %.6g 倍）'
      % (E_th_W / 1e6, E_den / E_th_W, E_th_W / E_den))
print('  功率平衡定出的 ELM 频率 f_ELM = P_SOL/dW_ELM = %.6g Hz' % f_ELM)
print('  ELM 周期 T_ELM = %.6g s = %.6g ms' % (T_ELM, T_ELM * 1e3))
print('  每 1000 个 ELM 的时间 = %.6g s' % (1000 * T_ELM))
print('  与 tau_E 的对照：tau_E=%.3g s -> f_ELM tau_E = %.6g；'
      'W_ped/tau_E = %.6g MW' % (tau_E, f_ELM * tau_E, W_ped / tau_E / 1e6))
print('  ELM 通道的功率 = f_ELM dW_ELM = %.6g MW（= P_SOL，自洽）'
      % (f_ELM * dW_ELM / 1e6))

# ---------------------------------------------------------------- [AD.5]
print('')
print('[AD.5] 缓解：弹丸步调与共振磁扰动')
f_pace = 20.0
dW_pace = P_SOL / f_pace
E_pace = dW_pace / A_div
q_pace = E_pace / t_ELM
E_ok = 0.30e6
f_min = P_SOL / (E_ok * A_div)
print('  弹丸步调 f=%g Hz -> 每次 dW = P_SOL/f = %.6g MJ' % (f_pace, dW_pace / 1e6))
print('  能量密度 = %.6g MJ/m^2（未缓解 %.6g MJ/m^2 的 %.4g 倍）'
      % (E_pace / 1e6, E_den / 1e6, E_pace / E_den))
print('  峰值热流 = %.6g GW/m^2' % (q_pace / 1e9))
print('  要把能量密度压到 %.3g MJ/m^2 以下，频率须 >= %.6g Hz' % (E_ok / 1e6, f_min))
for n_r in (2, 3):
    x = math.sqrt((n_r - q_0) / (q_95 - q_0))
    print('  共振面 q=%d：r/a = %.6g -> r = %.6g m（台基顶内 %.6g cm）'
          % (n_r, x, x * a, (r_ped - x * a) * 1e2))
m_pel = N_pellet * m_atom
dn_pel = N_pellet / V_ped
print('  弹丸：%.3g 个原子 -> 质量 %.6g kg = %.6g mg' % (N_pellet, m_pel, m_pel * 1e6))
print('  全部沉积在 V_ped=%.6g m^3 -> dn = %.6g m^-3 = %.6g %% n_ped'
      % (V_ped, dn_pel, dn_pel / n_ped * 100))

# ---------------------------------------------------------------- [AD.6]
print('')
print('[AD.6] 边界与 SOL：热流通道宽度（接口到附录 N）')
lam_q_mm = 0.63 * B_theta ** -1.19
lam_q = lam_q_mm * 1e-3
lam_t = f_exp * lam_q
A_out = 2.0 * math.pi * R * lam_t
q_peak = f_out * P_SOL / A_out
print('  Eich 定标 lambda_q = 0.63 B_theta^-1.19 = %.6g mm' % lam_q_mm)
print('  靶板足迹 lambda_t = f_exp lambda_q = %.6g m = %.6g mm'
      % (lam_t, lam_t * 1e3))
print('  外侧靶板湿面积 A_out = 2 pi R lambda_t = %.6g m^2' % A_out)
print('  未缓解峰值热流 q_peak = f_out P_SOL/A_out = %.6g W/m^2 = %.6g MW/m^2'
      % (q_peak, q_peak / 1e6))
print('  与材料上限 %.3g MW/m^2 之比 = %.6g' % (q_limit / 1e6, q_peak / q_limit))
f_rad_need = 1.0 - q_limit / q_peak
print('  只靠辐射：需要辐射掉 %.6g 的 P_SOL' % f_rad_need)
for widen, frad in ((5.0, 0.90), (5.0, 0.95), (10.0, 0.90)):
    q = f_out * (1.0 - frad) * P_SOL / (A_out * widen)
    print('  展宽 %.3g 倍 + 辐射 %.3g -> q = %.6g MW/m^2' % (widen, frad, q / 1e6))
print('  与附录 X 的壁面韧致辐射热流 55.1547 kW/m^2 之比 = %.6g'
      % (q_peak / 55.1547e3))

# ---------------------------------------------------------------- [AD.7]
print('')
print('[AD.7] 六道练习解答')
# 1) 辅助加热 50 MW
P_SOL1 = P_alpha + 50.0e6 - P_rad_core
print('  1) P_aux=50 MW: P_SOL = %.6g MW，P_SOL/P_th = %.6g，余量 %.6g MW'
      % (P_SOL1 / 1e6, P_SOL1 / P_th, (P_SOL1 - P_th) / 1e6))
# 2) 密度减半
P_th2 = C_th * 0.5 ** ex_n * B ** ex_B * S_wall ** ex_S * 1e6
print('  2) n_e=5e19 m^-3: P_th = %.6g MW（%.6g 倍），f_G = %.6g'
      % (P_th2 / 1e6, P_th2 / P_th, 0.5 * f_G))
# 3) 磁场 4 T
P_th3 = C_th * 1.0 ** ex_n * 4.0 ** ex_B * S_wall ** ex_S * 1e6
print('  3) B=4 T: P_th = %.6g MW（%.6g 倍）；B_theta 不变（I_p 不变）'
      % (P_th3 / 1e6, P_th3 / P_th))
# 4) ELM 3% W_th
dW4 = 0.03 * W_th
f4 = P_SOL / dW4
E4 = dW4 / A_div
print('  4) dW=3%% W_th: dW = %.6g MJ，f = %.6g Hz，能量密度 = %.6g MJ/m^2'
      % (dW4 / 1e6, f4, E4 / 1e6))
# 5) n=2 共振面（已算，重列）
x2 = math.sqrt((2 - q_0) / (q_95 - q_0))
print('  5) q=2 共振面 r/a = %.6g，r = %.6g m，距离最外磁面 %.6g cm'
      % (x2, x2 * a, (a - x2 * a) * 1e2))
# 6) 台基 0.8 n_e / 2 keV
p6 = 2.0 * 0.8 * n_e * 2.0 * keV
beta_p6 = 2.0 * mu0 * p6 / B_theta ** 2
D6 = C_w * beta_p6 ** ex_w * R
alpha6 = 2.0 * mu0 * R * q_95 ** 2 * (p6 / D6) / B ** 2
print('  6) n_ped=0.8 n_e, T_ped=2 keV: p_ped = %.6g kPa，beta_p,ped = %.6g，'
      'Delta_ped = %.6g cm，alpha_MHD = %.6g'
      % (p6 / 1e3, beta_p6, D6 * 1e2, alpha6))
print('')
print('=' * 78)
print('以上即附录 AD 正文与门禁 P36 可以引用的全部数值。')
print('=' * 78)
