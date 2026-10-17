# -*- coding: utf-8 -*-
"""附录 AI（惯性约束内爆与点火）的唯一数值来源。

只 print，不写文件。用法：
    python tmp\\appAI_calc.py | Out-File -Encoding utf8 tmp\\appAI_calc.txt

五类数值分三级：
  现算级 —— 本脚本打印（正文与门禁只允许引用本文件的数）
  声明级 —— 脚本顶部输入或明确标注为量级/文献值的常数
  文献级 —— 只作定性对照，不参与任何数值断言
"""
import math

# ---------- CODATA 2018 ----------
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
u = 1.66053906660e-27
hbar = 1.054571817e-34
kB = 1.380649e-23
eps0 = 8.8541878128e-12
sigma = 5.670374419e-8
a_rad = 4.0 * sigma / c
keV = 1.602176634e-16
MeV = 1.602176634e-13

print('=' * 72)
print('附录 AI 惯性约束内爆与点火：数值来源（tmp/appAI_calc.py）')
print('=' * 72)

# ---------- DECLARED 输入 ----------
E_drv = 2.0e6            # J    驱动器（3ω 激光）能量
eta_ce = 0.85            # -    激光 -> 黑腔 X 射线转换效率（声明）
tau_drv = 4.0e-9         # s    驱动脉宽
R_hohl = 1.5e-3          # m    黑腔半径
L_hohl = 5.0e-3          # m    黑腔长度
rho_ice = 255.0          # kg/m^3  DT 冰初始密度
m_fuel = 2.82743e-6      # kg   燃料质量（与附录 Y 的 2.8274 mg 同一点）
d_abl = 1.5e-4           # m    烧蚀层厚度
rho_CH = 1050.0          # kg/m^3  CH 烧蚀层密度
v_ex = 3.30e5            # m/s  烧蚀排出速度（声明为量级值）
rho_f = 200.0e3          # kg/m^3 压缩后燃料密度（与附录 Y 同点）
rho_hot = 30.0e3         # kg/m^3 热斑密度（声明）
T_hot_keV = 5.0          # keV   热斑温度（声明）
r_hot = 4.0e-5           # m     热斑半径（声明）
rhoR_alpha = 0.3         # g/cm^2 3.5 MeV α 在 DT 中的射程（声明，文献量级）
q_DT = 3.396e14          # J/kg  D-T 比能（与附录 Y 同一点）
m_DT = 2.5 * u           # kg    平均离子质量
A_T = 0.85               # -     RT 阿特伍德数（声明）
g_acc = 1.0e14           # m/s^2 内爆加速度（声明为量级值）
v_abl = 3.0e4            # m/s   烧蚀速度（声明为量级值）
delta0 = 20.0e-9         # m     初始表面粗糙度（声明）
A21_E_alpha = 3.5 * MeV  # J     α 粒子能量
Z_eff = 1.0              # -     DT 的有效电荷数

print('\n[DECLARED] 输入')
print('  E_drv = %.6g J (= %.4g MJ), eta_ce = %.6g, tau_drv = %.6g s'
      % (E_drv, E_drv / 1e6, eta_ce, tau_drv))
print('  黑腔 R = %.6g mm, L = %.6g mm' % (R_hohl * 1e3, L_hohl * 1e3))
print('  DT 冰 rho = %.6g kg/m^3, 燃料质量 = %.6g kg (= %.6g mg)'
      % (rho_ice, m_fuel, m_fuel * 1e6))
print('  烧蚀层厚度 = %.6g um, rho_CH = %.6g kg/m^3' % (d_abl * 1e6, rho_CH))
print('  v_ex = %.6g m/s = %.6g km/s（声明量级）' % (v_ex, v_ex / 1e3))
print('  压缩后 rho = %.6g kg/m^3；热斑 rho = %.6g kg/m^3, T = %.6g keV, r = %.6g um'
      % (rho_f, rho_hot, T_hot_keV, r_hot * 1e6))
print('  α 射程 rhoR = %.6g g/cm^2（声明，文献量级）' % rhoR_alpha)
print('  RT: A_T = %.6g, g = %.6g m/s^2, v_abl = %.6g m/s, delta0 = %.6g nm'
      % (A_T, g_acc, v_abl, delta0 * 1e9))

# ================= AI.1 驱动：黑腔里的辐射温度 =================
print('\n' + '-' * 72)
print('[AI.1] 驱动：黑腔里的辐射温度')
print('-' * 72)
A_hohl = 2.0 * math.pi * R_hohl * L_hohl + 2.0 * math.pi * R_hohl ** 2
V_hohl = math.pi * R_hohl ** 2 * L_hohl
E_X = eta_ce * E_drv
F_X = E_X / (A_hohl * tau_drv)
T_rad_K = (F_X / sigma) ** 0.25
T_rad_eV = T_rad_K * kB / e
P_rad = a_rad * T_rad_K ** 4 / 3.0
u_rad = a_rad * T_rad_K ** 4
E_cav = u_rad * V_hohl
F_SB = sigma * T_rad_K ** 4
print('  黑腔内表面积 A = %.6g m^2 (= %.6g cm^2)' % (A_hohl, A_hohl * 1e4))
print('  黑腔体积 V = %.6g m^3 (= %.6g cm^3)' % (V_hohl, V_hohl * 1e6))
print('  X 射线能量 E_X = %.6g J (= %.4g MJ)' % (E_X, E_X / 1e6))
print('  壁面平均通量 F_X = %.6g W/m^2' % F_X)
print('  辐射温度 T_rad = %.6g K = %.6g eV' % (T_rad_K, T_rad_eV))
print('  辐射压 P_rad = a T^4/3 = %.6g Pa = %.6g Mbar' % (P_rad, P_rad / 1e11))
print('  辐射能密度 u = a T^4 = %.6g J/m^3' % u_rad)
print('  腔内辐射能量 uV = %.6g J（远小于 E_X：辐射场是流动的，不是储存的）' % E_cav)
print('  sigma T^4 = %.6g W/m^2（与 F_X 的口径一致性检查）' % F_SB)

# ================= AI.2 内爆：火箭模型与滞止 =================
print('\n' + '-' * 72)
print('[AI.2] 内爆：火箭模型与滞止')
print('-' * 72)
R_ice = (3.0 * m_fuel / (4.0 * math.pi * rho_ice)) ** (1.0 / 3.0)
R_caps = R_ice + d_abl
V_abl = 4.0 / 3.0 * math.pi * (R_caps ** 3 - R_ice ** 3)
m_abl = rho_CH * V_abl
m0 = m_fuel + m_abl
ratio = m0 / m_fuel
v_imp = v_ex * math.log(ratio)
E_kin = 0.5 * m_fuel * v_imp ** 2
eta_hyd = E_kin / E_drv
t_imp = R_ice / v_imp
a_acc = v_imp ** 2 / (2.0 * R_ice)
P_stag = rho_f * v_imp ** 2
R_f = (3.0 * m_fuel / (4.0 * math.pi * rho_f)) ** (1.0 / 3.0)
CR_r = R_ice / R_f
CR_v = rho_f / rho_ice
print('  燃料冰球半径 R_ice = %.6g m = %.6g mm' % (R_ice, R_ice * 1e3))
print('  靶外半径 R_caps = %.6g m = %.6g mm' % (R_caps, R_caps * 1e3))
print('  烧蚀层质量 m_abl = %.6g kg = %.6g mg' % (m_abl, m_abl * 1e6))
print('  初始质量 m0 = %.6g kg = %.6g mg' % (m0, m0 * 1e6))
print('  质量比 m0/m_f = %.6g，ln = %.6g' % (ratio, math.log(ratio)))
print('  内爆速度 v_imp = %.6g m/s = %.6g km/s' % (v_imp, v_imp / 1e3))
print('  内爆动能 E_kin = %.6g J = %.6g kJ' % (E_kin, E_kin / 1e3))
print('  流体力学效率 E_kin/E_drv = %.6g (= %.6g %%)' % (eta_hyd, eta_hyd * 100))
print('  内爆时间 t_imp = %.6g s = %.6g ns' % (t_imp, t_imp * 1e9))
print('  平均加速度 a = %.6g m/s^2 = %.6g 倍 g0' % (a_acc, a_acc / 9.8))
print('  滞止压力估计 rho v^2 = %.6g Pa = %.6g Gbar' % (P_stag, P_stag / 1e14))
print('  压缩后燃料半径 R_f = %.6g m = %.6g um' % (R_f, R_f * 1e6))
print('  半径压缩比 R_ice/R_f = %.6g；密度压缩比 = %.6g' % (CR_r, CR_v))
print('  压缩后 rhoR = %.6g g/cm^2' % (rho_f * R_f / 10.0))

# ================= AI.3 绝热与熵 =================
print('\n' + '-' * 72)
print('[AI.3] 绝热与熵：为什么压缩要「冷」')
print('-' * 72)
n_i_cold = rho_f / m_DT
n_e_cold = n_i_cold * Z_eff
P_F_cold = (3.0 * math.pi ** 2) ** (2.0 / 3.0) * hbar ** 2 * n_e_cold ** (5.0 / 3.0) / (5.0 * me)
n_i_hot = rho_hot / m_DT
n_e_hot = n_i_hot * Z_eff
P_F_hot = (3.0 * math.pi ** 2) ** (2.0 / 3.0) * hbar ** 2 * n_e_hot ** (5.0 / 3.0) / (5.0 * me)
T_cold_eV = 100.0
P_th_cold = n_e_cold * (T_cold_eV * e)
P_th_hot = n_e_hot * (T_hot_keV * keV)
alpha_cold = P_th_cold / P_F_cold
alpha_hot = P_th_hot / P_F_hot
ratio_shock = (5.0 / 3.0 + 1.0) / (5.0 / 3.0 - 1.0)
print('  压缩燃料 n_i = %.6g m^-3，n_e = %.6g m^-3' % (n_i_cold, n_e_cold))
print('  费米压 P_F(rho=200 g/cm^3) = %.6g Pa = %.6g Gbar' % (P_F_cold, P_F_cold / 1e14))
print('  冷燃料热压（100 eV）= %.6g Pa，绝热因子 alpha = %.6g' % (P_th_cold, alpha_cold))
print('  热斑 n_e = %.6g m^-3，P_F = %.6g Pa，热压（5 keV）= %.6g Pa'
      % (n_e_hot, P_F_hot, P_th_hot))
print('  热斑绝热因子 alpha = %.6g（非简并，理想气体）' % alpha_hot)
print('  强激波压缩比极限 (gamma+1)/(gamma-1) = %.6g' % ratio_shock)

# ================= AI.4 流体不稳定性 =================
print('\n' + '-' * 72)
print('[AI.4] 流体不稳定性：瑞利-泰勒与烧蚀稳定化')
print('-' * 72)
lam = 1.0e-5
k = 2.0 * math.pi / lam
gamma_RT = math.sqrt(A_T * k * g_acc)
tau_RT = 1.0 / gamma_RT
N_efold = gamma_RT * t_imp
k_c = A_T * g_acc / v_abl ** 2
lam_c = 2.0 * math.pi / k_c
k_max = A_T * g_acc / (2.0 * v_abl ** 2)
lam_max = 2.0 * math.pi / k_max
gamma_max = A_T * g_acc / (2.0 * v_abl)
delta_f = delta0 * math.exp(gamma_max * t_imp)
frac_R = delta_f / R_f
print('  取 lambda = %.6g um：k = %.6g m^-1' % (lam * 1e6, k))
print('  经典 RT 增长率 gamma = sqrt(A k g) = %.6g s^-1，特征时间 %.6g ns'
      % (gamma_RT, tau_RT * 1e9))
print('  内爆期间 e 倍增长数 = %.6g（无稳定化则必然破碎）' % N_efold)
print('  烧蚀截止波数 k_c = %.6g m^-1 -> lambda_c = %.6g um' % (k_c, lam_c * 1e6))
print('  最不稳定 k_max = %.6g m^-1 -> lambda_max = %.6g um，gamma_max = %.6g s^-1'
      % (k_max, lam_max * 1e6, gamma_max))
print('  最危险模的 e 倍增长数（用 gamma_max）= %.6g' % (gamma_max * t_imp))
print('  初始 %.6g nm 的粗糙度 -> 末期 %.6g m = %.6g um（占最终半径 %.6g %%）'
      % (delta0 * 1e9, delta_f, delta_f * 1e6, frac_R * 100))
delta_allow = 5.0e-6 / math.exp(gamma_max * t_imp)
print('  若允许末期扰动不超过 5 um，则初始粗糙度上限 = %.6g m = %.6g nm'
      % (delta_allow, delta_allow * 1e9))

# ================= AI.5 热斑与点火 =================
print('\n' + '-' * 72)
print('[AI.5] 热斑与点火：alpha 自加热对轫致损失')
print('-' * 72)
V_hot = 4.0 / 3.0 * math.pi * r_hot ** 3
rhoR_hot = rho_hot * r_hot / 10.0
E_hot = 1.5 * P_th_hot * V_hot
c_s_hot = math.sqrt(2.0 * T_hot_keV * keV / m_DT)
tau_dwell = r_hot / c_s_hot
n_D = n_i_hot / 2.0
n_T = n_D
E_alpha = 3.5 * MeV


def sv_model(T_keV):
    return 1.1e-24 * T_keV ** 2


def P_alpha(nD, nT, T_keV):
    return nD * nT * sv_model(T_keV) * E_alpha


def P_brems(n_e, T_keV):
    return 1.69e-38 * Z_eff * n_e ** 2 * math.sqrt(T_keV * 1000.0)


P_a5 = P_alpha(n_D, n_T, T_hot_keV)
P_b5 = P_brems(n_e_hot, T_hot_keV)
E_a_dwell = P_a5 * V_hot * tau_dwell
E_b_dwell = P_b5 * V_hot * tau_dwell
print('  热斑体积 V = %.6g m^3；rhoR = %.6g g/cm^2' % (V_hot, rhoR_hot))
print('  热斑热压 P = %.6g Pa = %.6g Gbar' % (P_th_hot, P_th_hot / 1e14))
print('  热斑热能 E_hot = 1.5 P V = %.6g J = %.6g kJ' % (E_hot, E_hot / 1e3))
print('  热斑声速 c_s = %.6g m/s = %.6g km/s' % (c_s_hot, c_s_hot / 1e3))
print('  约束（滞止）时间 r/c_s = %.6g s = %.6g ps' % (tau_dwell, tau_dwell * 1e12))
print('  <sigma v>(5 keV) = %.6g m^3/s（本书模型 1.1e-24 T^2）' % sv_model(T_hot_keV))
print('  alpha 功率密度 = %.6g W/m^3；轫致功率密度 = %.6g W/m^3' % (P_a5, P_b5))
print('  比值 P_alpha/P_brems(5 keV) = %.6g' % (P_a5 / P_b5))
print('  一个滞止时间里 alpha 沉积能 = %.6g J = %.6g kJ（与 E_hot 比值 %.6g）'
      % (E_a_dwell, E_a_dwell / 1e3, E_a_dwell / E_hot))
print('  同期轫致损失能 = %.6g J = %.6g kJ' % (E_b_dwell, E_b_dwell / 1e3))
print('  温度扫描（alpha 与轫致的功率比）：')
for T in (1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0):
    pa = P_alpha(n_D, n_T, T)
    pb = P_brems(n_e_hot, T)
    print('    T = %.6g keV: P_a/P_br = %.6g' % (T, pa / pb))
print('  注意：本书 <sigma v> ~ 1.1e-24 T^2 在低 T 端偏高，'
      '模型给出的交叉温度约 2.4 keV，而文献的理想点火温度约 4.3 keV（文献级对照）')
print('  α 射程：rhoR = %.6g g/cm^2 -> 在 %.6g g/cm^3 中 %.6g um；在 %.6g g/cm^3 中 %.6g um'
      % (rhoR_alpha, 200.0, rhoR_alpha / 200.0 * 1e4, 30.0, rhoR_alpha / 30.0 * 1e4))

# ================= AI.6 燃烧与增益 =================
print('\n' + '-' * 72)
print('[AI.6] 燃烧与增益：把三本账合起来')
print('-' * 72)
t_burn = R_f / c_s_hot
E_fus = m_fuel * q_DT
f_burn = 3.0 / (3.0 + 6.0)
E_burn = f_burn * E_fus
G = E_burn / E_drv
G_wall = G * 0.10
E_alpha_tot = 0.1989 * E_burn
print('  燃料声波穿越时间 t_burn = R_f/c_s = %.6g s = %.6g ns（燃耗的时间上限）'
      % (t_burn, t_burn * 1e9))
print('  燃料总聚变能 E_fus = %.6g J = %.6g MJ' % (E_fus, E_fus / 1e6))
print('  燃耗份额 f = rhoR/(rhoR+6) = %.6g (= %.6g %%)' % (f_burn, f_burn * 100))
print('  实际释放 E_burn = %.6g J = %.6g MJ' % (E_burn, E_burn / 1e6))
print('  靶增益 G = E_burn/E_drv = %.6g' % G)
print('  驱动器效率 10%% 时的墙插增益 = %.6g' % G_wall)
print('  alpha 携带的能量份额 3.5/17.6 = %.6g -> E_alpha = %.4g MJ'
      % (0.1989, E_alpha_tot / 1e6))
print('  文献对照（不参与断言）：NIF 2022 年 2.05 MJ 激光能量 -> 约 3.15 MJ 聚变产额，靶增益约 1.54')

# ================= AI.7 练习解答 =================
print('\n' + '-' * 72)
print('[AI.7 练习解答]')
print('-' * 72)
# 1 驱动器降到 1 MJ
E_X1 = eta_ce * 1.0e6
F1 = E_X1 / (A_hohl * tau_drv)
T1 = (F1 / sigma) ** 0.25
T1_eV = T1 * kB / e
print('  1) E_drv = 1 MJ: E_X = %.6g J, F = %.6g W/m^2, T_rad = %.6g eV（比 2 MJ 低 %.6g 倍）'
      % (E_X1, F1, T1_eV, T_rad_eV / T1_eV))
# 2 烧蚀层 200 um
d_abl2 = 2.0e-4
R_caps2 = R_ice + d_abl2
V_abl2 = 4.0 / 3.0 * math.pi * (R_caps2 ** 3 - R_ice ** 3)
m_abl2 = rho_CH * V_abl2
ratio2 = (m_fuel + m_abl2) / m_fuel
v_imp2 = v_ex * math.log(ratio2)
E_kin2 = 0.5 * m_fuel * v_imp2 ** 2
print('  2) 烧蚀层 200 um: m_abl = %.6g mg, m0/m_f = %.6g, v_imp = %.6g km/s, '
      'E_kin = %.6g kJ, eta_hyd = %.6g %%'
      % (m_abl2 * 1e6, ratio2, v_imp2 / 1e3, E_kin2 / 1e3, E_kin2 / E_drv * 100))
# 3 v_ex = 250 km/s
v_ex3 = 2.5e5
v_imp3 = v_ex3 * math.log(ratio)
E_kin3 = 0.5 * m_fuel * v_imp3 ** 2
print('  3) v_ex = 250 km/s: v_imp = %.6g km/s, E_kin = %.6g kJ, eta_hyd = %.6g %% '
      '（动能降 %.6g 倍）'
      % (v_imp3 / 1e3, E_kin3 / 1e3, E_kin3 / E_drv * 100, E_kin / E_kin3))
# 4 粗糙度容限
d4 = 5.0e-6 / math.exp(gamma_max * t_imp)
print('  4) 允许末期 5 um: 初始粗糙度上限 = %.6g nm；20 nm 的初始粗糙度会长到 %.6g um'
      % (d4 * 1e9, delta_f * 1e6))
# 5 热斑半径翻倍
r_hot5 = 8.0e-5
V_hot5 = 4.0 / 3.0 * math.pi * r_hot5 ** 3
E_hot5 = 1.5 * P_th_hot * V_hot5
rhoR_hot5 = rho_hot * r_hot5 / 10.0
tau5 = r_hot5 / c_s_hot
print('  5) 热斑 r = 80 um（rho、T 不变）: rhoR = %.6g g/cm^2（仍差 0.3 阈值 %.6g 倍）、'
      'E_hot = %.6g kJ、滞止时间 %.6g ps、alpha 沉积 %.6g kJ'
      % (rhoR_hot5, 0.3 / rhoR_hot5, E_hot5 / 1e3, tau5 * 1e12, P_a5 * V_hot5 * tau5 / 1e3))
# 6 增益 100 需要多少 rhoR
G6 = 100.0
f6 = G6 * E_drv / E_fus
rhoR6 = 6.0 * f6 / (1.0 - f6)
print('  6) 要 G = 100: f = %.6g (= %.6g %%), rhoR = %.6g g/cm^2（附录 Y 的 3 g/cm^2 已有余量）'
      % (f6, f6 * 100, rhoR6))

print('\n' + '=' * 72)
print('以上为附录 AI 可以引用的全部数值（其余只能定性引用）。')
print('=' * 72)
