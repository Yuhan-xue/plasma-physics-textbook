# -*- coding: utf-8 -*-
r"""附录 AL（非线性等离子体物理）数值唯一来源：纯 print，不写任何文件。

用法：python tmp\appAL_calc.py | Out-File -Encoding utf8 tmp\appAL_calc.txt
工作点与附录 M/S/U/V/W/X/Y/Z/AD/AE/AI/AJ/AK 同一台（D1）：
  n_e = 1e20 m^-3, T_e = 9.194 keV, B = 5.3 T, a = 2.0 m, R = 6.2 m, m_i = 2.5 u
外加：
  ECRH：f0 = 140 GHz, P_RF = 1.0 MW, 光斑半径 r = 2.0 cm
  激光（ICF 冕区）：λ_L = 0.35 μm, n_e,corona = 1e26 m^-3, T_e,corona = 1 keV
"""
import math

# ---------------- 常数（CODATA 2018） ----------------
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
u = 1.66053906660e-27
eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6
kB = 1.380649e-23
hbar = 1.054571817e-34
keV = 1e3 * e

# ---------------- 声明输入（DECLARED） ----------------
n_e = 1e20
T_e_eV = 9194.0
B0 = 5.3
a = 2.0
R0 = 6.2
m_i = 2.5 * mp                    # 与附录 AE 的 c_s 口径一致：m_i = 2.5 m_p

f0 = 140e9                      # ECRH 频率
P_RF = 1.0e6                    # ECRH 功率
r_spot = 2.0e-2                 # ECRH 光斑半径
lam_L = 0.35e-6                 # 三倍频激光波长
n_corona = 1e26                 # ICF 冕区电子密度
T_corona_eV = 1.0e3             # ICF 冕区电子温度
P_laser = 5.0e14                # NIF 量级总功率（W，= 500 TW）
E_drv = 2.0e6                   # 驱动器能量（J）
COEF_PC = 17.4                  # 自聚焦临界功率系数 GW（声明量级，与附录 Y 同口径）
lam_De_ref = 7.12806e-5         # 与附录 AE 同一德拜长度（m）

T_e = T_e_eV * e                # J
v_the = math.sqrt(T_e / me)
w_pe = math.sqrt(n_e * e * e / (eps0 * me))
w_ce = e * B0 / me
c_s = math.sqrt(T_e / (m_i / 2.5 * 2.5))  # 见下：用附录 AE 的 c_s = sqrt(T_e/m_i)
c_s = math.sqrt(T_e / m_i)

print('=' * 78)
print('附录 AL 数值（非线性等离子体物理）——唯一数值来源 tmp/appAL_calc.py')
print('=' * 78)
print('[DECLARED] n_e=%.6g m^-3  T_e=%.6g keV  B=%.6g T  a=%.6g m  R=%.6g m  m_i=%.6g m_p'
      % (n_e, T_e_eV / 1e3, B0, a, R0, m_i / mp))
print('[DECLARED] ECRH: f0=%.6g GHz  P=%.6g MW  r=%.6g cm'
      % (f0 / 1e9, P_RF / 1e6, r_spot * 1e2))
print('[DECLARED] 激光: lambda=%.6g um  n_e,corona=%.6g m^-3  T_e,corona=%.6g keV  P_laser=%.6g TW'
      % (lam_L * 1e6, n_corona, T_corona_eV / 1e3, P_laser / 1e12))
print('[DECLARED] 自聚焦系数 %.6g (n_c/n_e) GW；德拜长度基准与附录 AE 同' % COEF_PC)
print()

# ================= AL.1 非线性从哪来：三个判据 =================
print('[AL.1] 非线性从哪来：三个无量纲数')
w0 = 2 * math.pi * f0
A_spot = math.pi * r_spot ** 2
S_RF = P_RF / A_spot
E_RF = math.sqrt(2 * S_RF / (eps0 * c))
v_osc = e * E_RF / (me * w0)
print('  ECRH 角频率 omega0 = %.6g rad/s' % w0)
print('  光斑面积 A = %.6g m^2 = %.6g cm^2' % (A_spot, A_spot * 1e4))
print('  功率密度 S = %.6g W/m^2 = %.6g MW/cm^2' % (S_RF, S_RF / 1e10))
print('  微波电场 E_RF = %.6g V/m = %.6g kV/cm' % (E_RF, E_RF / 1e5))
print('  电子抖动速度 v_osc = eE/(m omega) = %.6g m/s = %.6g c' % (v_osc, v_osc / c))
print('  热速度 v_the = sqrt(T_e/m_e) = %.6g m/s' % v_the)
print('  判据一：v_osc/v_the = %.6g（平方 %.6g）' % (v_osc / v_the, (v_osc / v_the) ** 2))
U_p_RF = me * v_osc ** 2 / 4
print('  有质动力势 U_p = m v_osc^2/4 = %.6g J = %.6g eV' % (U_p_RF, U_p_RF / e))
print('  U_p/T_e = %.6g' % (U_p_RF / T_e))
# 激光
I_L = 1e15 * 1e4             # W/m^2（1e15 W/cm^2）
I_L2 = 1e16 * 1e4
for tag, I in (('1e15 W/cm^2', I_L), ('1e16 W/cm^2', I_L2)):
    a0 = math.sqrt((I / 1e4) * (lam_L * 1e6) ** 2 / 1.37e18)
    Up = 9.33e-14 * (I / 1e4) * (lam_L * 1e6) ** 2
    print('  激光 %s: a0 = %.6g, U_p = %.6g eV, U_p/T_e(corona 1 keV) = %.6g'
          % (tag, a0, Up, Up / T_corona_eV))
n_c_035 = w0 ** 2  # placeholder to avoid unused warning
n_c_laser = me * eps0 * (2 * math.pi * c / lam_L) ** 2 / e ** 2
n_c_140 = me * eps0 * w0 ** 2 / e ** 2
print('  140 GHz 截止密度 n_c = %.6g m^-3；n_e/n_c = %.6g' % (n_c_140, n_e / n_c_140))
print('  0.35 um 截止密度 n_c = %.6g m^-3；n_corona/n_c = %.6g'
      % (n_c_laser, n_corona / n_c_laser))
print()

# ================= AL.2 孤子：色散与非线性平衡 =================
print('[AL.2] 孤子：色散与非线性平衡')
M_sonic = 1.1
A_sol = 3.0 * (M_sonic - 1.0)          # KdV：δn/n0 = 3(M-1)
width = 2.0 * lam_De_ref / math.sqrt(A_sol)
v_sol = M_sonic * c_s
print('  离子声马赫数 M = %.6g -> 相对振幅 dn/n0 = 3(M-1) = %.6g（%.6g %%）'
      % (M_sonic, A_sol, A_sol * 100))
print('  声速 c_s = sqrt(T_e/m_i) = %.6g m/s = %.6g km/s' % (c_s, c_s / 1e3))
print('  孤子宽度 Delta = 2 lambda_De/sqrt(A) = %.6g m = %.6g um'
      % (width, width * 1e6))
print('  宽度/德拜长度 = %.6g' % (width / lam_De_ref))
print('  孤子速度 = M c_s = %.6g m/s = %.6g km/s' % (v_sol, v_sol / 1e3))
print('  特征时间 Delta/c_s = %.6g s = %.6g ps' % (width / c_s, width / c_s * 1e12))
# 有质动力密度坑
def dent(Ef):
    W = eps0 * Ef ** 2 / 2
    return W / (4 * n_e * T_e)

E_probe = 1e5
print('  E = %.6g V/m 时相对密度坑 dn/n = W/(4 nT) = %.6g' % (E_probe, dent(E_probe)))
E_1pct = math.sqrt(0.01 * 4 * n_e * T_e * 2 / eps0)
S_1pct = eps0 * c * E_1pct ** 2 / 2
print('  要挖出 1%% 密度坑：E = %.6g V/m，功率密度 S = %.6g W/m^2 = %.6g W/cm^2'
      % (E_1pct, S_1pct, S_1pct / 1e4))
print('  对应 2 cm 光斑的功率 P = %.6g W = %.6g GW' % (S_1pct * A_spot, S_1pct * A_spot / 1e9))
print()

# ================= AL.3 参量不稳定性：三波共振与阈值 =================
print('[AL.3] 参量不稳定性：三波共振与阈值')
wk0_over_c = w0 / c
k0_O = wk0_over_c * math.sqrt(1 - n_e / n_c_140)
print('  O 模波数 k0 = (omega0/c) sqrt(1-n_e/n_c) = %.6g m^-1' % k0_O)
G1_over_w1 = 1.0e-3
G2_over_w2 = 1.0e-2
thr_ratio = math.sqrt(4 * G1_over_w1 * G2_over_w2)
print('  阈值条件 (v_osc/v_the)^2 > 4 (Gamma1/omega1)(Gamma2/omega2)')
print('  取 Gamma1/omega1 = %.6g、Gamma2/omega2 = %.6g -> 阈值 v_osc/v_the = %.6g'
      % (G1_over_w1, G2_over_w2, thr_ratio))
print('  工作点 v_osc/v_the = %.6g -> 低于阈值 %.6g 倍'
      % (v_osc / v_the, thr_ratio / (v_osc / v_the)))
print('  阈值功率 P_th = %.6g MW（工作点 %.6g MW）' % (P_RF * (thr_ratio / (v_osc / v_the)) ** 2 / 1e6, P_RF / 1e6))
# 激光侧：以同一套阻尼比的阈值强度
v_the_cor = math.sqrt(T_corona_eV * e / me)
w_L = 2 * math.pi * c / lam_L
v_osc_thr = thr_ratio * v_the_cor
I_thr = (v_osc_thr * me * w_L / e) ** 2 * eps0 * c / 2
print('  冕区 v_the = %.6g m/s；同一阻尼比下阈值 v_osc = %.6g m/s' % (v_the_cor, v_osc_thr))
print('  对应激光强度阈值 I_thr = %.6g W/m^2 = %.6g W/cm^2' % (I_thr, I_thr / 1e4))
# 后向受激拉曼散射（SRS）匹配：k_L = k0 + k_s，k_s = (w0 - w_L)/c * sqrt(1-n/n_c)
w0_L = 2 * math.pi * c / lam_L
w_pe_cor = math.sqrt(n_corona * e * e / (eps0 * me))
k0_cor = (2 * math.pi / lam_L) * math.sqrt(1 - n_corona / n_c_laser)
lam_De_cor = math.sqrt(eps0 * T_corona_eV * e / (n_corona * e * e))


def wL_cor(k):
    return math.sqrt(w_pe_cor ** 2 + 3 * k ** 2 * v_the_cor ** 2)


def h_srs(k):
    return k - k0_cor - (w0_L - wL_cor(k)) / c * math.sqrt(1 - n_corona / n_c_laser)


lo2, hi2 = 1.0e5, 1.0e9
for _ in range(300):
    mid = 0.5 * (lo2 + hi2)
    if h_srs(lo2) * h_srs(mid) <= 0:
        hi2 = mid
    else:
        lo2 = mid
k_L = 0.5 * (lo2 + hi2)
w_Lc = wL_cor(k_L)
k_s = k_L - k0_cor
w_s = w0_L - w_Lc
print('  冕区泵浦：omega0 = %.6g rad/s，k0 = %.6g m^-1（O 模）' % (w0_L, k0_cor))
print('  后向 SRS 匹配解：朗缪尔波 k_L = %.6g m^-1（lambda_L = %.6g nm），omega_L = %.6g rad/s'
      % (k_L, 2 * math.pi / k_L * 1e9, w_Lc))
print('  k_L lambda_De = %.6g（冕区 lambda_De = %.6g nm）；omega_L/omega0 = %.6g'
      % (k_L * lam_De_cor, lam_De_cor * 1e9, w_Lc / w0_L))
print('  散射光 k_s = %.6g m^-1（lambda_s = %.6g nm），omega_s = %.6g rad/s（f_s = %.6g THz）'
      % (k_s, 2 * math.pi / k_s * 1e9, w_s, w_s / 2 / math.pi / 1e12))
print('  校验 k_s c/sqrt(1-n/n_c) - omega_s = %.6g' %
      (k_s * c / math.sqrt(1 - n_corona / n_c_laser) - w_s))
print()

# ================= AL.4 波崩塌与自聚焦 =================
print('[AL.4] 波崩塌与自聚焦')
Pc_tok = COEF_PC * (n_c_140 / n_e) * 1e9
Pc_cor = COEF_PC * (n_c_laser / n_corona) * 1e9
print('  临界功率 P_c = %.6g (n_c/n_e) GW' % COEF_PC)
print('  托卡马克（140 GHz，n_e=1e20）：n_c/n_e = %.6g -> P_c = %.6g W = %.6g GW'
      % (n_c_140 / n_e, Pc_tok, Pc_tok / 1e9))
print('  相对工作点 ECRH 功率的倍数 = %.6g' % (Pc_tok / P_RF))
print('  ICF 冕区（0.35 um，n_e=1e26）：n_c/n_e = %.6g -> P_c = %.6g W = %.6g TW'
      % (n_c_laser / n_corona, Pc_cor, Pc_cor / 1e12))
print('  与驱动器峰值功率 %.6g TW 之比 = %.6g' % (P_laser / 1e12, P_laser / Pc_cor))
print('  单束（192 束）功率 %.6g TW 与 P_c 之比 = %.6g'
      % (P_laser / 192 / 1e12, P_laser / 192 / Pc_cor))
lam_De_cor = math.sqrt(eps0 * T_corona_eV * e / (n_corona * e * e))
print('  冕区德拜长度 = %.6g m = %.6g nm' % (lam_De_cor, lam_De_cor * 1e9))
E_cor_1pct = math.sqrt(0.01 * 4 * n_corona * T_corona_eV * e * 2 / eps0)
print('  冕区挖 1%% 密度坑所需电场 = %.6g V/m，强度 = %.6g W/cm^2'
      % (E_cor_1pct, (eps0 * c * E_cor_1pct ** 2 / 2) / 1e4))
print('  caviton 尺寸（约 10 lambda_De）核心 = %.6g m = %.6g um'
      % (10 * lam_De_ref, 10 * lam_De_ref * 1e6))
print()

# ================= AL.5 粒子俘获与非线性朗道阻尼 =================
print('[AL.5] 粒子俘获：什么时候非线性跑赢阻尼')
k03 = 0.3 / lam_De_ref
w03 = w_pe * math.sqrt(1 + 3 * 0.3 ** 2)
gam03 = 0.0200206 * w03
print('  k lambda_De = 0.3 -> k = %.6g m^-1, omega = %.6g rad/s, gamma = %.6g s^-1'
      % (k03, w03, gam03))
print('  线性阻尼时间 1/gamma = %.6g s = %.6g ps' % (1 / gam03, 1 / gam03 * 1e12))
for Ew in (1.0e3, 1.0e4):
    w_b = k03 * math.sqrt(e * Ew / me)
    print('  E = %.6g V/m -> 俘获频率 omega_b = %.6g rad/s, 俘获时间 2pi/omega_b = %.6g ps'
          % (Ew, w_b, 2 * math.pi / w_b * 1e12))
    print('      俘获速度宽度 2 omega_b/k = %.6g m/s（= %.6g v_the），omega_b/gamma = %.6g'
          % (2 * w_b / k03, (2 * w_b / k03) / v_the, w_b / gam03))
E_crit = me * gam03 ** 2 / (k03 ** 2 * e)
print('  非线性与线性相等的场强 E* = m gamma^2/(k^2 e) = %.6g V/m' % E_crit)
k05 = 0.5 / lam_De_ref
w05 = w_pe * math.sqrt(1 + 3 * 0.5 ** 2)
gam05 = 0.151387 * w05
E_crit5 = me * gam05 ** 2 / (k05 ** 2 * e)
print('  k lambda_De = 0.5 -> k = %.6g m^-1, gamma = %.6g s^-1, E* = %.6g V/m'
      % (k05, gam05, E_crit5))
print('  E* 之比（0.5 对 0.3）= %.6g' % (E_crit5 / E_crit))
print()

# ================= AL.6 相干结构：双电层与相空间空洞 =================
print('[AL.6] 相干结构：双电层、相空间空洞')
L_dl = 100 * lam_De_ref
drop = T_e_eV
print('  双电层：厚度 100 lambda_De = %.6g m = %.6g mm' % (L_dl, L_dl * 1e3))
print('  电势降 ~ T_e/e = %.6g V -> 平均电场 = %.6g V/m = %.6g kV/cm'
      % (drop, drop / L_dl, drop / L_dl / 1e5))
Ew = 1e4
w_b = k03 * math.sqrt(e * Ew / me)
print('  相空间空洞：E = %.6g V/m 时宽度 Delta v = %.6g m/s = %.6g v_the'
      % (Ew, 2 * w_b / k03, (2 * w_b / k03) / v_the))
print('  空洞寿命 ~ 2pi/omega_b = %.6g ps（与 %.6g ps 的线性阻尼时间可比）'
      % (2 * math.pi / w_b * 1e12, 1 / gam03 * 1e12))
print('  崩塌 caviton 电场上限（波破裂）：E_WB = m_e c omega_pe/e = %.6g V/m = %.6g MV/m（核心 n_e=1e20）'
      % (me * c * w_pe / e, me * c * w_pe / e / 1e6))
print()

# ================= AL.7 练习解答 =================
print('[AL.7] 六道练习解答')
# 1) ECRH 2 MW
S2 = 2 * P_RF / A_spot
E2 = math.sqrt(2 * S2 / (eps0 * c))
v2 = e * E2 / (me * w0)
print('  1) P=2 MW: S=%.6g W/m^2, E=%.6g V/m, v_osc/v_the=%.6g（较 1 MW 涨 sqrt2 倍）'
      % (S2, E2, v2 / v_the))
print('     仍低于阈值：阈值功率 = %.6g MW' % (P_RF * (thr_ratio / (v_osc / v_the)) ** 2 / 1e6))
# 2) 激光 1e16
a0_16 = math.sqrt(1e16 * 0.35 ** 2 / 1.37e18)
Up_16 = 9.33e-14 * 1e16 * 0.35 ** 2
print('  2) I=1e16 W/cm^2, lambda=0.35 um: a0=%.6g, U_p=%.6g eV, U_p/T_e(1 keV)=%.6g'
      % (a0_16, Up_16, Up_16 / T_corona_eV))
# 3) 冕区密度 1e27
Pc_cor2 = COEF_PC * (n_c_laser / 1e27) * 1e9
print('  3) n_e,corona=1e27 m^-3: n_c/n_e=%.6g, P_c=%.6g W=%.6g GW, 与 500 TW 之比=%.6g'
      % (n_c_laser / 1e27, Pc_cor2, Pc_cor2 / 1e9, P_laser / Pc_cor2))
# 4) k lambda = 0.5 的 E*
print('  4) k lambda_De=0.5: E* = %.6g V/m（是 0.3 时的 %.6g 倍）' % (E_crit5, E_crit5 / E_crit))
# 5) 冕区 1% 密度坑
print('  5) 冕区 n=1e26, T=1 keV 挖 1%% 坑: E=%.6g V/m, 强度=%.6g W/cm^2, 相对 1e15 W/cm^2 = %.6g'
      % (E_cor_1pct, (eps0 * c * E_cor_1pct ** 2 / 2) / 1e4,
         ((eps0 * c * E_cor_1pct ** 2 / 2) / 1e4) / 1e15))
# 6) 后向 SRS 匹配解
print('  6) 后向 SRS：k_L=%.6g m^-1（%.6g nm），k_L lambda_De=%.6g，散射光 %.6g nm，f_s=%.6g THz'
      % (k_L, 2 * math.pi / k_L * 1e9, k_L * lam_De_cor, 2 * math.pi / k_s * 1e9, w_s / 2 / math.pi / 1e12))
print()
print('[AL.8] 激光侧派生量（供正文表格与门禁引用）')
E_L15 = math.sqrt(2 * I_L / (eps0 * c))
E_L16 = math.sqrt(2 * I_L2 / (eps0 * c))
v_L15 = e * E_L15 / (me * w0_L)
v_L16 = e * E_L16 / (me * w0_L)
print('  I = 1e15 W/cm^2 = %.6g W/m^2 时 E = %.6g V/m = %.6g GV/m'
      % (I_L, E_L15, E_L15 / 1e9))
print('       v_osc = eE/(m omega_L) = %.6g m/s = %.6g c' % (v_L15, v_L15 / c))
print('  I = 1e16 W/cm^2 时 E = %.6g V/m，v_osc = %.6g m/s = %.6g c'
      % (E_L16, v_L16, v_L16 / c))
print('  功率密度之比 S_laser/S_RF = %.6g' % (I_L / S_RF))
print('  有质动力之比 U_p(laser)/U_p(RF) = %.6g' % (114.292 / (U_p_RF / e)))
print('  激光角频率 omega_L = %.6g rad/s' % w0_L)
print()

print()
print('=' * 78)
# ---------------- [AL.9] 补充派生量（AL.7 解答引用） ----------------
P_th = P_RF * (thr_ratio / (v_osc / v_the)) ** 2
S2 = 2.0 * P_RF / A_spot
E2 = math.sqrt(2.0 * S2 / (eps0 * c))
v2 = e * E2 / (me * w0)
print('')
print('[AL.9] 补充派生量')
print('  练习 1：P=2 MW 时 E = %.6g V/m，v_osc = %.6g m/s，v_osc/v_the = %.6g'
      % (E2, v2, v2 / v_the))
print('  练习 1：阈值功率/实际功率（1 MW）= %.6g' % (P_th / P_RF))
print('  练习 3：P_c(1e27)/50 MW = %.6g' % (Pc_cor2 / 50e6))
print('  练习 2：U_p(1e16)/T_e(9194 eV) = %.6g' % (Up_16 / T_e_eV))
print('  练习 2：U_p(1e16)/T_e(1 keV) = %.6g' % (Up_16 / T_corona_eV))
# ---------------- [AL.10] 补充派生量（正文表格引用） ----------------
print('')
print('[AL.10] 补充派生量')
print('  caviton 尺寸（冕区，10 lambda_De）= %.6g m = %.6g nm'
      % (10 * lam_De_cor, 10 * lam_De_cor * 1e9))
print('  德拜长度（芯部，与附录 AE 同）= %.6g m = %.6g um'
      % (lam_De_ref, lam_De_ref * 1e6))
E_dl = (T_e_eV / (100.0 * lam_De_ref))
print('  双电层场强 T_e/(100 lambda_De) = %.6g V/m = %.6g MV/m = %.6g kV/cm'
      % (E_dl, E_dl / 1e6, E_dl / 1e5))
print('以上为可以引用的全部数值')
print('=' * 78)
