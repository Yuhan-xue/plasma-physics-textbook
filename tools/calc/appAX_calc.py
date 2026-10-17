# -*- coding: utf-8 -*-
"""
附录 AX（束流物理与束流品质）的数值唯一来源。
正文 src/appendix_ax_beams.tex 与门禁 check_P56_appendix_ax_beams()
只允许引用本脚本打印出来的数字；本脚本不引用任何外部结论数值。

工作点 D1（与附录 M/S/U/V/W/AW 同一口径）：
  n_e = 1e20 m^-3、T_e = T_i = 9.194 keV、B = 5.3 T、a = 2.0 m、R = 6.2 m、
  V = 832.2 m^3、P_fus = 556.5 MW、P_alpha = 111.3 MW、lnΛ = 17、Z_eff = 2。

声明输入（中性束注入的束流参数）：
  氘束 E_b = 100 keV、总电流 I_b = 40 A、束半径 r_b = 0.15 m；
  引出小孔（beamlet）半径 r_bl = 5.0 mm、发散角 θ_bl = 0.30°、小孔数 N_bl = 900；
  引出间隙 d = 0.02 m；漂移距离 L_drift = 10 m；Twiss β* = 2.0 m、α* = 0；
  中和效率 η_neut = 0.55；气体靶线密度 n_gl = 1e20 m^-2、剥离截面 σ_strip = 1e-20 m^2；
  中性束有效电离截面 σ_eff = 1e-20 m^2（与附录 M 同口径）。
只算经典慢化与冷束极限，不含波纹/磁岛/不稳定性导致的额外损失。
"""
import math

# ---------- 物理常数（CODATA 2018） ----------
e = 1.602176634e-19
eps0 = 8.8541878128e-12
m_e = 9.1093837015e-31
amu = 1.66053906660e-27
m_D = 2.01410177812 * amu
m_i = 2.5 * amu
kB = 1.380649e-23
c = 2.99792458e8
keV = 1.0e3 * e
MeV = 1.0e6 * e

# ---------- D1 工作点 ----------
n = 1.0e20
T_keV = 9.194
T = T_keV * keV
B = 5.3
a = 2.0
R = 6.2
V = 832.2
P_fus = 556.5e6
P_alpha = 111.3e6
lnL = 17.0
Z_eff = 2.0
n_i = n

# ---------- 声明束流参数 ----------
E_b = 100.0 * keV
I_b = 40.0
r_b = 0.15
r_bl = 5.0e-3
th_bl = math.radians(0.30)
N_bl = 900
d_gap = 0.02
L_drift = 10.0
beta_star = 2.0
eta_neut = 0.55
n_gl = 1.0e20
sig_strip = 1.0e-20
sig_eff = 1.0e-20

out = []


def p(fmt, *args):
    s = fmt % args if args else fmt
    out.append(s)
    print(s)


def head(t):
    p('')
    p('=' * 72)
    p('== ' + t)
    p('=' * 72)


# ================================================================
head('0. 声明输入与工作点（正文引用的声明值）')
# ================================================================
p('E_b = %.1f keV、I_b = %.1f A、r_b = %.2f m', E_b / keV, I_b, r_b)
p('r_bl = %.1f mm、theta_bl = %.2f deg = %.6e rad、N_bl = %d', r_bl * 1e3,
  math.degrees(th_bl), th_bl, N_bl)
p('d_gap = %.3f m、L_drift = %.1f m、beta* = %.1f m', d_gap, L_drift, beta_star)
p('eta_neut = %.2f、n_gl = %.0e m^-2、sigma_strip = %.0e m^2、sigma_eff = %.0e m^2',
  eta_neut, n_gl, sig_strip, sig_eff)
p('D1: n_e = %.0e m^-3、T_e = %.3f keV、B = %.1f T、a = %.1f m、R = %.1f m、V = %.1f m^3',
  n, T_keV, B, a, R, V)
p('D1: P_fus = %.1f MW、P_alpha = %.1f MW、lnL = %.1f、Z_eff = %.1f、m_i = %.1f u',
  P_fus / 1e6, P_alpha / 1e6, lnL, Z_eff, m_i / amu)

# ================================================================
head('1. 束流运动学：速度、回旋半径、束内密度')
# ================================================================
v_b = math.sqrt(2.0 * E_b / m_D)
beta = v_b / c
gamma = 1.0 / math.sqrt(1.0 - beta * beta)
rho_c = m_D * v_b / (e * B)
f_c = e * B / (2.0 * math.pi * m_D)
A_b = math.pi * r_b * r_b
j_b = I_b / A_b
n_b = I_b / (e * v_b * A_b)
p('m_D                  = %.6e kg', m_D)
p('v_b                  = %.6e m/s = %.5f c', v_b, beta)
p('gamma                = %.8f', gamma)
p('rho_c (5.3 T)        = %.6f m = %.4f cm', rho_c, rho_c * 100)
p('f_c (5.3 T)          = %.6e Hz = %.4f MHz', f_c, f_c / 1e6)
p('A_b = pi r_b^2       = %.6f m^2', A_b)
p('j_b = I_b/A_b        = %.6f A/m^2', j_b)
p('n_b = I_b/(e v_b A_b)= %.6e m^-3', n_b)
p('n_b/n_e              = %.6e', n_b / n)
p('每小孔电流 I_b/N_bl   = %.6e A', I_b / N_bl)
A_bl = math.pi * r_bl * r_bl
p('小孔面积 pi r_bl^2    = %.6e m^2', A_bl)
p('小孔电流密度 j_bl     = %.6f A/m^2（与 j_b 同值即填充因子 1）', (I_b / N_bl) / A_bl)
p('填充因子 N_bl*(r_bl/r_b)^2 = %.6f', N_bl * (r_bl / r_b) ** 2)

# ================================================================
head('2. 发射度与亮度：束流品质的两个刻度')
# ================================================================
eps_bl = r_bl * th_bl
eps_n = beta * gamma * eps_bl
sig_x = math.sqrt(eps_bl * beta_star)
sig_p = math.sqrt(eps_bl / beta_star)
B_bright = I_b / (4.0 * math.pi ** 2 * eps_bl * eps_bl)
B_bl = (I_b / N_bl) / (4.0 * math.pi ** 2 * eps_bl * eps_bl)
T_eff = m_D * v_b * v_b * th_bl * th_bl / (2.0 * kB)
p('θ_bl                 = %.6e rad = %.4f deg', th_bl, math.degrees(th_bl))
p('ε_bl = r_bl*θ_bl     = %.6e m rad = %.4f mm mrad', eps_bl, eps_bl * 1e3)
p('ε_n = βγε_bl         = %.6e m rad', eps_n)
p('σ_x = sqrt(εβ*)      = %.6e m = %.4f mm', sig_x, sig_x * 1e3)
p("σ' = sqrt(ε/β*)      = %.6e rad = %.4f mrad", sig_p, sig_p * 1e3)
p('B = I/(4π²ε²)        = %.6e A/(m^2 rad^2)', B_bright)
p('B_bl（单小孔）        = %.6e A/(m^2 rad^2)', B_bl)
p('等效束温度 T_eff      = %.4f K = %.6f eV', T_eff, T_eff * kB / e)
p('T_eff/E_b            = %.6e（「冷束」的定量含义）', (T_eff * kB) / E_b)
p('4ε（4-rms 口径）      = %.6e m rad = %.4f mm mrad', 4.0 * eps_bl, 4.0 * eps_bl * 1e3)

# ================================================================
head('3. Twiss 与漂移：束包络长多快')
# ================================================================
beta_L = beta_star + L_drift * L_drift / beta_star
sig_L = math.sqrt(eps_bl * beta_L)
grow = sig_L / sig_x
p('β(L=10 m) = β* + L²/β* = %.6f m', beta_L)
p('σ(10 m)              = %.6e m = %.4f mm', sig_L, sig_L * 1e3)
p('束斑增长比 σ(10)/σ(0) = %.6f', grow)
p('发散角对应的几何增长   = L*θ_bl = %.6f m', L_drift * th_bl)

# ================================================================
head('4. 空间电荷与派文斯：等效电位差、Child--Langmuir 极限')
# ================================================================
K_perv = I_b / (E_b / e) ** 1.5
I_A = 4.0 * math.pi * eps0 * m_D * c ** 3 / e
K_gen = 2.0 * I_b / (I_A * beta ** 3 * gamma ** 3)
E_r = j_b / (eps0 * v_b)
dV = E_r * r_b
p('K = I/V^{3/2}        = %.6e A/V^{1.5}', K_perv)
p('I_A（阿尔芬电流）     = %.6e A = %.3f MA', I_A, I_A / 1e6)
p('K_gen = 2I/(I_A β³γ³) = %.6f', K_gen)
p('E_r = j/(ε0 v)       = %.6e V/m', E_r)
p('ΔV = E_r r_b         = %.6e V = %.4f MV', dV, dV / 1e6)
p('ΔV/V_b               = %.6f 倍', dV / (E_b / e))
j_CL = (4.0 * eps0 / 9.0) * math.sqrt(2.0 * e / m_D) * (E_b / e) ** 1.5 / d_gap ** 2
I_CL = j_CL * A_b
A_need = I_b / j_CL
r_need = math.sqrt(A_need / math.pi)
p('j_CL（d=0.02 m）      = %.6f A/m^2', j_CL)
p('I_CL（引出面 A_b）    = %.6f A（= 需求的 %.4f 倍）', I_CL, I_CL / I_b)
p('引出 40 A 所需面积     = %.6e m^2 = %.4f cm^2', A_need, A_need * 1e4)
p('对应半径              = %.6f m = %.4f cm', r_need, r_need * 100)
grow_unneut = 0.5 * (K_gen / sig_x) * L_drift ** 2
K_gen_max = 2.0 * 0.01 * sig_x / L_drift ** 2
f_need = 1.0 - K_gen_max / K_gen
p('未中和时 10 m 包络膨胀   = %.6e m（= %.4f km，完全不可用）', grow_unneut, grow_unneut / 1e3)
p('若要求漂移 10 m 内只涨 0.01 m，允许的 K_gen 上限 = %.6e', K_gen_max)
p('所需中和度 f_n ≥     = %.7f（即 %.5f%%）', f_need, f_need * 100)

# ================================================================
head('5. 中和与残余离子：束流功率账')
# ================================================================
P_b = I_b * (E_b / e)
P_neut = eta_neut * P_b
P_ion = (1.0 - eta_neut) * P_b
att = math.exp(-n_gl * sig_strip)
lam = 1.0 / (n * sig_eff)
shine = math.exp(-R / lam)
P_abs = P_neut * (1.0 - shine)
w_abs = P_abs / V
p('P_b = I_b*E_b        = %.6e W = %.4f MW', P_b, P_b / 1e6)
p('P_neut = η P_b       = %.6e W = %.4f MW', P_neut, P_neut / 1e6)
p('P_ion = (1-η) P_b    = %.6e W = %.4f MW（须由残余离子偏转板带走）', P_ion, P_ion / 1e6)
p('中和靶透过率 e^{-nσL}= %.6f', att)
p('λ = 1/(n_e σ_eff)    = %.6f m', lam)
p('穿透份额 e^{-R/λ}     = %.6e = %.6f%%', shine, shine * 100)
p('被吸收份额 1-e^{-R/λ}  = %.6e = %.6f%%', 1.0 - shine, (1.0 - shine) * 100)
p('被吸收功率            = %.6e W = %.4f MW', P_abs, P_abs / 1e6)
p('功率密度 P_abs/V      = %.4f W/m^3', w_abs)

# ================================================================
head('6. 束流在等离子体里停在哪里：临界能量与慢化时间')
# ================================================================
v_Te = math.sqrt(2.0 * T / m_e)
x_e = v_b / v_Te
G_exact = math.erf(x_e) - (2.0 * x_e / math.sqrt(math.pi)) * math.exp(-x_e * x_e)
G_asym = 4.0 * x_e ** 3 / (3.0 * math.sqrt(math.pi))
ratio_G = G_exact / G_asym
K_drag = (16.0 * math.sqrt(math.pi) / 3.0) * n * e ** 4 * lnL / \
    ((4.0 * math.pi * eps0) ** 2 * m_e * v_Te ** 3)
K_eff = K_drag * ratio_G
tau_v = m_D / K_eff
tau_E = m_D / (2.0 * K_eff)
C_i = 4.0 * math.pi * n * Z_eff * e ** 4 * lnL / ((4.0 * math.pi * eps0) ** 2 * m_i)
v_c = (C_i / K_eff) ** (1.0 / 3.0)
E_c = 0.5 * m_D * v_c * v_c
E_c_rule = 14.8 * 2.0 * T_keV
F_e = K_eff * v_b
F_i = C_i / (v_b * v_b)
f_ion = F_i / (F_i + F_e)
path = v_b * tau_v
n_orbit = path / (2.0 * math.pi * R)
p('v_Te                 = %.6e m/s', v_Te)
p('x = v_b/v_Te         = %.6f（x<<1：束流是「慢束」，阻力在线性段）', x_e)
p('G(x) 精确 / 渐近      = %.8f / %.8f，比值 %.6f', G_exact, G_asym, ratio_G)
p('K_drag               = %.6e kg/s', K_drag)
p('K_eff = K_drag*比值   = %.6e kg/s', K_eff)
p('τ_v = m_D/K_eff      = %.6f s', tau_v)
p('τ_E = τ_v/2          = %.6f s', tau_E)
p('C_i（含 Z_eff 因子）   = %.6e', C_i)
p('v_c = (C_i/K_eff)^1/3= %.6e m/s', v_c)
p('E_c = ½ m_D v_c²     = %.6e J = %.4f keV = %.4f T_e', E_c, E_c / keV, E_c / T)
p('E_c 的 14.8 A T_e 定标 = %.4f keV（与上式之比 %.4f）', E_c_rule, E_c_rule / (E_c / keV))
p('E_b = 100 keV vs E_c  : E_b/E_c = %.6f（<1：阻力以离子为主）', E_b / E_c)
p('电子阻力 F_e = K_eff v = %.6e（相对单位）', F_e)
p('离子阻力 F_i = C_i/v²  = %.6e（相对单位）', F_i)
p('离子阻力/电子阻力      = %.6f = (E_c/E_b)^{3/2}', F_i / F_e)
p('能量给离子的份额       = %.6f = %.4f%%', f_ion, f_ion * 100)
p('慢化路程 v_b τ_v      = %.6e m = %.2f km', path, path / 1e3)
p('环向周长 2πR          = %.6f m', 2.0 * math.pi * R)
p('相当于绕环            = %.2f 圈', n_orbit)
p('中性束射程 λ = %.4f m vs 慢化路程 %.4e m，相差 %.3e 倍', lam, path, path / lam)

# ================================================================
head('7. 束-等离子体耦合：冷束双流模与稳定化门槛')
# ================================================================
w_pe = math.sqrt(n * e * e / (eps0 * m_e))
w_pb = math.sqrt(n_b * e * e / (eps0 * m_D))
f_pb = w_pb / (2.0 * math.pi)
gam_ts = (math.sqrt(3.0) / 2.0) * (n_b / (2.0 * n)) ** (1.0 / 3.0) * w_pe
t_grow = 1.0 / gam_ts
lam_De = v_Te / (math.sqrt(2.0) * w_pe)
k_res = w_pe / v_b
k_lam = k_res * lam_De
gam_landau = math.sqrt(math.pi / 8.0) * (1.0 / (k_lam ** 3)) * \
    math.exp(-1.0 / (2.0 * k_lam ** 2) - 1.5)
k_dv = w_pe * th_bl
transit = R / v_b
p('ω_pe                 = %.6e rad/s', w_pe)
p('ω_pb                 = %.6e rad/s = %.6f MHz', w_pb, f_pb / 1e6)
p('γ_ts = (√3/2)(n_b/2n_e)^{1/3} ω_pe = %.6e rad/s', gam_ts)
p('冷束增长时间 1/γ_ts   = %.6e s = %.4f ns', t_grow, t_grow * 1e9)
p('λ_De                 = %.6e m', lam_De)
p('共振波数 k = ω_pe/v_b = %.6e m^-1', k_res)
p('k λ_De               = %.6f（≫1：电子Landau 阻尼很弱）', k_lam)
p('Landau 阻尼率 γ_L/ω_pe   = %.6e', gam_landau)
p('束流展宽门槛 k Δv = ω_pe θ_bl = %.6e rad/s', k_dv)
p('γ_ts / (k Δv)        = %.6f（≈1 量级：冷束极限被展宽抵掉大半）', gam_ts / k_dv)
p('穿越时间 R/v_b        = %.6e s = %.4f μs（含 %.1f 个冷束 e 折）',
  transit, transit * 1e6, transit * gam_ts)

# ================================================================
head('8. 估算练习（6 题）的解答数值')
# ================================================================
# 练习 1：80 keV 氘束的 v、β、ρ_c（B=5.3 T）
E1 = 80.0 * keV
v1 = math.sqrt(2.0 * E1 / m_D)
rho1 = m_D * v1 / (e * B)
p('[练习 1] v = %.6e m/s、β = %.6f、ρ_c = %.6f m = %.4f cm',
  v1, v1 / c, rho1, rho1 * 100)
p('[练习 1] E/E_b = %.6f = %.4f%%、ρ_c/ρ_c(AX.1) = %.6f = %.4f%%（降 %.4f%%）',
  E1 / E_b, E1 / E_b * 100, rho1 / rho_c, rho1 / rho_c * 100, (1.0 - rho1 / rho_c) * 100)
# 练习 2：小孔 r=4 mm、θ=0.25° 的 ε 与 ε_n
eps2 = 4.0e-3 * math.radians(0.25)
th2 = math.radians(0.25)
p('[练习 2] θ = %.6e rad、ε = %.6e m rad = %.7f mm mrad、ε_n = %.6e m rad、σ_x(β*=2 m) = %.4f mm',
  th2, eps2, eps2 * 1e3, beta * gamma * eps2, math.sqrt(eps2 * 2.0) * 1e3)
p('[练习 2] ε/ε(AX.2) = %.6f（降 %.4f%%）、σ_x/σ_x(AX.2) = %.6f（降 %.4f%%）、r_bl/r_bl(AX.2) = %.6f = %.4f%%、θ/θ(AX.2) = %.6f = %.4f%%',
  eps2 / eps_bl, (1.0 - eps2 / eps_bl) * 100, math.sqrt(eps2 * 2.0) / sig_x,
  (1.0 - math.sqrt(eps2 * 2.0) / sig_x) * 100,
  4.0 / 5.0, 4.0 / 5.0 * 100, 0.25 / 0.30, 0.25 / 0.30 * 100)
# 练习 3：I=25 A、r_b=0.12 m 的 j、E_r、ΔV/V
I3, r3 = 25.0, 0.12
j3 = I3 / (math.pi * r3 * r3)
E3 = j3 / (eps0 * v_b)
p('[练习 3] j = %.4f A/m^2、E_r = %.6e V/m、ΔV/V_b = %.4f',
  j3, E3, E3 * r3 / (E_b / e))
p('[练习 3] I/I_b = %.6f = %.4f%%、r_b/r_b(AX.4) = %.6f = %.4f%%',
  I3 / I_b, I3 / I_b * 100, r3 / r_b, r3 / r_b * 100)
# 练习 4：V=120 kV、d=0.025 m 的 j_CL 与 I_CL（A_b=pi*0.15^2）
V4, d4 = 120.0e3, 0.025
j4 = (4.0 * eps0 / 9.0) * math.sqrt(2.0 * e / m_D) * V4 ** 1.5 / d4 ** 2
p('[练习 4] j_CL = %.4f A/m^2、I_CL(0.15 m) = %.4f A', j4, j4 * A_b)
p('[练习 4] V/V(AX.4) = %.6f = %.4f%%、d/d(AX.4) = %.6f = %.4f%%、I_CL/I_CL(AX.4) = %.6f = %.4f%%',
  V4 / (E_b / e), V4 / (E_b / e) * 100, d4 / d_gap, d4 / d_gap * 100,
  j4 * A_b / I_CL, j4 * A_b / I_CL * 100)
# 练习 5：束能 300 keV 的 x、阻力分工与 E_c 比较
E5 = 300.0 * keV
v5 = math.sqrt(2.0 * E5 / m_D)
F5e = K_eff * v5
F5i = C_i / (v5 * v5)
p('[练习 5] v = %.6e m/s、x = v/v_Te = %.6f、F_i/F_e = %.6f、给电子份额 = %.4f%%（300 keV > E_c：先加热电子）',
  v5, v5 / v_Te, F5i / F5e, 100.0 * F5e / (F5i + F5e))
# 练习 6：η=0.45 时的中性/离子功率
P6n = 0.45 * P_b
P6i = 0.55 * P_b
p('[练习 6] η=0.45：P_neut = %.4f MW、P_ion = %.4f MW', P6n / 1e6, P6i / 1e6)
p('[练习 6] 相对 η=0.55：中性功率少 %.4f MW（%.4f%%）、离子功率多 %.4f MW（%.4f%%）',
  (P_neut - P6n) / 1e6, (P_neut - P6n) / P_neut * 100,
  (P6i - P_ion) / 1e6, (P6i - P_ion) / P_ion * 100)

with open('tmp/appAX_calc.txt', 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(out) + '\n')
print('')
print('[写出] tmp/appAX_calc.txt  共 %d 行' % (len(out) + 1))
