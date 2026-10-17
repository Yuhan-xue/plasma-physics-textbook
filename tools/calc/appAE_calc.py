# -*- coding: utf-8 -*-
"""附录 AE（动理学与回旋动理学）的数值唯一来源 —— 只 print，不写文件。

用法：
    python tmp\\appAE_calc.py | Out-File -Encoding utf8 tmp\\appAE_calc.txt

工作点沿用 D1（与附录 M/S/U/V/W/X/Y/Z/AA/AB/AC/AD 同口径）：
    n_e = 1e20 m^-3, T_e = T_i = 9.194 keV, B = 5.3 T, R = 6.2 m, a = 2.0 m,
    m_i = 2.5 m_p, Z_eff = 2, q95 = 3.5, tau_E = 3.7 s, lnLambda = 17

速度约定（与既有附录一致，务必写进正文通知读者）：
    v_the = sqrt(T_e/m_e)          —— 与附录 AC 的 v_the 同口径（1/e 速度）
    v_Ti  = sqrt(2 T_i/m_i)        —— 与附录 AD 的 rho_i 同口径（最概然速度）
    c_s   = sqrt(T_e/m_i)          —— 附录 AD 的离子声速；T_i = T_e 时 v_Ti = sqrt(2) c_s
"""
import math

# ---------------------------------------------------------------- 物理常数
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
eps0 = 8.8541878128e-12
keV = 1e3 * e
MeV = 1e6 * e

# ---------------------------------------------------------------- DECLARED
n_e = 1.0e20          # m^-3
T_e_eV = 9.194e3      # eV
T_i_eV = 9.194e3      # eV
B0 = 5.3              # T
R0 = 6.2              # m
a0 = 2.0              # m
m_i = 2.5 * mp        # kg
Z_eff = 2.0
q95 = 3.5
tau_E = 3.7           # s（附录 M/S/Z 口径）
lnL = 17.0
tau_ei = 8.92144e-5   # s（与附录 AC 打印值一致：3.44e5 T^1.5/(n_cm3 lnL Z_eff)）
tau_ii = tau_ei * math.sqrt(m_i / me)   # s（量级估计：离子-离子动量交换）
m_alpha = 4.0 * mp
E_alpha = 3.5 * MeV

T_e = T_e_eV * e
T_i = T_i_eV * e

print('=' * 78)
print('附录 AE 数值来源（工作点 D1；只 print 不写文件）')
print('=' * 78)
print('[DECLARED] 输入与声明值：')
print('  n_e = %.6g m^-3，T_e = T_i = %.6g eV，B = %.6g T' % (n_e, T_e_eV, B0))
print('  R = %.6g m，a = %.6g m，m_i = %.6g m_p，Z_eff = %.6g' % (R0, a0, m_i / mp, Z_eff))
print('  q95 = %.6g，tau_E = %.6g s，lnLambda = %.6g' % (q95, tau_E, lnL))
print('  tau_ei = %.6g s（与附录 AC 同口径），tau_ii = tau_ei sqrt(m_i/m_e) = %.6g s'
      % (tau_ei, tau_ii))
print('  alpha 出生能 = %.6g MeV，m_alpha = %.6g m_p' % (E_alpha / MeV, m_alpha / mp))
print('  速度约定：v_the = sqrt(T_e/m_e)、v_Ti = sqrt(2T_i/m_i)、c_s = sqrt(T_e/m_i)')
print('')

# ---------------------------------------------------------------- AE.1
v_the = math.sqrt(T_e / me)
v_Ti = math.sqrt(2.0 * T_i / m_i)
c_s = math.sqrt(T_e / m_i)
lam_De = math.sqrt(eps0 * T_e / (n_e * e * e))
N_D = n_e * (4.0 / 3.0) * math.pi * lam_De ** 3
w_pe = math.sqrt(n_e * e * e / (eps0 * me))
w_ce = e * B0 / me
w_ci = e * B0 / m_i
rho_e = v_the / w_ce
rho_i = v_Ti / w_ci

print('-' * 78)
print('[AE.1] 速度空间与长度刻度')
print('  v_the = sqrt(T_e/m_e) = %.6g m/s = %.6g c' % (v_the, v_the / c))
print('  v_Ti  = sqrt(2T_i/m_i) = %.6g m/s = %.6g km/s' % (v_Ti, v_Ti / 1e3))
print('  c_s   = sqrt(T_e/m_i)  = %.6g m/s = %.6g km/s（= v_Ti/sqrt(2)）' % (c_s, c_s / 1e3))
print('  v_Ti/v_the = %.6g' % (v_Ti / v_the))
print('  lambda_De = %.6g m = %.6g um' % (lam_De, lam_De * 1e6))
print('  德拜球内粒子数 (4pi/3) n lambda_De^3 = %.6g' % N_D)
print('  德拜立方 n lambda_De^3 = %.6g（与附录 AC 的 N_D 同口径）' % (n_e * lam_De ** 3))
print('  omega_pe = %.6g rad/s = %.6g Hz' % (w_pe, w_pe / (2 * math.pi)))
print('  omega_ce = %.6g rad/s = %.6g Hz' % (w_ce, w_ce / (2 * math.pi)))
print('  omega_ci = %.6g rad/s = %.6g Hz' % (w_ci, w_ci / (2 * math.pi)))
print('  rho_e = v_the/omega_ce = %.6g m = %.6g um' % (rho_e, rho_e * 1e6))
print('  rho_i = v_Ti/omega_ci  = %.6g m = %.6g mm' % (rho_i, rho_i * 1e3))
print('  rho_i/rho_e = %.6g（= sqrt(2 m_i/m_e)）' % (rho_i / rho_e))
print('  a/rho_i = %.6g，R/rho_i = %.6g' % (a0 / rho_i, R0 / rho_i))
print('  (rho_i/a)^2 = %.6g' % ((rho_i / a0) ** 2))
print('  rho_i/a = %.6g' % (rho_i / a0))
print('  a/lambda_De = %.6g，lambda_De/a = %.6g' % (a0 / lam_De, lam_De / a0))
print('')

# ---------------------------------------------------------------- AE.2
nu_ei = 1.0 / tau_ei
lam_mfp = v_the * tau_ei
L_conn = 2.0 * math.pi * R0
eps = a0 / R0
eps15 = eps ** 1.5
nu_star = nu_ei * q95 * R0 / (eps15 * v_the)
f_trap = math.sqrt(2.0 * eps / (1.0 + eps))
delta_b = rho_i * q95 / math.sqrt(eps)
w_b = v_Ti * math.sqrt(eps) / (q95 * R0)
f_b = w_b / (2.0 * math.pi)
T_b = 1.0 / f_b

print('-' * 78)
print('[AE.2] 无碰撞判据、捕获粒子与香蕉轨道')
print('  nu_ei = 1/tau_ei = %.6g Hz' % nu_ei)
print('  omega_pe tau_ei = %.6g，omega_ce tau_ei = %.6g，omega_ci tau_ei = %.6g'
      % (w_pe * tau_ei, w_ce * tau_ei, w_ci * tau_ei))
print('  lambda_mfp = v_the tau_ei = %.6g m = %.6g km' % (lam_mfp, lam_mfp / 1e3))
print('  连接长度 L = 2 pi R = %.6g m' % L_conn)
print('  lambda_mfp/L = %.6g（热电子绕环 %0.1f 圈才碰一次）' % (lam_mfp / L_conn, lam_mfp / L_conn))
print('  epsilon = a/R = %.6g，epsilon^1.5 = %.6g，q95 R = %.6g m' % (eps, eps15, q95 * R0))
print('  碰撞率 nu* = nu_ei q R/(epsilon^1.5 v_the) = %.6g（< 1 即香蕉区）' % nu_star)
print('  捕获粒子比例 f_t = sqrt(2eps/(1+eps)) = %.6g%%' % (f_trap * 100))
print('  1/sqrt(epsilon) = %.6g（新经典对经典的抬高因子）' % (1.0 / math.sqrt(eps)))
print('  香蕉宽度 Delta_b = rho_i q95/sqrt(eps) = %.6g m = %.6g cm' % (delta_b, delta_b * 100))
print('  Delta_b/rho_i = %.6g，Delta_b/a = %.6g%%' % (delta_b / rho_i, delta_b / a0 * 100))
print('  弹跳频率 omega_b = v_Ti sqrt(eps)/(q95 R) = %.6g /s' % w_b)
print('  弹跳频率 f_b = %.6g Hz，弹跳周期 T_b = %.6g s = %.6g us'
      % (f_b, T_b, T_b * 1e6))
print('')

# ---------------------------------------------------------------- AE.3
def landau_langmuir(kld):
    return -math.sqrt(math.pi / 8.0) * kld ** -3 * math.exp(-1.0 / (2.0 * kld ** 2) - 1.5)

def landau_ion_acoustic(te_over_ti, mi_over_me):
    first = math.sqrt(1.0 / mi_over_me)
    second = (1.0 / te_over_ti) ** 1.5 * math.exp(-te_over_ti / 2.0 - 1.5)
    return -math.sqrt(math.pi / 8.0) * (first + second)

kld = [0.1, 0.3, 0.5, 1.0]
print('-' * 78)
print('[AE.3] 朗道阻尼：谁在共振，阻尼多大')
for x in kld:
    print('  Langmuir k lambda_De = %.6g : gamma/omega = %.6g（%.6g%%）'
          % (x, landau_langmuir(x), -landau_langmuir(x) * 100))
k03 = 0.3
w_k03 = w_pe * math.sqrt(1.0 + 3.0 * k03 ** 2)
k_03 = k03 / lam_De
lam_k03 = 2.0 * math.pi / k_03
v_ph = w_k03 / k_03
ratio = v_ph / v_the
print('  k lambda_De = 0.3 时：omega = %.6g rad/s，k = %.6g m^-1，波长 = %.6g m = %.6g mm'
      % (w_k03, k_03, lam_k03, lam_k03 * 1e3))
print('  相速度 v_ph = %.6g m/s = %.6g v_the' % (v_ph, ratio))
print('  共振因子 exp(-v_ph^2/(2 v_the^2)) = %.6g' % math.exp(-ratio ** 2 / 2.0))
tail_ph = math.erfc(ratio / math.sqrt(2.0)) + math.sqrt(2.0 / math.pi) * ratio * math.exp(-ratio ** 2 / 2.0)
print('  速度超过 v_ph 的电子比例 = %.6g（朗道阻尼只吃这一小撮）' % tail_ph)
for r in [1.0, 10.0]:
    g = landau_ion_acoustic(r, m_i / me)
    print('  离子声 T_e/T_i = %.6g : gamma/omega = %.6g（%.6g%%）' % (r, g, -g * 100))
k_ia = 0.1 / lam_De
lam_ia = 2.0 * math.pi / k_ia
w_ia = k_ia * c_s
print('  离子声 k lambda_De = 0.1 ：k = %.6g m^-1，波长 = %.6g m = %.6g mm' % (k_ia, lam_ia, lam_ia * 1e3))
print('  离子声 omega = k c_s = %.6g rad/s，f = %.6g Hz = %.6g MHz'
      % (w_ia, w_ia / (2 * math.pi), w_ia / (2 * math.pi) / 1e6))
print('')

# ---------------------------------------------------------------- AE.4
def bessel_i0(x):
    s, term, k = 1.0, 1.0, 0
    while True:
        k += 1
        term *= (x * x / 4.0) / (k * k)
        s += term
        if term < 1e-16 or k > 200:
            return s

def gamma0(b):
    return bessel_i0(b) * math.exp(-b)

b_list = [0.09, 1.0, 4.0]
print('-' * 78)
print('[AE.4] 回旋动理学尺度：k_perp rho_i 与有限拉莫尔半径因子')
kb_rho_i = 0.3
k_perp = kb_rho_i / rho_i
lam_perp = 2.0 * math.pi / k_perp
print('  k_perp rho_i = %.6g : k_perp = %.6g m^-1，lambda_perp = %.6g m = %.6g cm'
      % (kb_rho_i, k_perp, lam_perp, lam_perp * 100))
print('  同一波长下 k_perp rho_e = %.6g（电子几乎感觉不到这个尺度）' % (k_perp * rho_e))
for b in b_list:
    print('  b = (k_perp rho_i)^2 = %.6g : Gamma_0(b) = I_0(b) e^-b = %.6g，1-Gamma_0 = %.6g'
          % (b, gamma0(b), 1.0 - gamma0(b)))
chi_gB = rho_i ** 2 * c_s / a0
chi_B = T_e / (16.0 * e * B0)
chi_eff = a0 ** 2 / (2.0 * tau_E)
print('  回旋玻姆扩散系数 chi_gB = rho_i^2 c_s/a = %.6g m^2/s' % chi_gB)
print('  玻姆扩散系数 chi_B = T_e/(16 e B) = %.6g m^2/s' % chi_B)
print('  chi_B/chi_gB = %.6g' % (chi_B / chi_gB))
print('  由 tau_E = a^2/(2 chi) 反推 chi_eff = %.6g m^2/s' % chi_eff)
print('  chi_eff/chi_gB = %.6g（该工作点比回旋玻姆「好」%.6g 倍）'
      % (chi_eff / chi_gB, chi_gB / chi_eff))
print('  chi_B/chi_eff = %.6g' % (chi_B / chi_eff))
print('')

# ---------------------------------------------------------------- AE.5
def tail_speed(u):
    """3 维麦克斯韦分布中 v > u v_th 的比例，v_th = sqrt(T/m)。"""
    return math.erfc(u) + (2.0 / math.sqrt(math.pi)) * u * math.exp(-u * u)

def tail_energy(x):
    """3 维麦克斯韦分布中 E > x T 的比例。"""
    r = math.sqrt(x)
    return math.erfc(r) + (2.0 / math.sqrt(math.pi)) * r * math.exp(-x)

print('-' * 78)
print('[AE.5] 速度分布：尾部到底有多少粒子')
for u in [1.0, 2.0, 3.0, 4.0]:
    print('  v > %.6g v_the : 比例 = %.6g（%.6g%%），n_e 中 %.6g m^-3'
          % (u, tail_speed(u), tail_speed(u) * 100, n_e * tail_speed(u)))
T_eV = T_e_eV
for E_keV in [10.0, 64.0, 100.0]:
    x = E_keV * 1e3 / T_eV
    print('  E > %.6g keV = %.6g T : 比例 = %.6g' % (E_keV, x, tail_energy(x)))
ratio_alpha = math.sqrt(2.0 * E_alpha / m_alpha) / v_Ti
print('  3.5 MeV alpha：v = %.6g m/s = %.6g v_Ti' % (math.sqrt(2.0 * E_alpha / m_alpha),
                                                    ratio_alpha))
print('  把它当热尾：exp(-u^2) = %.3g（完全不可积分的尾部，必须单独处理）'
      % math.exp(-ratio_alpha ** 2))
rho_alpha = m_alpha * math.sqrt(2.0 * E_alpha / m_alpha) / (2.0 * e * B0)
print('  3.5 MeV alpha 的回旋半径 = %.6g m = %.6g cm（是 rho_i 的 %.6g 倍）'
      % (rho_alpha, rho_alpha * 100, rho_alpha / rho_i))
for Rm in [2.0, 5.0]:
    sin_tc = math.sqrt(1.0 - 1.0 / Rm)
    frac = 1.0 - math.sqrt(1.0 - 1.0 / Rm)
    print('  镜比 R_m = %.6g : 损失锥角 = %.6g 度，各向同性分布中落入损失锥的比例 = %.6g%%'
          % (Rm, math.degrees(math.asin(sin_tc)), frac * 100))
print('')

# ---------------------------------------------------------------- AE.6
Kn = lam_mfp / a0
lam_mfp_i = v_Ti * tau_ii
Kn_i = lam_mfp_i / a0
print('-' * 78)
print('[AE.6] 从动理学到流体：矩、闭合与各向异性输运')
print('  Knudsen 数（电子）Kn = lambda_mfp/a = %.6g' % Kn)
print('  离子平均自由程 lambda_mfp,i = v_Ti tau_ii = %.6g m = %.6g km' % (lam_mfp_i, lam_mfp_i / 1e3))
print('  Knudsen 数（离子）Kn_i = %.6g' % Kn_i)
print('  平行/垂直热导之比量级：kappa_perp/kappa_par ~ 1/(omega_ce tau_e)^2 = %.6g'
      % (1.0 / (w_ce * tau_ei) ** 2))
print('  离子：1/(omega_ci tau_ii)^2 = %.6g' % (1.0 / (w_ci * tau_ii) ** 2))
print('')

# ---------------------------------------------------------------- AE.7
print('=' * 78)
print('[AE.7] 练习解答（6 道）')
print('-' * 78)
# 练习 1
n1 = 1.0e19
T1 = 1.0e3 * e
Ti1 = 1.0e3 * e
lam1 = math.sqrt(eps0 * T1 / (n1 * e * e))
ND1 = n1 * (4.0 / 3.0) * math.pi * lam1 ** 3
vTi1 = math.sqrt(2.0 * Ti1 / m_i)
rho_i1 = vTi1 / w_ci
print('  练习 1  n=1e19 m^-3、T=1 keV、B=5.3 T：lambda_De = %.6g m = %.6g um，'
      % (lam1, lam1 * 1e6))
print('          德拜立方 N_D = %.6g（附录 AC 口径），德拜球 = %.6g' % (n1 * lam1 ** 3, ND1))
print('          v_Ti = %.6g m/s，rho_i = %.6g m = %.6g mm，a/rho_i = %.6g'
      % (vTi1, rho_i1, rho_i1 * 1e3, a0 / rho_i1))
# 练习 2
k2 = 1.0 / rho_i
lam2 = 2.0 * math.pi / k2
print('  练习 2  k_perp rho_i = 1：k_perp = %.6g m^-1，lambda_perp = %.6g m = %.6g cm，'
      % (k2, lam2, lam2 * 100))
print('          Gamma_0(1) = %.6g，抑制掉 %.6g%%' % (gamma0(1.0), (1.0 - gamma0(1.0)) * 100))
print('          同一波长下 k_perp rho_e = %.6g（= 1/%.6g）' % (k2 * rho_e, rho_i / rho_e))
# 练习 3
t3 = tail_speed(3.0)
print('  练习 3  v > 3 v_the：比例 = %.6g（%.6g%%），1 m^3 里 %.6g 个电子'
      % (t3, t3 * 100, n_e * t3))
# 练习 4
tau4 = tau_ei * 10.0
lam4 = v_the * tau4
print('  练习 4  n=1e19 m^-3（tau ∝ 1/n）：tau_ei = %.6g s，lambda_mfp = %.6g m = %.6g km，'
      % (tau4, lam4, lam4 / 1e3))
print('          lambda_mfp/L = %.6g' % (lam4 / L_conn))
# 练习 5
db5 = rho_i * 1.7 / math.sqrt(eps)
print('  练习 5  q = 1.7 时 Delta_b = %.6g m = %.6g cm（q = 3.5 时为 %.6g cm，比 %.6g 倍）'
      % (db5, db5 * 100, delta_b * 100, q95 / 1.7))
# 练习 6
v_alpha = math.sqrt(2.0 * E_alpha / m_alpha)
print('  练习 6  3.5 MeV alpha：v = %.6g m/s = %.6g v_Ti，热尾比例 exp(-u^2) = %.3g，'
      % (v_alpha, v_alpha / v_Ti, math.exp(-(v_alpha / v_Ti) ** 2)))
print('          rho_alpha = %.6g m = %.6g cm = %.6g rho_i'
      % (rho_alpha, rho_alpha * 100, rho_alpha / rho_i))
print('=' * 78)
print('（本文件由 tmp/appAE_calc.py 生成；正文与门禁 P37 引用的数字都必须出现在上面。）')
