# -*- coding: utf-8 -*-
"""附录 AT 数值基线：冷等离子体色散、截止与共振、消逝与隧穿、伯恩斯坦波、加热落点。
只 print，不写文件（输出重定向到 tmp/appAT_calc.txt 作为唯一数值来源）。
工作点与附录 M/S/U/V/W/X/Y/AI/AL/AM/AO/AP/AQ/AR/AS 同口径：
n_e0 = 1e20 m^-3、T_e = 9.194 keV、B0 = 5.3 T、R0 = 6.2 m、a = 2.0 m、D-T 平均离子 m_i = 2.5 u。
"""
import math

e = 1.602176634e-19
m_e = 9.1093837015e-31
m_p = 1.67262192369e-27
m_D = 2.01410177812 * 1.66053906660e-27
eps0 = 8.8541878128e-12
c = 2.99792458e8
u = 1.66053906660e-27
keV = 1.602176634e-16

n_e0 = 1.0e20
T_e = 9.194 * keV
B0 = 5.3
R0 = 6.2
a = 2.0
m_i = 2.5 * u
Z = 1.0
f_rf = 140.0e9
f_ICRH = 80.8e6
e_per_m = e * e / (eps0 * m_e)
e_per_mi = Z * e * e / (eps0 * m_i)

omega = 2 * math.pi * f_rf
omega_ce0 = e * B0 / m_e
omega_pe0 = math.sqrt(n_e0 * e_per_m)
k0 = omega / c

print('==== 声明输入 ====')
print('n_e0 = %.10g m^-3, T_e = %.10g keV = %.10g J, B0 = %.10g T, R0 = %.10g m, a = %.10g m'
      % (n_e0, T_e / keV, T_e, B0, R0, a))
print('m_i = %.10g u（D-T 平均）, Z = %.10g, f_RF = %.10g Hz, f_ICRH = %.10g Hz'
      % (m_i / u, Z, f_rf, f_ICRH))
print('omega = %.10g rad/s（%.6g GHz）, k0 = omega/c = %.10g m^-1, lambda0 = %.10g mm'
      % (omega, f_rf / 1e9, k0, 2 * math.pi / k0 * 1e3))
print('e^2/(eps0 m_e) = %.10g, Z e^2/(eps0 m_i) = %.10g' % (e_per_m, e_per_mi))
print('密度剖面 n_e(r) = n_e0 (1 - (r/a)^2)；磁场 B(r) = B0 R0/(R0 + r)')

v_te = math.sqrt(T_e / m_e)          # 单分量宽度 sigma_ve，用于 3k^2 v_te^2 色散项
v_te_mp = math.sqrt(2.0 * T_e / m_e)  # 热速度，回旋半径用它
rho_e = v_te_mp / omega_ce0
v_ti = math.sqrt(T_e / m_i)
omega_ci0 = Z * e * B0 / m_i
rho_i = v_ti / omega_ci0


def B_of_r(r):
    return B0 * R0 / (R0 + r)


def n_e_of_r(r):
    return n_e0 * max(1.0 - (r / a) ** 2, 0.0)


def omega_ce_of_r(r):
    return e * B_of_r(r) / m_e


def omega_pe_of_r(r):
    return math.sqrt(n_e_of_r(r) * e_per_m)


def stix_RLSP(r, om):
    X = omega_pe_of_r(r) ** 2 / om ** 2
    Y = omega_ce_of_r(r) / om
    R = 1 - X / (1 + Y)
    L = 1 - X / (1 - Y) if abs(1 - Y) > 1e-14 else float('nan')
    S = 0.5 * (R + L)
    P = 1 - X
    return X, Y, R, L, S, P


def n2_X(r, om):
    X, Y, R, L, S, P = stix_RLSP(r, om)
    if S == 0 or math.isnan(S):
        return float('nan')
    return R * L / S


print()
print('==== AT.1 特征频率与冷等离子体介电张量 ====')
omega_ce = omega_ce0
omega_pe = omega_pe0
omega_pi = math.sqrt(n_e0 * e_per_mi)
omega_ci = Z * e * B0 / m_i
print('omega_pe = %.10g rad/s  f_pe = %.10g GHz' % (omega_pe, omega_pe / 2 / math.pi / 1e9))
print('omega_ce = %.10g rad/s  f_ce = %.10g GHz' % (omega_ce, omega_ce / 2 / math.pi / 1e9))
print('omega_pi = %.10g rad/s  f_pi = %.10g GHz' % (omega_pi, omega_pi / 2 / math.pi / 1e9))
print('omega_ci = %.10g rad/s  f_ci = %.10g MHz' % (omega_ci, omega_ci / 2 / math.pi / 1e6))
print('sigma_ve = sqrt(T_e/m_e) = %.10g m/s, v_te = sqrt(2T_e/m_e) = %.10g m/s'
      % (v_te, v_te_mp))
print('rho_e = v_te/omega_ce = %.10g m = %.6g um' % (rho_e, rho_e * 1e6))
print('v_ti = %.10g m/s, rho_i = %.10g m = %.6g mm' % (v_ti, rho_i, rho_i * 1e3))
print('rho_i/rho_e = %.10g' % (rho_i / rho_e))

X0 = omega_pe ** 2 / omega ** 2
Y0 = omega_ce / omega
R0v = 1 - X0 / (1 + Y0)
L0v = 1 - X0 / (1 - Y0)
S0v = 0.5 * (R0v + L0v)
D0v = 0.5 * (R0v - L0v)
P0v = 1 - X0
print('X = omega_pe^2/omega^2 = %.10g, Y = omega_ce/omega = %.10g' % (X0, Y0))
print('R = 1 - X/(1+Y) = %.10g' % R0v)
print('L = 1 - X/(1-Y) = %.10g' % L0v)
print('S = (R+L)/2 = %.10g, D = (R-L)/2 = %.10g, P = 1-X = %.10g' % (S0v, D0v, P0v))
print('n_O^2 = P = %.10g（O 模）' % P0v)
print('n_X^2 = R L / S = %.10g（X 模，垂直传播）' % (R0v * L0v / S0v))


def n2_roots(theta_deg):
    """Stix: A n^4 - B n^2 + C = 0，返回两个根（可能为负或 None）。"""
    th = math.radians(theta_deg)
    s2, c2 = math.sin(th) ** 2, math.cos(th) ** 2
    A = S0v * s2 + P0v * c2
    B = R0v * L0v * s2 + P0v * S0v * (1 + c2)
    C = P0v * R0v * L0v
    disc = B * B - 4 * A * C
    if disc < 0:
        return None, None
    return (B + math.sqrt(disc)) / (2 * A), (B - math.sqrt(disc)) / (2 * A)


for th in (0.0, 30.0, 60.0, 90.0):
    n2a, n2b = n2_roots(th)
    if n2a is None:
        print('theta = %5.1f deg: 判别式 < 0（无实数解）' % th)
    else:
        print('theta = %5.1f deg: n^2 两根 = %.10g（快支）与 %.10g（慢支）' % (th, max(n2a, n2b), min(n2a, n2b)))

print()
print('==== AT.2 截止与共振 ====')
omega_R = 0.5 * (omega_ce + math.sqrt(omega_ce ** 2 + 4 * omega_pe ** 2))
omega_L = 0.5 * (-omega_ce + math.sqrt(omega_ce ** 2 + 4 * omega_pe ** 2))
print('n_e = n_e0 时的截止频率：f_R = %.10g GHz, f_L = %.10g GHz'
      % (omega_R / 2 / math.pi / 1e9, omega_L / 2 / math.pi / 1e9))
print('轴向 f_ce = %.10g GHz, 2 f_ce = %.10g GHz'
      % (omega_ce / 2 / math.pi / 1e9, 2 * omega_ce / 2 / math.pi / 1e9))
if omega > omega_ce:
    n_R = omega * (omega - omega_ce) / e_per_m
    print('140 GHz：n_R = %.10g m^-3' % n_R)
else:
    print('140 GHz：omega < omega_ce -> 右旋（R）截止密度为负，等离子体内不存在 R 截止层')
n_L = omega * (omega + omega_ce) / e_per_m
print('140 GHz：n_L = %.10g m^-3' % n_L)
n_O = omega ** 2 / e_per_m
print('140 GHz：n_c（O 模截止密度）= %.10g m^-3' % n_O)
omega_UH = math.sqrt(omega_pe ** 2 + omega_ce ** 2)
print('f_UH = sqrt(f_pe^2 + f_ce^2) = %.10g GHz' % (omega_UH / 2 / math.pi / 1e9))
omega_LH = math.sqrt(omega_pi ** 2 / (1 + omega_pe ** 2 / omega_ce ** 2) + omega_ci ** 2)
print('f_LH（S = 0，准垂直）= %.10g GHz' % (omega_LH / 2 / math.pi / 1e9))
print('omega_pi^2/omega_ce^2 = %.10g, omega_ci^2/omega_LH^2 = %.6g'
      % (omega_pi ** 2 / omega_ce ** 2, omega_ci ** 2 / omega_LH ** 2))
r_ec = B0 * R0 / (omega * m_e / e) - R0
print('EC 共振半径（omega = omega_ce(r)）r = %.10g m = %.6g a, n_e(r) = %.10g m^-3'
      % (r_ec, r_ec / a, n_e_of_r(r_ec)))
B_2w = 0.5 * omega * m_e / e
r_2w = B0 * R0 / B_2w - R0
print('140 GHz 的二次谐波：omega = 2 omega_ce -> B = %.10g T -> r = %.10g m（器壁 r = a = 2 m）'
      % (B_2w, r_2w))


def bisect(fn, lo, hi, target, n=200):
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if fn(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


r_uh = bisect(lambda r: omega_pe_of_r(r) ** 2 + omega_ce_of_r(r) ** 2, 0.0, a, omega ** 2)
print('UH 共振半径（omega_UH(r) = omega）r = %.10g m = %.6g a, n_e(r) = %.10g m^-3'
      % (r_uh, r_uh / a, n_e_of_r(r_uh)))
print('  该处 omega_ce = %.10g rad/s, omega_pe = %.10g rad/s'
      % (omega_ce_of_r(r_uh), omega_pe_of_r(r_uh)))

print()
print('==== AT.3 消逝层与隧穿 ====')
print('X 模（垂直传播）n^2 随半径：')
rows = []
for r in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.3, 1.310497953, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.85, 2.0):
    X, Y, Rv, Lv, Sv, Pv = stix_RLSP(r, omega)
    n2 = Rv * Lv / Sv if Sv else float('nan')
    rows.append((r, X, Y, Rv, Lv, Sv, n2))
    print('  r = %7.4f m: X = %.6g Y = %.6g R = %.6g L = %.6g S = %.6g n_X^2 = %.6g'
          % (r, X, Y, Rv, Lv, Sv, n2))

# 找 R = 0 与 L = 0 的层（这就是消逝带的两端）
def R_of(r):
    return stix_RLSP(r, omega)[2]


def L_of(r):
    return stix_RLSP(r, omega)[3]


# 在 r ∈ (0, a) 上扫 R、L 的符号变化
cross_R, cross_L, cross_S = [], [], []
prev = None
grid = [i * a / 2000.0 for i in range(2001)]
for r in grid:
    X, Y, Rv, Lv, Sv, Pv = stix_RLSP(r, omega)
    cur = (Rv > 0, Lv > 0, Sv > 0)
    if prev is not None:
        if cur[0] != prev[0]:
            cross_R.append(r)
        if cur[1] != prev[1]:
            cross_L.append(r)
        if cur[2] != prev[2]:
            cross_S.append(r)
    prev = cur
print('R 变号处：%s' % ', '.join('%.6g m' % x for x in cross_R))
print('L 变号处：%s' % ', '.join('%.6g m' % x for x in cross_L))
print('S 变号处（上混杂共振层）：%s' % ', '.join('%.6g m' % x for x in cross_S))

# 消逝带：n_X^2 < 0 的区间
evan_ranges = []
start = None
for r in grid:
    n2 = n2_X(r, omega)
    bad = (not math.isnan(n2)) and n2 < 0
    if bad and start is None:
        start = r
    if (not bad) and start is not None:
        evan_ranges.append((start, r))
        start = None
if start is not None:
    evan_ranges.append((start, a))
print('n_X^2 < 0 的区间：%s'
      % '; '.join('%.6g -- %.6g m' % (x, y) for x, y in evan_ranges))
k0_ = k0
for (x1, x2) in evan_ranges:
    N = 4000
    acc = 0.0
    for i in range(N):
        rr = x1 + (x2 - x1) * (i + 0.5) / N
        n2 = n2_X(rr, omega)
        acc += math.sqrt(max(-n2, 0.0))
    integ = acc * (x2 - x1) / N
    tau = 2 * k0_ * integ
    print('  消逝带 %.6g -- %.6g m：宽 %.6g m，积分 ∫|n_X| dr = %.6g m，'
          '隧穿指数 2 k0 ∫|n| dr = %.6g，T = exp(-...) = %.6g'
          % (x1, x2, x2 - x1, integ, tau, math.exp(-tau)))
    n_avg = integ / (x2 - x1)
    print('    平均 |n_X| = %.6g（等效：宽 %.6g m、均匀 |n| 的势垒）' % (n_avg, x2 - x1))
    lam0_ = 2 * math.pi / k0
    print('    消逝带宽度 %.6g m = %.6g 个真空波长（lambda0 = %.6g mm）'
          % (x2 - x1, (x2 - x1) / lam0_, lam0_ * 1e3))

for label, rr in (('L 截止层', 1.578), ('UH 共振层', 1.310497953)):
    print('  %s r = %.6g m: B = %.10g T, n_e = %.10g m^-3, f_ce = %.6g GHz, f_pe = %.6g GHz, f_UH = %.6g GHz'
          % (label, rr, B_of_r(rr), n_e_of_r(rr), omega_ce_of_r(rr) / 2 / math.pi / 1e9,
             omega_pe_of_r(rr) / 2 / math.pi / 1e9,
             math.sqrt(omega_pe_of_r(rr) ** 2 + omega_ce_of_r(rr) ** 2) / 2 / math.pi / 1e9))


def f_R_layer(r):
    oc = omega_ce_of_r(r)
    if omega <= oc:
        return -1.0
    return n_e_of_r(r) - omega * (omega - oc) / e_per_m


_lo, _hi = r_ec + 0.05, a
if f_R_layer(_hi) < 0 < f_R_layer(_lo):
    r_R = bisect(lambda r: -f_R_layer(r), _lo, _hi, 0.0)
    Yr = omega_ce_of_r(r_R) / omega
    print('  R 截止层（n_e(r) = n_R(r)）r = %.10g m = %.6g a: n_e = %.10g m^-3, B = %.10g T, Y = %.6g'
          % (r_R, r_R / a, n_e_of_r(r_R), B_of_r(r_R), Yr))
    print('    该层 O-X 最优 n_par = %.6g（角度 %.6g 度）'
          % (math.sqrt(Yr / (1 + Yr)), math.degrees(math.asin(math.sqrt(Yr / (1 + Yr))))))
    print('    距回旋共振层 %.6g m，距上混杂层 %.6g m（%.6g 倍）'
          % (r_R - r_ec, r_uh - r_R, (r_uh - r_R) / (r_R - r_ec)))
else:
    print('  该剖面上没有 R 截止层：f_R_layer(a) = %.6g, f_R_layer(%.6g) = %.6g'
          % (f_R_layer(_hi), _lo, f_R_layer(_lo)))

print()
print('==== AT.4 伯恩斯坦波（动理学修正）====')
for kperp_rho in (0.1, 0.3, 0.5):
    kperp = kperp_rho / rho_e
    om = math.sqrt(omega_UH ** 2 + 3 * kperp ** 2 * v_te ** 2)
    print('  EBW k_perp rho_e = %.3g: k_perp = %.10g m^-1, f = %.10g GHz'
          '（比 f_UH 高 %.6g GHz）, lambda_perp = %.6g mm, v_g = %.6g m/s'
          % (kperp_rho, kperp, om / 2 / math.pi / 1e9,
             (om - omega_UH) / 2 / math.pi / 1e9, 2 * math.pi / kperp * 1e3,
             3 * kperp * v_te ** 2 / om))
for kperp_rho in (0.1, 0.3, 0.5):
    kperp = kperp_rho / rho_i
    om = math.sqrt(omega_LH ** 2 + 3 * kperp ** 2 * v_ti ** 2)
    print('  IBW k_perp rho_i = %.3g: k_perp = %.10g m^-1, f = %.10g GHz'
          % (kperp_rho, kperp, om / 2 / math.pi / 1e9))
kperp_1ghz = math.sqrt(((omega_UH + 2 * math.pi * 1e9) ** 2 - omega_UH ** 2) / (3 * v_te ** 2))
print('  EBW 频率上移 1 GHz 所需 k_perp = %.10g m^-1（k_perp rho_e = %.6g，lambda_perp = %.6g mm）'
      % (kperp_1ghz, kperp_1ghz * rho_e, 2 * math.pi / kperp_1ghz * 1e3))

print()
print('==== AT.5 共振条件与加热落点 ====')
f_cH = Z * e * B0 / m_p / 2 / math.pi
f_cD = Z * e * B0 / m_D / 2 / math.pi
f_cHe3 = Z * e * B0 / (3 * u) / 2 / math.pi
print('B0 = 5.3 T 时的回旋频率：H %.10g MHz, D %.10g MHz, He-3 %.10g MHz'
      % (f_cH / 1e6, f_cD / 1e6, f_cHe3 / 1e6))
print('D 二次谐波 2 f_cD = %.10g MHz；与 H 基频之比 = %.10g（理想 m_D = 2 m_p 时正好 1）'
      % (2 * f_cD / 1e6, 2 * f_cD / f_cH))
print('ICRH 取 %.6g MHz：' % (f_ICRH / 1e6))
B_H = 2 * math.pi * f_ICRH * m_p / (Z * e)
B_D2 = 2 * math.pi * f_ICRH * m_D / (2 * Z * e)
print('  H 基频共振层：B = %.10g T -> r = %.10g m' % (B_H, B0 * R0 / B_H - R0))
print('  D 二次谐波层：B = %.10g T -> r = %.10g m（两者相差 %.6g m）'
      % (B_D2, B0 * R0 / B_D2 - R0, abs((B0 * R0 / B_D2 - R0) - (B0 * R0 / B_H - R0))))

print()
print('==== AT.6 色散面与 CMA 位置 ====')
for f in (100e9, 140e9, 170e9, 200e9):
    om = 2 * math.pi * f
    X = omega_pe ** 2 / om ** 2
    Y = omega_ce / om
    Rv = 1 - X / (1 + Y)
    Lv = 1 - X / (1 - Y)
    Sv = 0.5 * (Rv + Lv)
    Pv = 1 - X
    print('  f = %6.1f GHz: X = %.6g, Y = %.6g, n_O^2 = %.6g, n_X^2 = %.6g'
          % (f / 1e9, X, Y, Pv, Rv * Lv / Sv))

print()
print('==== AT.7 估算练习解答 ====')
f1 = 170e9
om1 = 2 * math.pi * f1
n_uh_170 = (om1 ** 2 - omega_ce ** 2) / e_per_m
print('[练1] 170 GHz：n_UH = (omega^2 - omega_ce^2)/(e^2/eps0 m_e) = %.10g m^-3（%.6g n_e0）'
      % (n_uh_170, n_uh_170 / n_e0))
r_uh_170 = bisect(lambda r: omega_pe_of_r(r) ** 2 + omega_ce_of_r(r) ** 2, 0.0, a, om1 ** 2)
print('      170 GHz 的 UH 半径 r = %.10g m = %.6g a' % (r_uh_170, r_uh_170 / a))
n_L_170 = om1 * (om1 + omega_ce) / e_per_m
print('[练2] 170 GHz 的 L 截止密度 n_L = %.10g m^-3' % n_L_170)
B_70 = 2 * math.pi * 70e9 * m_e / e
print('[练3] 70 GHz 对应 B = %.10g T -> r = %.10g m（在等离子体外）' % (B_70, B0 * R0 / B_70 - R0))
print('[练4] 目标透过率 T = 1%%（2 k0 |n| L = ln 100 = %.6g）：' % math.log(100.0))
for nbar in (0.3, 0.5, 1.0):
    L_req = math.log(100.0) / (2 * k0 * nbar)
    print('      |n| = %.3g 时允许的消逝层厚度 L = %.10g m = %.6g mm' % (nbar, L_req, L_req * 1e3))
print('      本装置的消逝带指数 2 k0 ∫|n| dr = %.6g（T = exp(-2792.83)，完全截止）' % (2 * k0 * 0.475913))
Y_cut = 1.0
n_par_opt = math.sqrt(Y_cut / (1 + Y_cut))
print('[练5] O-X 最优 n_par = sqrt(Y/(1+Y))（Y = 1，即在 R 截止层）= %.10g（入射角 %.6g 度）'
      % (n_par_opt, math.degrees(math.asin(n_par_opt))))
Y_axis = omega_ce / omega
n_par_axis = math.sqrt(Y_axis / (1 + Y_axis))
print('      轴上 Y = %.6g -> n_par = %.6g（角度 %.6g 度）'
      % (Y_axis, n_par_axis, math.degrees(math.asin(n_par_axis))))
print('[练6] k0 L_n（L_n = 0.02 m）= %.6g；波长 lambda0 = %.6g mm'
      % (k0 * 0.02, 2 * math.pi / k0 * 1e3))
