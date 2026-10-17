# -*- coding: utf-8 -*-
"""附录 R（尘埃等离子体）数值唯一来源。

把附录 R 正文与门禁里出现的每一个数字都算出来，输出到 tmp/appR_calc.txt。
只做独立复算，不引用任何二手数字。
"""
import io
import math
import os
import sys

# ------------------------------------------------------------------ 常数
E = 1.602176634e-19          # C
EPS0 = 8.8541878128e-12      # F/m
KB = 1.380649e-23            # J/K
ME = 9.1093837015e-31        # kg
U = 1.66053906660e-27        # kg
M_AR = 39.948 * U            # kg
M_D = 2.013553212745 * U     # kg
M_H = 1.007825032 * U        # kg
G0 = 9.80665                 # m/s^2
AU = 1.495978707e11          # m
C_LIGHT = 2.99792458e8       # m/s
G_NEWTON = 6.67430e-11       # m^3/(kg s^2)
L_SUN = 3.828e26             # W
M_SUN = 1.98847e30           # kg
EV_J = E                     # 1 eV = E joule
COUL_NM = E / (4 * math.pi * EPS0) * 1e9   # e^2/(4 pi eps0) in eV*nm

OUT = []


def say(fmt, *a):
    OUT.append(fmt % a if a else fmt)


def debye(te_ev, ne):
    return math.sqrt(EPS0 * te_ev * EV_J / (ne * E ** 2))


def z_oml(te_ev, ti_ev, m_i):
    """解 OML 悬浮电位：exp(z) = sqrt(m_i/(tau m_e)) / (1 + z/tau)。"""
    tau = ti_ev / te_ev
    lo, hi = 1e-6, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = math.exp(mid) * (1.0 + mid / tau) - math.sqrt(m_i / (tau * ME))
        if f < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def zd_of(a_m, te_ev, z):
    return 4.0 * math.pi * EPS0 * a_m * z * te_ev / E


def v_eth(te_ev, m):
    return math.sqrt(8.0 * te_ev * EV_J / (math.pi * m))


def m_grain(a_m, rho):
    return 4.0 / 3.0 * math.pi * a_m ** 3 * rho


def delta_spacing(n_d):
    return (3.0 / (4.0 * math.pi * n_d)) ** (1.0 / 3.0)


def gamma_coul(z, delta_m, kt_ev):
    return z ** 2 * COUL_NM / (delta_m * 1e9 * kt_ev)


# ------------------------------------------------------------------ 工况
CASES = [
    dict(name='C1 刻蚀工艺腔尘埃', gas='Ar', m_i=M_AR, te=3.0, ti=300.0 * KB / EV_J,
         ne=1.0e17, a=0.5e-6, rho=2000.0, nd=1.0e12, td=0.0259),
    dict(name='C2 实验室尘埃晶体', gas='Ar', m_i=M_AR, te=3.0, ti=300.0 * KB / EV_J,
         ne=1.0e15, a=4.8e-6, rho=1510.0, nd=1.0e9, td=0.03),
    dict(name='C3 聚变边界尘埃', gas='D', m_i=M_D, te=15.0, ti=15.0,
         ne=2.0e19, a=1.0e-6, rho=2000.0, nd=1.0e10, td=1.0),
    dict(name='C4 星际尘埃', gas='H', m_i=M_H, te=1.0, ti=100.0 * KB / EV_J,
         ne=1.0e6, a=0.1e-6, rho=2000.0, nd=1.0e-6, td=0.01),
]

say('COUL_NM = e^2/(4 pi eps0) = %.6f eV*nm' % COUL_NM)
say('')

for c in CASES:
    ld = debye(c['te'], c['ne'])
    z = z_oml(c['te'], c['ti'], c['m_i'])
    phi = z * c['te']
    zd = zd_of(c['a'], c['te'], z)
    qd = zd * E
    md = m_grain(c['a'], c['rho'])
    area = 4.0 * math.pi * c['a'] ** 2
    ve = v_eth(c['te'], ME)
    ie0 = math.pi * c['a'] ** 2 * E * c['ne'] * ve
    tch = qd / ie0
    delta = delta_spacing(c['nd'])
    havnes = c['nd'] * zd / c['ne']
    ld_eff = ld / math.sqrt(1.0 + havnes)
    g_delta = gamma_coul(zd, delta, c['td'])
    g_a = gamma_coul(zd, c['a'], c['td'])
    kappa = delta / ld
    screen = (1.0 + kappa + 0.5 * kappa ** 2) * math.exp(-kappa)
    cda = math.sqrt(zd ** 2 * c['nd'] * c['ti'] * EV_J / (c['ne'] * md))
    f_da = cda / 1.0e-3
    e_need = md * G0 / qd
    f_coul = (qd ** 2 / (4.0 * math.pi * EPS0)) / delta ** 2
    f_grav = md * G0
    say('== %s (%s, T_e=%.2f eV, n_e=%.3g m^-3, a=%.2f um, rho=%.0f, n_d=%.3g m^-3)'
        % (c['name'], c['gas'], c['te'], c['ne'], c['a'] * 1e6, c['rho'], c['nd']))
    say('  lambda_D        = %.6e m = %.4f um' % (ld, ld * 1e6))
    say('  a/lambda_D      = %.6f' % (c['a'] / ld))
    say('  a/lambda_D (sci)= %.4e' % (c['a'] / ld))
    say('  m_d/m_i         = %.6e' % (md / c['m_i']))
    say('  omega_pd        = %.6e rad/s（= sqrt(n_d Z_d^2 e^2/(eps0 m_d))）'
        % math.sqrt(c['nd'] * zd ** 2 * E ** 2 / (EPS0 * md)))
    say('  z = e|phi|/Te   = %.6f' % z)
    say('  phi_d           = %.6f V' % phi)
    say('  Z_d             = %.6e e' % zd)
    say('  Q_d             = %.6e C' % qd)
    say('  m_d             = %.6e kg' % md)
    say('  4 pi a^2        = %.6e m^2' % area)
    say('  v_eth,e         = %.6e m/s' % ve)
    say('  I_e0            = %.6e A' % ie0)
    say('  t_charge=Q/I_e0 = %.6e s = %.4f ns' % (tch, tch * 1e9))
    say('  Delta (spacing) = %.6e m = %.4f um' % (delta, delta * 1e6))
    say('  P_Havnes        = %.6e' % havnes)
    say('  lambda_eff      = %.4f um' % (ld_eff * 1e6))
    say('  Gamma(Delta)    = %.6e' % g_delta)
    say('  Gamma(a)        = %.6e' % g_a)
    say('  Gamma_a/Gamma_D = %.6f' % (g_a / g_delta))
    say('  kappa=Delta/ld  = %.6f' % kappa)
    say('  Yukawa factor   = %.6f' % screen)
    say('  Gamma_eff       = %.6e' % (g_delta * screen))
    say('  c_DA            = %.6e m/s = %.4f cm/s' % (cda, cda * 100))
    say('  f_DA (1 mm)     = %.6f Hz' % f_da)
    say('  E_need=m g/Q    = %.6f V/m' % e_need)
    say('  F_coul/F_grav   = %.6f' % (f_coul / f_grav))
    say('')

# ------------------------------------------------------------------ 标度律
say('== Z_d 随半径与 T_e 的标度（Ar, T_i = 300 K）==')
z_ar = z_oml(3.0, 300.0 * KB / EV_J, M_AR)
for a_um in (0.1, 0.5, 1.0, 4.8):
    say('  a=%.1f um, T_e=3.0 eV -> Z_d=%.6e, Q_d=%.6e C'
        % (a_um, zd_of(a_um * 1e-6, 3.0, z_ar), zd_of(a_um * 1e-6, 3.0, z_ar) * E))
say('  Z_d / (a/um) at 3 eV = %.6e' % zd_of(1e-6, 3.0, z_ar))
say('  Z_d / (a/um)/(Te/eV) = %.6e' % (zd_of(1e-6, 3.0, z_ar) / 3.0))
say('  z_Ar = %.6f（与密度无关：n_i=n_e 时密度在 OML 平衡里约掉）' % z_ar)
say('')

# ------------------------------------------------------------------ 射频比较
say('== 充电时间 vs 射频周期 ==')
tch_c1 = zd_of(0.5e-6, 3.0, z_ar) * E / (math.pi * (0.5e-6) ** 2 * E * 1e17 * v_eth(3.0, ME))
for f_hz, lab in ((13.56e6, '13.56 MHz（常用射频）'), (1.0e4, '10 kHz（附录 P 的 DBD 驱动）')):
    period = 1.0 / f_hz
    say('  %s: T=%.6e s = %.4f ns, t_charge/T=%.6f' % (lab, period, period * 1e9, tch_c1 / period))
say('')

# ------------------------------------------------------------------ Havnes 临界
say('== Havnes 参数 P=1 的临界尘埃密度 n_d = n_e/Z_d ==')
for lab, ne, zd in (('刻蚀工艺腔', 1.0e17, zd_of(0.5e-6, 3.0, z_ar)),
                    ('实验室晶体', 1.0e15, zd_of(4.8e-6, 3.0, z_ar)),
                    ('聚变边界', 2.0e19, 28894.6)):
    say('  %s: n_d(P=1) = %.6e m^-3' % (lab, ne / zd))
say('')

# ------------------------------------------------------------------ 鞘层
say('== 鞘层里的悬浮（附录 P 的 100 V / 0.633 mm 鞘层）==')
s_sheath = 0.633e-3
v_sheath = 100.0
e_mean = v_sheath / s_sheath
e_wall = 5.0 / 3.0 * e_mean
say('  s=%.6f m, V_s=%.1f V -> E_mean=%.6e V/m, E_wall(Child)=(5/3)E_mean=%.6e V/m'
    % (s_sheath, v_sheath, e_mean, e_wall))
for lab, a_um, rho, te in (('0.5 um 工艺尘', 0.5, 2000.0, 3.0), ('4.8 um 晶体尘', 4.8, 1510.0, 3.0)):
    md = m_grain(a_um * 1e-6, rho)
    qd = zd_of(a_um * 1e-6, te, z_ar) * E
    e_need = md * G0 / qd
    x_over_s = (e_need * s_sheath / ((5.0 / 3.0) * v_sheath)) ** 1.5
    say('  %s: E_need=%.6f V/m, E_need/E_mean=%.6f, E_need/E_wall=%.6f, x/s=%.6e, x=%.6e m'
        % (lab, e_need, e_need / e_mean, e_need / e_wall, x_over_s, x_over_s * s_sheath))
say('')

# ------------------------------------------------------------------ 聚变尘埃
say('== 聚变边界尘埃（D 等离子体，T_e=T_i=15 eV）==')
z_d2 = z_oml(15.0, 15.0, M_D)
say('  z = %.6f, phi_d = %.6f V' % (z_d2, z_d2 * 15.0))
say('  Z_d(a=1 um) = %.6e, Q_d = %.6e C' % (zd_of(1e-6, 15.0, z_d2), zd_of(1e-6, 15.0, z_d2) * E))
m1 = m_grain(1e-6, 2000.0)
say('  m_d(1 um, 2000 kg/m3) = %.6e kg' % m1)
for mass_g in (1.0, 100.0):
    ssa = 3.0 / (2000.0 * 1e-6)
    area = mass_g * 1e-3 * ssa
    n_t = area * 1.0e21
    say('  尘埃总量 %.0f g: 比表面积 %.1f m^2/kg, 总表面积 %.4f m^2, 若面密度 1e21 m^-2 则氚 %.6e 个 = %.6f g'
        % (mass_g, ssa, area, n_t, n_t / 6.02214076e23 * 3.0160492))
say('')

# ------------------------------------------------------------------ 天体尘埃
say('== 空间尘埃：光电充电与辐射压 ==')
for phi_v in (1.0, 2.0, 5.0):
    say('  a=0.1 um, phi=+%.1f V -> Z_d = %.6e' % (phi_v, 4.0 * math.pi * EPS0 * 0.1e-6 * phi_v / E))
A_beta = 3.0 * L_SUN / (16.0 * math.pi * C_LIGHT * G_NEWTON * M_SUN)
a_crit = A_beta / 2000.0
say('  beta = A/(rho a), A = 3L/(16 pi c G M) = %.6e kg/m^2' % A_beta)
for a_um, rho in ((0.1, 2000.0), (1.0, 2000.0), (1.0, 3000.0)):
    say('  a=%.1f um, rho=%.0f -> beta(1 AU) = %.6f' % (a_um, rho, A_beta / (rho * a_um * 1e-6)))
say('  rho=2000 kg/m^3 时的临界半径 a_crit = %.6e m = %.4f um' % (a_crit, a_crit * 1e6))
say('')

# ------------------------------------------------------------------ 尘埃声波对照
say('== 尘埃声波与离子声速对照（C1）==')
c1 = CASES[0]
ld1 = debye(c1['te'], c1['ne'])
cs_ion = math.sqrt(c1['te'] * EV_J / M_AR)
zd1 = zd_of(c1['a'], c1['te'], z_ar)
md1 = m_grain(c1['a'], c1['rho'])
cda1 = math.sqrt(zd1 ** 2 * c1['nd'] * c1['ti'] * EV_J / (c1['ne'] * md1))
say('  c_s(Ar, T_e=3 eV)  = %.6f m/s' % cs_ion)
say('  c_DA               = %.6e m/s = %.4f cm/s' % (cda1, cda1 * 100))
say('  c_s/c_DA           = %.6e' % (cs_ion / cda1))
say('  omega_pi(Ar)       = %.6e rad/s' % math.sqrt(c1['ne'] * E ** 2 / (EPS0 * M_AR)))
say('  omega_pi/omega_pd  = %.6e'
    % (math.sqrt(c1['ne'] * E ** 2 / (EPS0 * M_AR))
       / math.sqrt(c1['nd'] * zd1 ** 2 * E ** 2 / (EPS0 * md1))))
for lam in (1e-3, 1e-2):
    say('  lambda=%.0f mm -> f = %.6f Hz' % (lam * 1e3, cda1 / lam))
say('')

# ------------------------------------------------------------------ 量级跨度
say('== 四个工况的量级跨度与量级判据（写正文时直接引用）==')
nes = [c['ne'] for c in CASES]
lds = [debye(c['te'], c['ne']) for c in CASES]
say('  n_e 范围 %s -> %s m^-3，跨度 %.2f 个数量级'
    % ('%.3g' % min(nes), '%.3g' % max(nes), math.log10(max(nes) / min(nes))))
say('  lambda_D 范围 %.4f um -> %.2f m，跨度 %.2f 个数量级'
    % (min(lds) * 1e6, max(lds), math.log10(max(lds) / min(lds))))
for a_um in (0.01, 10.0):
    say('  m_d/m_Ar（rho=2000, a=%.2f um）= %.4e' % (a_um, m_grain(a_um * 1e-6, 2000.0) / M_AR))
_zd_c4 = zd_of(0.1e-6, 1.0, z_oml(1.0, CASES[3]['ti'], CASES[3]['m_i']))
say('  P_Havnes: C1 %.4e（Z_d=%.1f）, C4 %.4e（Z_d=%.1f），相差 %.2f 个数量级'
    % (1e12 * zd_of(0.5e-6, 3.0, z_ar) / 1e17, zd_of(0.5e-6, 3.0, z_ar), 1e-6 * _zd_c4 / 1e6, _zd_c4,
       math.log10((1e12 * zd_of(0.5e-6, 3.0, z_ar) / 1e17) / (1e-6 * _zd_c4 / 1e6))))
say('  0.751240 g 氚 = %.6f mol（M_T = 3.0160492 g/mol）' % (0.751240 / 3.0160492))
say('')

# ------------------------------------------------------------------ 练习答案汇总
say('== 六道练习答案（与正文一致）==')
zd_05 = zd_of(0.5e-6, 3.0, z_ar)
say('  练习1: phi_d=%.6f V, Z_d=%.6e, Q_d=%.6e C' % (z_ar * 3.0, zd_05, zd_05 * E))
say('  练习2: t_charge=%.4f ns, 13.56 MHz 周期=%.4f ns, 比值=%.6f'
    % (tch_c1 * 1e9, 1e9 / 13.56e6, tch_c1 * 13.56e6))
delta1 = delta_spacing(1e12)
say('  练习3: Delta=%.4f um, Gamma(Delta)=%.6e, Gamma(a)=%.6e, 比值=%.4f'
    % (delta1 * 1e6, gamma_coul(zd_05, delta1, 0.0259), gamma_coul(zd_05, 0.5e-6, 0.0259),
       gamma_coul(zd_05, 0.5e-6, 0.0259) / gamma_coul(zd_05, delta1, 0.0259)))
say('  练习4: n_d(P=1)=%.6e m^-3, lambda_eff(1e12)=%.4f um'
    % (1e17 / zd_05, debye(3.0, 1e17) / math.sqrt(1.0 + 1e12 * zd_05 / 1e17) * 1e6))
say('  练习5: c_DA=%.4f cm/s, f(1 mm)=%.6f Hz' % (cda1 * 100, cda1 / 1e-3))
md_05 = m_grain(0.5e-6, 2000.0)
say('  练习6: E_need=%.6f V/m; beta(1 um, 2000)=%.6f; a_crit=%.4f um'
    % (md_05 * G0 / (zd_05 * E), A_beta / 2000.0 / 1e-6, a_crit * 1e6))
say('')
say('以上数字是正文与门禁的唯一来源。')

text = '\n'.join(OUT) + '\n'
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = os.path.join(root, 'tmp', 'appR_calc.txt')
io.open(path, 'w', encoding='utf-8', newline='\n').write(text)
sys.stdout.write(text)
sys.stdout.write('\n[written] %s (%d bytes)\n' % (path, len(text.encode('utf-8'))))
