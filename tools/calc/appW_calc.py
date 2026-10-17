# -*- coding: utf-8 -*-
"""
附录 W（快粒子与 alpha 粒子物理）的数值唯一来源。
正文 src/appendix_w_fastparticle.tex 与门禁 check_P27_appendix_w_fastparticle()
只允许引用本脚本打印出来的数字；本脚本不引用任何外部结论数值。

工作点 D1（与附录 S/U/V 同一口径）：n=1e20 m^-3（n_D=n_T=n/2）、T_e=T_i=9.194 keV、
B=5.3 T、a=2.0 m、R=6.2 m、I_p=15 MA、P_fus=556.5 MW、P_alpha=111.3 MW、V=832.2 m^3。
简化声明：T_i=T_e（书中工作点口径）；只算经典慢化，不含波纹/磁岛/不稳定性导致的额外损失。
"""
import math

# ---------- 物理常数（CODATA 2018） ----------
e = 1.602176634e-19
eps0 = 8.8541878128e-12
m_e = 9.1093837015e-31
amu = 1.66053906660e-27
m_alpha = 4.001506179127 * amu          # 6.6446573e-27 kg
kB = 1.380649e-23
mu0 = 4.0e-7 * math.pi
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
Ip_MA = 15.0
P_fus = 556.5e6
P_alpha = 111.3e6
V = 832.2
lnL = 17.0
Z_eff = 2.0
n_i = n                     # 离子总密度（D+T）
m_i = 2.5 * amu             # 平均离子质量（D-T 各半）
m_D = 2.01410177812 * amu   # 氘核质量（算回旋半径/回旋频率时用单一species更清楚）
Z_f = 2.0                   # alpha 电荷数

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
head('1. alpha 粒子出生：能量、速度、回旋半径、出生率')
# ================================================================
E_a0 = 3.5 * MeV
v_a0 = math.sqrt(2.0 * E_a0 / m_alpha)
rho_a = m_alpha * v_a0 / (Z_f * e * B)
S_a = P_alpha / E_a0
p('E_alpha0            = %.6e J  = %.1f MeV', E_a0, E_a0 / MeV)
p('v_alpha0            = %.5e m/s = %.5f c', v_a0, v_a0 / c)
p('rho_alpha           = %.6f m = %.2f cm', rho_a, rho_a * 100)
p('S_alpha (总出生率)   = %.5e /s', S_a)
p('S_alpha/V (体出生率) = %.5e /(m^3 s)', S_a / V)
p('热离子回旋半径对比   : rho_i(D,9.194keV) = %.6e m = %.2f mm',
  math.sqrt(2.0 * m_D * T) / (e * B), math.sqrt(2.0 * m_D * T) / (e * B) * 1e3)

# ================================================================
head('2. 慢化：电子阻力（线性段）、时间尺度、临界能量')
# ================================================================
v_Te = math.sqrt(2.0 * T / m_e)
x_e = v_a0 / v_Te
def Hfun(x):
    return math.erf(x) - (2.0 * x / math.sqrt(math.pi)) * math.exp(-x * x)
H_exact = Hfun(x_e)
H_asym = 4.0 * x_e ** 3 / (3.0 * math.sqrt(math.pi))
p('v_Te (电子热速度)    = %.6e m/s', v_Te)
p('x = v_a0/v_Te        = %.6f  (x<<1 -> 线性阻力段)', x_e)
p('H(x) 精确            = %.8f', H_exact)
p('H(x) 渐近 4x^3/3sqrt(pi) = %.8f  (比值 %.4f)', H_asym, H_exact / H_asym)
# 电子阻力在 x<<1 时的线性系数 K：F = -K v
K = (16.0 * math.sqrt(math.pi) / 3.0) * n * Z_f ** 2 * e ** 4 * lnL / \
    ((4.0 * math.pi * eps0) ** 2 * m_e * v_Te ** 3)
tau_v = m_alpha / K
tau_E = m_alpha / (2.0 * K)
p('K (F = -K v)         = %.6e kg/s', K)
p('tau_v = m/K          = %.6f s  (速度 e 折时间)', tau_v)
p('tau_E = m/(2K)       = %.6f s  (能量 e 折时间 = tau_v/2)', tau_E)
# 交叉验证：Spitzer/NRL 电子-离子能量均分时间（取质子 Z=1 作参照；alpha 因
# m/Z^2 = m_p/1 与质子同量级，故可直接对比）
m_p = 1.007276466621 * amu
tau_eq_ref = (3.0 / (8.0 * math.sqrt(2.0 * math.pi))) * m_p * (4.0 * math.pi * eps0) ** 2 * \
             T ** 1.5 / (math.sqrt(m_e) * n * 1.0 * e ** 4 * lnL)
p('交叉验证 tau_eq = 3/(8 sqrt(2pi)) m_p (4 pi eps0)^2 T^1.5/(sqrt(m_e) n e^4 lnL)，T 为能量[J]')
p('                 = %.6f s（质子参照；alpha 的 m/Z^2 = m_p/1，可直接对比）', tau_eq_ref)
p('             与 tau_E = m_alpha/(2K) = %.6f s 之比 = %.3f（应 ~1，量级校验）',
  tau_E, tau_eq_ref / tau_E)
# 离子阻力（alpha 对背景 D/T 离子）：H(x_i) ~ 1
v_thi = math.sqrt(2.0 * T / m_i)
C_i = 4.0 * math.pi * n_i * Z_f ** 2 * e ** 4 * lnL / ((4.0 * math.pi * eps0) ** 2 * m_i)
vc3 = C_i / K
v_c = vc3 ** (1.0 / 3.0)
E_c = 0.5 * m_alpha * v_c ** 2
p('v_th,i               = %.6e m/s', v_thi)
('')
p('C_i                  = %.6e', C_i)
p('v_c^3 = C_i/K        = %.6e', vc3)
p('v_c (临界速度)        = %.6e m/s   (v_c/v_Te = %.5f)', v_c, v_c / v_Te)
p('E_c (临界能量)        = %.6e J = %.1f keV = %.2f T_e', E_c, E_c / keV, E_c / T)
t_therm = tau_v * math.log(v_a0 / v_thi)
p('出生->热化总时间      = tau_v*ln(v_a0/v_th,i) = %.4f s  (%.2f 个 e 折)',
  t_therm, math.log(v_a0 / v_thi))
p('慢化路径长度 v_a0*tau_v = %.5e m = %.0f km', v_a0 * tau_v, v_a0 * tau_v / 1e3)
p('环向周长 2piR        = %.4f m -> 相当于 %.0f 次环向绕行', 2 * math.pi * R,
  v_a0 * tau_v / (2 * math.pi * R))

# ================================================================
head('3. alpha 加热落在哪：捕获、香蕉宽度、轨道、局域性')
# ================================================================
eps = a / R                       # 边界处反纵横比（首次轨道损失取最坏情形）
sqrt_eps = math.sqrt(eps)
q_edge = 2.0
rho_theta = q_edge * rho_a / sqrt_eps    # 极向回旋半径 = 相对出生磁面的径向偏移
Delta_a = 2.0 * rho_theta                # 完整香蕉宽度（峰-峰）
p('eps = a/R            = %.6f', eps)
p('sqrt(eps)            = %.6f', sqrt_eps)
p('极向回旋半径 rho_th = q rho_a/sqrt(eps) (q=2) = %.6f m = %.2f cm = %.2f%% 的 a',
  rho_theta, rho_theta * 100, rho_theta / a * 100)
p('完整香蕉宽度 Delta_a = 2 rho_th = %.6f m = %.2f cm (%.2f%% 的 a)',
  Delta_a, Delta_a * 100, Delta_a / a * 100)
# 首次轨道损失：n ~ (1-rho^2), n^2 ~ (1-rho^2)^2，出生在 r > a - rho_theta 的 alpha
x_edge = (a - rho_theta) / a
def tail_frac(x0):
    def F(x):
        return x * x / 2.0 - x ** 4 / 2.0 + x ** 6 / 6.0
    return (F(1.0) - F(x0)) / F(1.0)
f_prompt = tail_frac(x_edge)
P_prompt = f_prompt * P_alpha
p('首次轨道损失边界 r/a  = %.6f', x_edge)
p('首次轨道损失份额      = %.6f = %.3f%%', f_prompt, f_prompt * 100)
p('对应功率损失          = %.4f MW', P_prompt / 1e6)

# ================================================================
head('4. 慢化分布与快 alpha 压强（数值积分）')
# ================================================================
NGRID = 400000
v1 = v_thi
v2 = v_a0
dv = (v2 - v1) / NGRID
I0 = 0.0
I1 = 0.0
I2 = 0.0
for k in range(NGRID + 1):
    vv = v1 + k * dv
    w = 0.5 if (k == 0 or k == NGRID) else 1.0
    f = 1.0 / (vv ** 3 + v_c ** 3)
    I0 += w * vv ** 2 * f * dv
    I1 += w * vv ** 4 * f * dv
    I2 += w * vv ** 0 * f * dv
E_mean = 0.5 * m_alpha * I1 / I0
n_fast = (S_a / V) * t_therm
p_alpha_fast = (2.0 / 3.0) * n_fast * E_mean
beta_fast = 2.0 * mu0 * p_alpha_fast / (B * B)
p_total = 2.0 * n * T
beta_tot = 2.0 * mu0 * p_total / (B * B)
p('慢化分布 f(v) ~ 1/(v^3+v_c^3)，积分区间 [v_th,i, v_a0]')
p('E_mean (平均能量)     = %.6e J = %.1f keV', E_mean, E_mean / keV)
p('n_fast = (S/V)*t_therm = %.5e m^-3', n_fast)
p('n_fast/n_e           = %.6f = %.4f%%', n_fast / n, n_fast / n * 100)
p('p_fast = (2/3) n_fast E_mean = %.5e Pa', p_alpha_fast)
p('p_th   = 2 n T       = %.5e Pa', p_total)
p('p_fast/p_th          = %.6f = %.3f%%', p_alpha_fast / p_total, p_alpha_fast / p_total * 100)
p('beta_fast            = %.6f = %.4f%%', beta_fast, beta_fast * 100)
p('beta_tot             = %.6f = %.4f%%', beta_tot, beta_tot * 100)
p('beta_fast/beta_tot   = %.6f = %.3f%%', beta_fast / beta_tot, beta_fast / beta_tot * 100)
# 高能尾巴：f(v) 在 v>v_c 时 ~ 1/v^3，给出 v>2v_c 的粒子份额与能量份额
num_hi = 0.0
den_hi = 0.0
num_hiE = 0.0
den_hiE = 0.0
for k in range(NGRID + 1):
    vv = v1 + k * dv
    w = 0.5 if (k == 0 or k == NGRID) else 1.0
    f = 1.0 / (vv ** 3 + v_c ** 3)
    dn = w * vv ** 2 * f * dv
    den_hi += dn
    den_hiE += dn * 0.5 * m_alpha * vv ** 2
    if vv > 2.0 * v_c:
        num_hi += dn
        num_hiE += dn * 0.5 * m_alpha * vv ** 2
p('v > 2 v_c 的粒子份额  = %.6f = %.3f%%', num_hi / den_hi, num_hi / den_hi * 100)
p('v > 2 v_c 的能量份额  = %.6f = %.3f%%', num_hiE / den_hiE, num_hiE / den_hiE * 100)

# ================================================================
head('5. 快粒子与 Alfven 波：特征频率（只给共振条件，不编阈值）')
# ================================================================
rho_m = n_i * m_i
v_A = B / math.sqrt(mu0 * rho_m)
p('rho_m = n_i m_i      = %.6e kg/m^3', rho_m)
p('v_A                  = %.6e m/s', v_A)
p('v_a0/v_A             = %.6f  (>1 -> 与 TAE 共振)', v_a0 / v_A)
for q in (1.0, 2.0, 3.0):
    p('f_TAE = v_A/(4 pi q R), q=%.0f -> %.1f kHz', q, v_A / (4 * math.pi * q * R) / 1e3)
w_d = E_a0 / (Z_f * e * B * R * R)
p('进动频率 omega_d = E/(Z e B R^2) = %.2f rad/s -> f = %.1f Hz', w_d, w_d / (2 * math.pi))
qq = 2.0
rr = 1.0
eps_r = rr / R
w_b = (v_a0 / (qq * R)) * math.sqrt(eps_r / 2.0)
p('捕获 alpha 弹跳频率(r=1m, q=2) = %.3e rad/s -> %.2f kHz', w_b, w_b / (2 * math.pi) / 1e3)
p('三个频率之比 f_d:f_b:f_TAE = 1 : %.1f : %.1f', (w_b / (2 * math.pi)) / (w_d / (2 * math.pi)),
  (v_A / (4 * math.pi * 2 * R)) / (w_d / (2 * math.pi)))
# 与热离子回旋频率、漂移频率对比
f_cD = e * B / (2 * math.pi * m_D)
p('氘离子回旋频率 f_cD   = %.4e Hz = %.2f MHz', f_cD, f_cD / 1e6)

# ================================================================
head('6. 快粒子损失：纹波、轨道、局部热斑')
# ================================================================
delta_ripple = 0.005
v_par_crit = math.sqrt(2.0 * delta_ripple)
p('纹波幅度 delta        = %.4f', delta_ripple)
p('纹波捕获条件 |v_par/v| < sqrt(2 delta) = %.4f  -> 损失锥立体角份额 ~ %.2f%%',
  v_par_crit, (1.0 - math.sqrt(1.0 - v_par_crit ** 2)) * 100)
p('简化模型取 1%%~3%% 的 alpha 落入纹波损失锥 -> %.3f~%.3f MW',
  0.01 * P_alpha / 1e6, 0.03 * P_alpha / 1e6)
# 局部热斑：首次轨道损失 1.3% 的功率落在多大面积上
A_strike = 2.0 * math.pi * R * 0.2     # 简化：沿极向 0.2 m 宽、环向一圈
p('首次轨道损失功率      = %.3f MW', P_prompt / 1e6)
p('简化落点面积 2piR*0.2m = %.4f m^2', A_strike)
p('对应热流              = %.4f MW/m^2', P_prompt / 1e6 / A_strike)

# ================================================================
head('7. 快粒子诊断的量级（NPA / FIDA 的可观测量）')
# ================================================================
p('n_fast/n_e            = %.5f%%  (中性粒子分析 NPA 的信号量级)',
  n_fast / n * 100)
p('E_mean/E_birth        = %.6f  (FIDA 谱重心相对出生能量)', E_mean / E_a0)
p('v_c/v_a0              = %.6f  (慢化谱的拐点位置)', v_c / v_a0)

# ================================================================
head('9. 印刷精度自检（正文与门禁 P27 只允许用这些值）')
# ================================================================
print('')
p('[1] E_alpha0  = 3.5 MeV = %.4e J', E_a0)
p('[2] v_alpha0  = %.4e m/s = %.4f c', v_a0, v_a0 / c)
p('[3] rho_alpha = %.4f m = %.2f cm', rho_a, rho_a * 100)
p('[4] S_alpha   = %.4e /s ; S_alpha/V = %.4e /(m^3 s)', S_a, S_a / V)
p('[5] rho_i(D)  = %.4e m = %.2f mm', math.sqrt(2.0 * m_D * T) / (e * B),
  math.sqrt(2.0 * m_D * T) / (e * B) * 1e3)
p('[6] v_Te      = %.4e m/s ; x_e = %.4f ; H/H_asym = %.4f', v_Te, x_e, H_exact / H_asym)
p('[7] K         = %.4e kg/s ; tau_v = %.4f s ; tau_E = %.4f s', K, tau_v, tau_E)
p('[8] tau_eq    = %.4f s (交叉验证, tau_eq/tau_E = %.3f)', tau_eq_ref, tau_eq_ref / tau_E)
p('[9] v_c       = %.4e m/s ; E_c = %.4e J = %.1f keV = %.2f T_e', v_c, E_c, E_c / keV, E_c / T)
p('[10] t_therm  = %.4f s ; 路径 = %.4e m = %.0f km ; 绕行 %.4e 次',
  t_therm, v_a0 * tau_v, v_a0 * tau_v / 1e3, v_a0 * tau_v / (2 * math.pi * R))
p('[11] rho_th    = %.4f m = %.2f cm = %.2f%% a ; Delta_a = 2 rho_th = %.4f m',
  rho_theta, rho_theta * 100, rho_theta / a * 100, Delta_a)
p('[12] 首次轨道损失 = %.4f%% = %.4f MW ; 边界 r/a = %.4f',
  f_prompt * 100, P_prompt / 1e6, x_edge)
p('[13] E_mean   = %.4e J = %.1f keV ; n_fast = %.4e m^-3 ; n_fast/n_e = %.4f%%',
  E_mean, E_mean / keV, n_fast, n_fast / n * 100)
p('[14] p_fast   = %.4e Pa ; p_th = %.4e Pa ; 比 = %.4f%%', p_alpha_fast, p_total,
  p_alpha_fast / p_total * 100)
p('[15] beta_fast= %.6f%% ; beta_tot = %.4f%% ; 比 = %.4f%%',
  beta_fast * 100, beta_tot * 100, beta_fast / beta_tot * 100)
p('[16] 尾巴: v>2v_c 粒子 %.4f%% 、能量 %.4f%%', num_hi / den_hi * 100, num_hiE / den_hiE * 100)
p('[17] v_A      = %.4e m/s ; v_a0/v_A = %.4f', v_A, v_a0 / v_A)
p('[18] f_TAE(q=2) = %.2f kHz ; f_d = %.2f Hz ; f_b = %.2f kHz',
  v_A / (4 * math.pi * 2 * R) / 1e3, w_d / (2 * math.pi), w_b / (2 * math.pi) / 1e3)
p('[19] f_cD     = %.4e Hz = %.2f MHz', f_cD, f_cD / 1e6)
p('[20] 纹波: sqrt(2 delta) = %.4f ; 取 1%%~3%% -> %.4f~%.4f MW',
  v_par_crit, 0.01 * P_alpha / 1e6, 0.03 * P_alpha / 1e6)
p('[21] 热斑: %.4f MW / %.4f m^2 = %.4f MW/m^2', P_prompt / 1e6, A_strike,
  P_prompt / 1e6 / A_strike)
p('[22] E_mean/E_birth = %.4f ; v_c/v_a0 = %.4f', E_mean / E_a0, v_c / v_a0)
p('')
p('本文件是附录 W 正文与门禁 P27 的唯一数值来源。')

f = 'tmp/appW_calc.txt'
import io
io.open(f, 'w', encoding='utf-8', newline='\n').write('\n'.join(out) + '\n')
print('')
print('已写出 %s（%d 行）' % (f, len(out) + 1))
