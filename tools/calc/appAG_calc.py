# -*- coding: utf-8 -*-
"""附录 AG（磁流体平衡与位形）的唯一数值来源。只打印，不写文件。

用法：python tmp\\appAG_calc.py | Out-File -Encoding utf8 tmp\\appAG_calc.txt

工作点沿用 D1 家族（与附录 M/S/U/V/W/X/Z/AD/AE/AF 同机）：
  R=6.2 m、a=2.0 m、B0=5.3 T、I_p=15 MA、kappa=1.7、delta=0.4、
  n_e=1e20 m^-3、T_e=T_i=9194 eV、m_i=2.5 m_p。
所有正文与门禁可引用的数字都从这里打印。
"""
import math

# ---------------------------------------------------------------- 常数
mu0 = 1.25663706212e-6
e = 1.602176634e-19
mp = 1.67262192369e-27
keV = 1.602176634e-16

# ---------------------------------------------------------------- [DECLARED] 输入
R = 6.2
a = 2.0
B0 = 5.3
Ip = 15.0e6
kappa = 1.7
delta = 0.4
n_e = 1.0e20
Te_eV = 9194.0
mi = 2.5 * mp
li = 1.0
q0 = 1.0
q95_target = 3.0
betaN_limit = 3.0
eps_alpha = 0.05

print('=== 附录 AG 磁流体平衡与位形：数值来源（唯一） ===')
print('[DECLARED] R=%.4g m  a=%.4g m  B0=%.4g T  I_p=%.6g MA  kappa=%.4g  delta=%.4g'
      % (R, a, B0, Ip / 1e6, kappa, delta))
print('[DECLARED] n_e=%.6g m^-3  T_e=T_i=%.6g eV  m_i=%.4g m_p  l_i=%.4g  q0=%.4g'
      % (n_e, Te_eV, mi / mp, li, q0))
print('[DECLARED] q95 目标=%.4g  beta_N 上限=%.4g  alpha 参考=%.4g' % (q95_target, betaN_limit, eps_alpha))
print('')

# ---------------------------------------------------------------- [AG.1] 磁面与安全因子
eps = a / R
Aspect = R / a
Btheta = mu0 * Ip / (2.0 * math.pi * a)
q_cyl = a * B0 / (R * Btheta)
q95 = q_cyl * (1.0 + kappa * kappa) / 2.0
shear95 = 2.0 * (q95 - q0) / q95
Phi_tor = B0 * math.pi * a * a * kappa
Phi_pol = Phi_tor / q95
Bpol_energy = Btheta * Btheta / (2.0 * mu0)
Btor_energy = B0 * B0 / (2.0 * mu0)
V_plasma = 2.0 * math.pi * math.pi * R * a * a * kappa
S_wall = 4.0 * math.pi * math.pi * R * a

print('[AG.1] 磁面、安全因子与位形几何')
print('  epsilon=a/R=%.6g' % eps)
print('  A=R/a=%.6g' % Aspect)
print('  B_theta(a)=mu0 I_p/(2 pi a)=%.6g T' % Btheta)
print('  B_theta/B0=%.6g = %.6g %%' % (Btheta / B0, Btheta / B0 * 100.0))
print('  (B0/B_theta)^2=%.6g（环向与极向比压之比）' % ((B0 / Btheta) ** 2))
print('  q_cyl=a B0/(R B_theta)=%.6g' % q_cyl)
print('  q95=q_cyl (1+kappa^2)/2=%.6g' % q95)
print('  磁剪切 s95=2(q95-q0)/q95=%.6g' % shear95)
print('  环向磁通 Phi_tor=B0 pi a^2 kappa=%.6g Wb' % Phi_tor)
print('  极向磁通 Phi_pol=Phi_tor/q95=%.6g Wb' % Phi_pol)
print('  极向场能量密度 B_theta^2/(2mu0)=%.6g J/m^3' % Bpol_energy)
print('  环向场能量密度 B0^2/(2mu0)=%.6g J/m^3' % Btor_energy)
print('  等离子体体积 V=2 pi^2 R a^2 kappa=%.6g m^3' % V_plasma)
print('  壁面积（圆柱近似）S=4 pi^2 R a=%.6g m^2' % S_wall)
print('  环向场沿大半径衰减：内中平面 B(R-a)=%.6g T、外中平面 B(R+a)=%.6g T'
      % (B0 * R / (R - a), B0 * R / (R + a)))
print('')

# ---------------------------------------------------------------- [AG.2] Grad-Shafranov
p_avg = 2.0 * n_e * Te_eV * e
beta_t = 2.0 * mu0 * p_avg / (B0 * B0)
beta_p = 2.0 * mu0 * p_avg / (Btheta * Btheta)
beta_N = beta_t * 100.0 * a * B0 / (Ip / 1e6)
F_profile = B0 * R
print('[AG.2] Grad--Shafranov 方程的三个系数')
print('  <p>=2 n_e T_e=%.6g Pa' % p_avg)
print('  beta_t=2 mu0 <p>/B0^2=%.6g = %.6g %%' % (beta_t, beta_t * 100.0))
print('  beta_p=2 mu0 <p>/B_theta^2=%.6g = %.6g %%' % (beta_p, beta_p * 100.0))
print('  beta_N=beta_t[%%] a B0/I_p[MA]=%.6g' % beta_N)
print('  压强梯度量级 |dp/dr|~<p>/a=%.6g Pa/m' % (p_avg / a))
print('  电流密度量级 <j>=I_p/(pi a^2 kappa)=%.6g A/m^2' % (Ip / (math.pi * a * a * kappa)))
print('  环向场函数 F=B0 R=%.6g T m' % F_profile)
print('  压强与磁压之比（轴心量级）<p>/(B0^2/2mu0)=%.6g' % (p_avg / (B0 * B0 / (2.0 * mu0))))
print('')

# ---------------------------------------------------------------- [AG.3] Shafranov 位移
Delta = (a * a / (2.0 * R)) * (beta_p + li / 2.0)
R_axis = R + Delta
R_in = R - a
R_out = R + a
gap_in = a - Delta
gap_out = a + Delta
ratio_gap = gap_out / gap_in
Delta_limit = a / 2.0
bp_limit = Delta_limit * 2.0 * R / (a * a) - li / 2.0
print('[AG.3] Shafranov 位移')
print('  Delta=(a^2/2R)(beta_p+l_i/2)=%.6g m' % Delta)
print('  Delta/a=%.6g = %.6g %%' % (Delta / a, Delta / a * 100.0))
print('  磁轴位置 R_axis=R+Delta=%.6g m' % R_axis)
print('  内外侧磁面间距：内侧 %.6g m、外侧 %.6g m，比值 %.6g' % (gap_in, gap_out, ratio_gap))
print('  轴心环向场 B0 R/R_axis=%.6g T' % (B0 * R / R_axis))
print('  Delta=a/2 需要 beta_p=%.6g（远超当前 %.6g）' % (bp_limit, beta_p))
print('  位移随 beta_p 线性：beta_p 翻倍时 Delta=%.6g m' % ((a * a / (2.0 * R)) * (2.0 * beta_p + li / 2.0)))
print('')

# ---------------------------------------------------------------- [AG.4] 位形参数
shape_1 = 1.0
shape_13 = (1.0 + 1.3 * 1.3) / 2.0
shape_17 = (1.0 + kappa * kappa) / 2.0
shape_20 = (1.0 + 2.0 * 2.0) / 2.0
Ip_circ = Ip / shape_17
Ip_20 = Ip * shape_20 / shape_17
delta_beta_gain = 1.2
print('[AG.4] 伸长比与三角形变买到了什么')
print('  形状因子 (1+kappa^2)/2：kappa=1 -> %.6g、kappa=1.3 -> %.6g、kappa=1.7 -> %.6g、kappa=2 -> %.6g'
      % (shape_1, shape_13, shape_17, shape_20))
print('  同样 q95 下环向电流上限之比（相对圆截面）=%.6g' % shape_17)
print('  同样 q95 下允许电流：kappa=1.3 -> %.6g MA、kappa=1.7 -> %.6g MA、kappa=2.0 -> %.6g MA'
      % (Ip * shape_13 / shape_17 / 1e6, Ip / 1e6, Ip * shape_20 / shape_17 / 1e6))
print('  圆截面要达到同样 q95 只能通 %.6g MA（现在是 %.6g MA）' % (Ip_circ / 1e6, Ip / 1e6))
print('  若 kappa 提到 2.0，同样 q95 可通 %.6g MA' % (Ip_20 / 1e6))
print('  三角形变增益（声明）delta 从 0 到 0.4：beta_N 上限 +%.4g %%' % ((delta_beta_gain - 1.0) * 100.0))
print('  同样 kappa 下体积随 a^2 R 线性：a 加 10%%、R 加 5%% -> V=%.6g m^3'
      % (2.0 * math.pi * math.pi * (R * 1.05) * (a * 1.1) ** 2 * kappa))
print('')

# ---------------------------------------------------------------- [AG.5] 压强极限
beta_t_troyon = betaN_limit * (Ip / 1e6) / (a * B0) / 100.0
pressure_troyon = beta_t_troyon * B0 * B0 / (2.0 * mu0)
frac_troyon = beta_t / beta_t_troyon
Ip_q3 = 2.0 * math.pi * a * a * B0 * shape_17 / (mu0 * R * q95_target)
nG = Ip_q3 / 1e6 / (math.pi * a * a) * 1e20
frac_G = n_e / nG
virial = eps / (q95 * q95)
print('[AG.5] 压强极限：beta 能推到多高')
print('  beta_N 上限 %.4g 对应 beta_t=%.6g = %.6g %%' % (betaN_limit, beta_t_troyon, beta_t_troyon * 100.0))
print('  对应 <p> 上限=%.6g Pa' % pressure_troyon)
print('  当前工作点用了极限的 %.6g %%' % (frac_troyon * 100.0))
print('  维里量级 eps/q95^2=%.6g = %.6g %%' % (virial, virial * 100.0))
print('  q95=%.4g 需要的电流 I_p=%.6g MA（当前 %.4g MA）' % (q95_target, Ip_q3 / 1e6, Ip / 1e6))
print('  该电流下 Greenwald 密度 n_G=%.6g m^-3，工作点 n_e/n_G=%.6g' % (nG, frac_G))
print('  压强梯度驱动的气球模参数 alpha=2 mu0 q95^2 R <p>/(a B0^2)=%.6g'
      % (2.0 * mu0 * q95 * q95 * R * p_avg / (a * B0 * B0)))
print('  与附录 AD 的 alpha_crit=2.5 相比：比值 %.6g'
      % (2.0 * mu0 * q95 * q95 * R * p_avg / (a * B0 * B0) / 2.5))
print('')

# ---------------------------------------------------------------- [AG.6] 从平衡到稳定性
Btheta_95 = Btheta * (a / (0.95 * a))
q_at_95 = q95
shear_at_95 = shear95
mercier_margin = 0.5 * shear95 * shear95
print('[AG.6] 从平衡到稳定性（只给量级，不给稳定性证明）')
print('  q95=%.6g、磁剪切 s95=%.6g' % (q_at_95, shear_at_95))
print('  Suydam 型量级 s95^2/8=%.6g' % (shear95 * shear95 / 8.0))
print('  Mercier 型量级 s95^2/2=%.6g' % mercier_margin)
print('  内电感 l_i=%.4g、边界 q 与轴心 q 之比 q95/q0=%.6g' % (li, q95 / q0))
print('  极向 beta 与内电感之和 beta_p+l_i/2=%.6g（决定位移）' % (beta_p + li / 2.0))
print('  与附录 AD 衔接：台基处 alpha/alpha_crit=1.31084（无自举）/0.819273（含自举）')
print('')

# ---------------------------------------------------------------- [AG.7] 练习解答
print('[AG.7] 六道练习的解答值')
n2, T2, a2, R2, B2, I2, k2 = 0.8e20, 8000.0, 1.5, 5.5, 5.3, 8.0e6, 1.6
e2 = a2 / R2
bt2 = mu0 * I2 / (2.0 * math.pi * a2)
q2 = a2 * B2 / (R2 * bt2) * (1.0 + k2 * k2) / 2.0
p2 = 2.0 * n2 * T2 * e
beta_t2 = 2.0 * mu0 * p2 / (B2 * B2)
beta_p2 = 2.0 * mu0 * p2 / (bt2 * bt2)
d2 = (a2 * a2 / (2.0 * R2)) * (beta_p2 + li / 2.0)
print('  1) 练习1：epsilon=%.6g、B_theta=%.6g T、q95=%.6g、beta_t=%.6g %%'
      % (e2, bt2, q2, beta_t2 * 100.0))
print('     环径比 R/a=%.6g' % (R2 / a2))
print('     beta_p=%.6g、Delta=%.6g m、Delta/a=%.6g %%' % (beta_p2, d2, d2 / a2 * 100.0))
print('  2) 同样 q95、kappa 1.0 -> 2.0 允许电流倍数 = %.6g' % (shape_20 / shape_1))
print('  3) Delta 使外侧间距变 %.6g 倍、内侧变 %.6g 倍' % ((a + Delta) / a, (a - Delta) / a))
V3 = 2.0 * math.pi * math.pi * 6.5 * 2.2 ** 2 * 1.8
V_ref = 832.205
print('  4) a=2.2、R=6.5、kappa=1.8 的体积 = %.6g m^3（是 %.6g m^3 的 %.6g 倍）'
      % (V3, V_ref, V3 / V_ref))
bt5 = mu0 * Ip / (2.0 * math.pi * a)
p5 = betaN_limit * (Ip / 1e6) / (a * B0) / 100.0 * B0 * B0 / (2.0 * mu0)
print('  5) beta_N=%.4g 对应 <p>=%.6g Pa（当前 %.6g Pa，比值 %.6g）'
      % (betaN_limit, p5, p_avg, p5 / p_avg))
Ip6 = 2.0 * math.pi * a * a * B0 * shape_17 / (mu0 * R * q95_target)
print('  6) q95=3.0 需要 I_p=%.6g MA（降 %.6g 倍）、n_G=%.6g m^-3、n_e/n_G=%.6g'
      % (Ip6 / 1e6, Ip / Ip6, nG, frac_G))
print('')

print('=== 可以引用的全部数值结束 ===')
