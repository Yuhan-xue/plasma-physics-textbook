# -*- coding: utf-8 -*-
"""附录 AF（等离子体湍流与反常输运）的数值唯一来源。

只打印，不写文件。用法：
    python tmp\\appAF_calc.py | Out-File -Encoding utf8 tmp\\appAF_calc.txt
所有正文与门禁里标出的数字都必须能在这里找到（或写进 [DECLARED] 段说明是输入/文献值）。
工作点与附录 M/S/U/V/W/X/Y/Z/AA/AB/AC/AD/AE 一致（D1）。
"""
import math

# ---------------------------------------------------------------- 常数（CODATA 2018）
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
eps0 = 8.8541878128e-12
keV = 1.602176634e-16
MeV = 1.602176634e-13

# ---------------------------------------------------------------- 声明输入（DECLARED）
n_e = 1.0e20          # m^-3
T_e_eV = 9.194e3      # eV
T_e = T_e_eV * e      # J
B = 5.3               # T
R = 6.2               # m
a = 2.0               # m
m_i = 2.5 * mp        # kg
Z_eff = 2.0
tau_E = 3.7           # s   （与附录 AD 同一工作点）
lnL = 17.0
S_wall = 489.532      # m^2 （与附录 X 同值）
W_th = 371.4e6        # J
P_alpha = 111.3e6     # W
L_n = 1.0             # m   密度梯度长度（声明）
L_T = 0.5             # m   温度梯度长度（声明）
kth_rho_s = 0.3       # 湍流谱峰无量纲波数（声明）
RLT_crit = 5.0        # 临界 R/L_T（文献量级，声明）
eta_i_crit = 1.5      # ITG 阈值（文献量级，声明）
zeta_mix = 1.0        # 混合长系数（声明；χ = ζ γ/k⊥²）
delta_n_rel = 0.01    # 涨落水平（声明）
E_r_zf = 53.2405e3    # V/m 带状流径向电场（取自附录 AD 的剪切口径）

print('=' * 78)
print('附录 AF 数值清单（等离子体湍流与反常输运）——唯一来源 tmp/appAF_calc.py')
print('=' * 78)
print('[DECLARED] n_e=%.6g m^-3  T_e=T_i=%.6g eV  B=%.6g T  R=%.6g m  a=%.6g m'
      % (n_e, T_e_eV, B, R, a))
print('[DECLARED] m_i=2.5 m_p  Z_eff=%.6g  tau_E=%.6g s  lnLambda=%.6g'
      % (Z_eff, tau_E, lnL))
print('[DECLARED] L_n=%.6g m  L_T=%.6g m  k_theta rho_s=%.6g  R/L_T 临界=%.6g'
      % (L_n, L_T, kth_rho_s, RLT_crit))
print('[DECLARED] eta_i 阈值=%.6g  混合长系数=%.6g  涨落水平=%.6g'
      % (eta_i_crit, zeta_mix, delta_n_rel))
print('[DECLARED] S=%.6g m^2  W_th=%.6g MJ  P_alpha=%.6g MW'
      % (S_wall, W_th / 1e6, P_alpha / 1e6))

# ---------------------------------------------------------------- AF.1 三个输运系数
omega_ci = e * B / m_i
c_s = math.sqrt(T_e / m_i)
rho_s = c_s / omega_ci
rho_i = math.sqrt(2.0 * T_e / m_i) / omega_ci
chi_gB_i = rho_i ** 2 * c_s / a
chi_gB_s = rho_s ** 2 * c_s / a
chi_B = T_e / (16.0 * e * B)
chi_eff = a ** 2 / (2.0 * tau_E)

print('')
print('[AF.1] 三个输运系数')
print('  离子声速度 c_s = %.6g m/s = %.6g km/s' % (c_s, c_s / 1e3))
print('  离子回旋频率 omega_ci = %.6g rad/s' % omega_ci)
print('  离子回旋半径 rho_i = %.6g m = %.6g mm' % (rho_i, rho_i * 1e3))
print('  离子声回旋半径 rho_s = %.6g m = %.6g mm' % (rho_s, rho_s * 1e3))
print('  rho_s/rho_i = %.6g   rho_i/rho_s = %.6g' % (rho_s / rho_i, rho_i / rho_s))
print('  回旋玻姆扩散（用 rho_i）= %.6g m^2/s' % chi_gB_i)
print('  回旋玻姆扩散（用 rho_s）= %.6g m^2/s' % chi_gB_s)
print('  玻姆扩散 chi_B = %.6g m^2/s' % chi_B)
print('  实验反推 chi_eff = a^2/(2 tau_E) = %.6g m^2/s' % chi_eff)
print('  chi_B/chi_gB(rho_i) = %.6g' % (chi_B / chi_gB_i))
print('  chi_B/chi_gB(rho_s) = %.6g' % (chi_B / chi_gB_s))
print('  chi_gB(rho_i)/chi_eff = %.6g（比回旋玻姆好这么多倍）' % (chi_gB_i / chi_eff))
print('  chi_gB(rho_s)/chi_eff = %.6g' % (chi_gB_s / chi_eff))
print('  chi_B/chi_eff = %.6g' % (chi_B / chi_eff))
print('  总损失功率 W_th/tau_E = %.6g MW' % (W_th / tau_E / 1e6))
print('  湍流热流密度 n T chi_eff/L_T = %.6g MW/m^2'
      % (n_e * T_e * chi_eff / L_T / 1e6))
print('  乘迎风面积 = %.6g MW' % (n_e * T_e * chi_eff / L_T * S_wall / 1e6))
print('  占总损失功率 %.6g' % (n_e * T_e * chi_eff / L_T * S_wall / (W_th / tau_E)))
print('  占 alpha 功率 %.6g' % (n_e * T_e * chi_eff / L_T * S_wall / P_alpha))

# ---------------------------------------------------------------- AF.2 漂移波
k_th = kth_rho_s / rho_s
lam_th = 2.0 * math.pi / k_th
omega_star = kth_rho_s * c_s / L_n
f_star = omega_star / (2.0 * math.pi)
v_dia = rho_s * c_s / L_n
gamma_max = kth_rho_s * c_s / L_n

print('')
print('[AF.2] 漂移波：梯度与磁场谁说了算')
print('  取 k_theta rho_s = %.6g 时的垂直波数 k_theta = %.6g m^-1' % (kth_rho_s, k_th))
print('  对应垂直波长 lambda_theta = %.6g m = %.6g cm' % (lam_th, lam_th * 100))
print('  抗磁漂移速度 v_dia = rho_s c_s/L_n = %.6g m/s = %.6g km/s'
      % (v_dia, v_dia / 1e3))
print('  漂移波频率 omega_star = k_theta rho_s c_s/L_n = %.6g rad/s' % omega_star)
print('  对应频率 f_star = %.6g Hz = %.6g kHz' % (f_star, f_star / 1e3))
print('  omega_star/omega_ci = %.6g' % (omega_star / omega_ci))
print('  线性增长率量级 gamma_max = %.6g s^-1' % gamma_max)
print('  特征时间 1/gamma_max = %.6g s = %.6g us'
      % (1.0 / gamma_max, 1e6 / gamma_max))
print('  一周环向距离 2 pi R = %.6g m' % (2 * math.pi * R))
print('  湍流涡在环向的个数 2 pi R/lambda_theta = %.6g'
      % (2 * math.pi * R / lam_th))

# ---------------------------------------------------------------- AF.3 ITG 与 TEM
eta_i = L_n / L_T
RLT = R / L_T
RLn = R / L_n
LT_crit = R / RLT_crit
excess = (RLT - RLT_crit) / RLT_crit
f_trap = 0.69843   # 与附录 AE 的捕获比例同口径（声明）

print('')
print('[AF.3] ITG 与 TEM：阈值在哪里')
print('  eta_i = L_n/L_T = %.6g' % eta_i)
print('  R/L_T = %.6g   R/L_n = %.6g' % (RLT, RLn))
print('  绝对值 |grad T|/T = 1/L_T = %.6g m^-1' % (1.0 / L_T))
print('  临界温度梯度长度 L_T,crit = R/(R/L_T)_crit = %.6g m' % LT_crit)
print('  实际梯度是临界值的 %.6g 倍' % (LT_crit / L_T))
print('  超临界幅度 (R/L_T - 临界)/临界 = %.6g' % excess)
print('  ITG 阈值 eta_i,crit = %.6g（声明）' % eta_i_crit)
print('  实际 eta_i 是阈值的 %.6g 倍' % (eta_i / eta_i_crit))
print('  捕获电子比例（附录 AE 口径）= %.6g' % f_trap)
print('  TEM 驱动权重（按捕获比例估）= %.6g' % f_trap)

# ---------------------------------------------------------------- AF.4 混合长与饱和
tau_c = 1.0 / gamma_max
dr_mix = 1.0 / k_th
chi_mix = zeta_mix * gamma_max / k_th ** 2
dr_eff = math.sqrt(chi_eff * tau_c)
print('')
print('[AF.4] 混合长与饱和：上界与真值')
print('  相关时间 tau_c = 1/gamma_max = %.6g s = %.6g us'
      % (tau_c, tau_c * 1e6))
print('  混合长步长 1/k_theta = %.6g m = %.6g mm' % (dr_mix, dr_mix * 1e3))
print('  混合长扩散 chi_mix = gamma_max/k_theta^2 = %.6g m^2/s' % chi_mix)
print('  chi_mix/chi_eff = %.6g（混合长高估的倍数）' % (chi_mix / chi_eff))
print('  由 chi_eff 反推的有效步长 dr = sqrt(chi_eff tau_c) = %.6g m = %.6g mm'
      % (dr_eff, dr_eff * 1e3))
print('  dr/rho_s = %.6g   dr/(1/k_theta) = %.6g' % (dr_eff / rho_s, dr_eff * k_th))
print('  dr_mix/rho_s = %.6g' % (dr_mix / rho_s))
print('  涨落水平 delta_n/n = %.6g（声明）' % delta_n_rel)
print('  对应的扩散系数量级 (delta_n/n)^2 dr_mix^2/tau_c = %.6g m^2/s'
      % (delta_n_rel ** 2 * dr_mix ** 2 / tau_c))

# ---------------------------------------------------------------- AF.5 带状流与剪切
v_zf = E_r_zf / B
dr_zf_rho = rho_i
gamma_E_rho = v_zf / dr_zf_rho
dr_zf_macro = a / 10.0
gamma_E_macro = v_zf / dr_zf_macro
print('')
print('[AF.5] 带状流与 E x B 剪切')
print('  带状流速度 v_zf = E_r/B = %.6g m/s = %.6g km/s' % (v_zf, v_zf / 1e3))
print('  径向尺度取 rho_i = %.6g mm 时剪切率 gamma_E = %.6g s^-1'
      % (dr_zf_rho * 1e3, gamma_E_rho))
print('  gamma_E/gamma_max = %.6g（大于 1 才能抑制）' % (gamma_E_rho / gamma_max))
print('  径向尺度取 a/10 = %.6g m 时 gamma_E = %.6g s^-1' % (dr_zf_macro, gamma_E_macro))
print('  gamma_E/gamma_max = %.6g（小于 1，抑制不了）' % (gamma_E_macro / gamma_max))
print('  剪切率与增长率相等时的临界尺度 = %.6g m = %.6g mm'
      % (v_zf / gamma_max, v_zf / gamma_max * 1e3))

# ---------------------------------------------------------------- AF.6 谱与级联
print('')
print('[AF.6] 谱与级联：输运由哪一段波数贡献')
k1, k2, k3 = 0.1, 0.5, 1.0
frac_band = math.log(k2 / k1) / math.log(k3 / k1)
print('  谱峰 k_theta rho_s = %.6g' % kth_rho_s)
print('  二维漂移波湍流：熵级串谱指数 -3、逆能量级串 -5/3（文献，声明）')
print('  三维 Kolmogorov 谱指数 -5/3（文献，声明）')
print('  若 S(k) ~ k^-3，则 k^2 S(k) ~ 1/k')
print('  0.1—0.5 段贡献占 0.1—1.0 段的比例 = %.6g' % frac_band)
print('  归一化到全谱（0.1—1.0）的剩余比例 = %.6g' % (1.0 - frac_band))
print('  波数带宽比 k2/k1 = %.6g  k3/k1 = %.6g' % (k2 / k1, k3 / k1))

# ---------------------------------------------------------------- AF.7 练习解答
print('')
print('[AF.7 练习解答]')
# 1) B = 10 T
B2 = 10.0
w_ci2 = e * B2 / m_i
rho_s2 = c_s / w_ci2
print('  1) B=%.6g T：omega_ci=%.6g rad/s  rho_s=%.6g m=%.6g mm  k_theta=%.6g m^-1'
      % (B2, w_ci2, rho_s2, rho_s2 * 1e3, kth_rho_s / rho_s2))
print('     chi_B=%.6g m^2/s  chi_gB(rho_s)=%.6g m^2/s  比值=%.6g'
      % (T_e / (16.0 * e * B2), rho_s2 ** 2 * c_s / a,
         (T_e / (16.0 * e * B2)) / (rho_s2 ** 2 * c_s / a)))
print('     chi_gB 随 B 缩小 %.6g 倍  chi_B 缩小 %.6g 倍'
      % (chi_gB_s / (rho_s2 ** 2 * c_s / a), chi_B / (T_e / (16.0 * e * B2))))
print('     omega_star 不变（固定 k_theta rho_s）= %.6g rad/s' % omega_star)
# 2) L_T = 0.3 m
LT2 = 0.3
print('  2) L_T=%.6g m：R/L_T=%.6g  超临界 %.6g 倍  L_T,crit=%.6g m  实际/临界=%.6g'
      % (LT2, R / LT2, (R / LT2 - RLT_crit) / RLT_crit, LT_crit, LT_crit / LT2))
# 3) 涨落与步长
print('  3) 若 delta_n/n=%.6g：chi_mix=%.6g m^2/s  有效步长 dr=%.6g mm  dr/rho_s=%.6g'
      % (delta_n_rel, chi_mix, dr_eff * 1e3, dr_eff / rho_s))
# 4) 定标
print('  4) B 翻倍：chi_gB 变 %.6g 倍、chi_B 变 %.6g 倍，两者之比变 %.6g 倍'
      % (0.25, 0.5, 2.0))
# 5) 要达到 chi_eff = 1 m^2/s
chi_t = 1.0
print('  5) 若要 chi_eff=%.6g m^2/s：tau_E=a^2/(2 chi_eff)=%.6g s'
      % (chi_t, a ** 2 / (2 * chi_t)))
print('     是现在 tau_E=%.6g s 的 %.6g 倍（更短，约束要变差）'
      % (tau_E, (a ** 2 / (2 * chi_t)) / tau_E))
# 6) 带状流抑制
print('  6) 带状流：v_zf=%.6g m/s  gamma_E(rho_i)=%.6g s^-1 = %.6g gamma_max'
      % (v_zf, gamma_E_rho, gamma_E_rho / gamma_max))
print('     抑制失效的临界尺度 = %.6g mm（比 rho_i=%.6g mm 大 %.6g 倍）'
      % (v_zf / gamma_max * 1e3, rho_i * 1e3, (v_zf / gamma_max) / rho_i))

print('')
print('=' * 78)
print('以上即 tmp/appAF_calc.py 可以引用的全部数值。')
print('=' * 78)
