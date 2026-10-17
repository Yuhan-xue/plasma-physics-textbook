# -*- coding: utf-8 -*-
"""附录 AH 数值唯一来源：原子过程与谱线（只 print，不写文件）。

用法：python tmp/appAH_calc.py | Out-File -Encoding utf8 tmp/appAH_calc.txt
"""
import math

# ---- CODATA 2018 ----
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
mu_atom = 1.66053906660e-27
h = 6.62607015e-34
hbar = 1.054571817e-34
kB = 1.380649e-23
eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6
alpha = 7.2973525693e-3
a0 = 5.29177210903e-11
Ry_eV = 13.605693122994
mec2_eV = 510998.95
muc2_eV = 931.49410242e6
mu_B = 9.2740100783e-24
sigma_T = 6.6524587321e-29
hc_eVnm = 1239.841984

print('=' * 78)
print('附录 AH 数值唯一来源：原子过程与谱线')
print('=' * 78)
print('[DECLARED] 输入与声明')
ne = 1.0e20           # m^-3 主工作点（与其它附录同口径）
Te = 9194.0           # eV
Ti = 9194.0           # eV
B0 = 5.3              # T
Zeff = 2.0
ne_edge = 1.0e20      # m^-3 边缘/偏滤器
Te_sol = 100.0        # eV
ne_low = 1.0e16       # m^-3 低温等离子体（附录 P 口径）
lam_Da = 656.28       # nm D_alpha（声明）
lam_HeII = 468.57     # nm He II n=4->3（声明）
A21_Da = 4.41e7       # s^-1 D_alpha 自发辐射（声明，文献值）
dlam_inst = 0.02      # nm 仪器展宽（声明）
Stark_Da = 0.053      # nm D_alpha Stark FWHM 在 n_e=1e20 m^-3（声明，文献经验值）
chi_H = 13.605693122994
chi_He1 = 24.587      # eV
chi_He2 = 54.418
chi_C1 = 11.260
chi_C2 = 24.383
mD_u = 2.014
mHe_u = 4.0026
mC_u = 12.011
mW_u = 183.84
print('  主工作点 n_e=%.6g m^-3, T_e=T_i=%.6g eV, B=%.6g T, Z_eff=%.6g' % (ne, Te, B0, Zeff))
print('  边缘点 n_e=%.6g m^-3, T_e=%.6g eV；低温点 n_e=%.6g m^-3' % (ne_edge, Te_sol, ne_low))
print('  谱线（声明）：D_alpha %.6g nm, He II %.6g nm；A21(D_alpha)=%.6g s^-1' % (lam_Da, lam_HeII, A21_Da))
print('  声明系数：仪器展宽 %.6g nm；D_alpha Stark FWHM %.6g nm @ n_e=1e20 m^-3' % (dlam_inst, Stark_Da))
print('  电离能（eV）：H %.6g, He %.6g/%.6g, C %.6g/%.6g' % (chi_H, chi_He1, chi_He2, chi_C1, chi_C2))

# ---- [AH.1] 能级与跃迁 ----
print('')
print('[AH.1] 能级与跃迁：谱线为什么是离散的')
Z_C = 6.0
E1_C = Ry_eV * Z_C ** 2
print('  类氢 C VI 基态能量       = %.6g eV' % E1_C)
dE_Lya = Ry_eV * Z_C ** 2 * (1.0 - 1.0 / 4.0)
dE_Ha = Ry_eV * Z_C ** 2 * (1.0 / 4.0 - 1.0 / 9.0)
dE_Lyb = Ry_eV * Z_C ** 2 * (1.0 - 1.0 / 9.0)
print('  C VI Ly-alpha 能量       = %.6g eV  -> lambda = %.6g nm' % (dE_Lya, hc_eVnm / dE_Lya))
print('  C VI Balmer-alpha 能量   = %.6g eV  -> lambda = %.6g nm' % (dE_Ha, hc_eVnm / dE_Ha))
print('  C VI Ly-beta 能量        = %.6g eV  -> lambda = %.6g nm' % (dE_Lyb, hc_eVnm / dE_Lyb))
print('  D_alpha 能量（由声明波长）= %.6g eV' % (hc_eVnm / lam_Da))
v1_C = Z_C * alpha * c
print('  C VI 基态电子速度        = %.6g m/s = %.6g c' % (v1_C, v1_C / c))
print('  相对论修正量级 (v/c)^2/2 = %.6g' % (0.5 * (v1_C / c) ** 2))
print('  C VI Ly-alpha / T_e(边缘) 比 = %.6g （边缘等离子体里碳不是类氢离子）' % (dE_Lya / Te_sol))
print('  氢基态玻尔半径 a0        = %.6g m' % a0)

# ---- [AH.2] Saha 电离平衡 ----
print('')
print('[AH.2] Saha 方程：电离度对温度有多敏感')
Saha_const = 2.0 * (2.0 * math.pi * me * e / h ** 2) ** 1.5
print('  Saha 常数 2*(2*pi*m_e*k*T/h^2)^{3/2} = %.6g * T[eV]^{3/2} m^-3' % Saha_const)


def ion_frac(T_eV, n_tot, chi, g_ratio=1.0):
    """纯氢样近似：n_e^2/(n_tot-n_e) = S(T)*g_ratio*exp(-chi/T)。

    用数值稳定形式 x = 2n/(1+sqrt(1+4n/A))，避免 A 很大时的相消误差。
    """
    A = Saha_const * g_ratio * T_eV ** 1.5 * math.exp(-chi / T_eV)
    ne_ = 2.0 * n_tot / (1.0 + math.sqrt(1.0 + 4.0 * n_tot / A))
    return ne_, ne_ / n_tot


print('  -- 氢：温度扫描（n=1e20 m^-3，过渡就发生在 0.5 eV 附近）')
for T in (0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 1.00):
    ne_, f = ion_frac(T, ne_edge, chi_H)
    print('  D, n=1e20 m^-3, T=%.6g eV -> n_e=%.6g m^-3, 电离度=%.6g (=%.6g%%)'
          % (T, ne_, f, 100.0 * f))
print('  -- 氢：密度扫描（T=0.5 eV，密度越低越容易电离）')
for n_tot in (1.0e18, 1.0e19, 1.0e20, 1.0e21):
    ne_, f = ion_frac(0.5, n_tot, chi_H)
    print('  D, T=0.5 eV, n=%.6g m^-3 -> 电离度=%.6g (=%.6g%%)' % (n_tot, f, 100.0 * f))
print('  -- 杂质的一次电离（n=1e20 m^-3）：碳的电离能比氢低，所以在同一温度下更电离')
for T in (0.5, 1.0, 2.0, 5.0):
    _, fC = ion_frac(T, ne_edge, chi_C1)
    _, fH = ion_frac(T, ne_edge, chi_H)
    print('  T=%.6g eV -> C 一次电离度=%.6g (=%.6g%%)；H 电离度=%.6g (=%.6g%%)'
          % (T, fC, 100.0 * fC, fH, 100.0 * fH))
for T in (1.0, 2.0, 5.0):
    _, fHe = ion_frac(T, ne_edge, chi_He1)
    print('  He 一次电离, T=%.6g eV, n=1e20 m^-3 -> 电离度=%.6g (=%.6g%%)'
          % (T, fHe, 100.0 * fHe))
_, fC2 = ion_frac(Te_sol, 1.0e19, chi_C2)
print('  C 二次电离, n=1e19 m^-3, T=%.6g eV -> 电离度=%.6g (=%.6g%%)' % (Te_sol, fC2, 100.0 * fC2))
for T, lbl in ((Te_sol, 'SOL 100 eV'), (1000.0, '1 keV'), (Te, '主工作点 9.194 keV')):
    A = Saha_const * T ** 1.5 * math.exp(-chi_H / T)
    _, f = ion_frac(T, ne, chi_H)
    print('  H, n=1e20 m^-3, T=%-18s -> 电离度=%.6g；中性份额 = %.6g'
          % (lbl, f, 1.0 / (1.0 + A / ne)))

# ---- [AH.3] 谱线强度与强度比 ----
print('')
print('[AH.3] 强度比测温：Ly-beta / Ly-alpha 对温度有多敏感')
dE_gap = dE_Lyb - dE_Lya
print('  两条线的激发能差         = %.6g eV' % dE_gap)
for T in (300.0, 500.0, 1000.0):
    R = math.exp(-dE_gap / T)
    print('  T=%.6g eV -> R = I(Ly-beta)/I(Ly-alpha) ∝ exp(-dE/T) = %.6g' % (T, R))
T_ref = 500.0
dlnT = 0.2
dlnR = dE_gap / T_ref * dlnT
print('  在 T=500 eV 附近：T 变 20%% -> R 变 %.6g%%（灵敏度 dlnR/dlnT = %.6g）'
      % (100.0 * dlnR, dE_gap / T_ref))
print('  日冕极限判据（声明）：n_e < %.6g m^-3 时辐射衰变快于碰撞退激发' % 1.0e20)
print('  D_alpha 自发辐射 A21      = %.6g s^-1（声明）' % A21_Da)

# ---- [AH.4] 多普勒展宽 ----
print('')
print('[AH.4] 多普勒展宽：线宽测离子温度')
k8ln2 = 8.0 * math.log(2.0)


def doppler_nm(lam_nm, T_eV, m_u):
    mc2 = m_u * muc2_eV
    return lam_nm * math.sqrt(k8ln2 * T_eV / mc2)


for T in (50.0, 100.0, 200.0):
    print('  C VI Ly-alpha, T_i=%.6g eV -> Delta-lambda = %.6g nm'
          % (T, doppler_nm(hc_eVnm / dE_Lya, T, mC_u)))
for T in (50.0, 100.0):
    print('  D_alpha, T_i=%.6g eV       -> Delta-lambda = %.6g nm' % (T, doppler_nm(lam_Da, T, mD_u)))
for T in (100.0,):
    d_C = doppler_nm(hc_eVnm / dE_Lya, T, mC_u)
    d_D = doppler_nm(lam_Da, T, mD_u)
    print('  同 T_i=%.6g eV：D_alpha 线宽 / C VI 线宽 = %.6g（∝1/sqrt(m) 与波长的共同结果）'
          % (T, d_D / d_C))
print('  仪器展宽 %.6g nm 与 C VI 100 eV 线宽之比 = %.6g'
      % (dlam_inst, dlam_inst / doppler_nm(hc_eVnm / dE_Lya, 100.0, mC_u)))
print('  He II %.6g nm, T_i=200 eV -> Delta-lambda = %.6g nm'
      % (lam_HeII, doppler_nm(lam_HeII, 200.0, mHe_u)))

# ---- [AH.5] 塞曼与斯塔克 ----
print('')
print('[AH.5] 塞曼分裂（测磁场）与斯塔克展宽（测密度）')
dnu_Z = 13.996e9 * B0
print('  B=%.6g T -> 塞曼频率分裂 dnu = 13.996 GHz/T * B = %.6g Hz = %.6g GHz' % (B0, dnu_Z, dnu_Z / 1e9))
for lam in (lam_Da, hc_eVnm / dE_Lya):
    dlam_Z = lam * 1e-9 * (lam * 1e-9) * dnu_Z / c
    print('  lambda=%.6g nm -> dlambda_Z = %.6g nm' % (lam, dlam_Z * 1e9))
dZ_Da = lam_Da * 1e-9 * (lam_Da * 1e-9) * dnu_Z / c * 1e9
dZ_C = (hc_eVnm / dE_Lya) * 1e-9 * ((hc_eVnm / dE_Lya) * 1e-9) * dnu_Z / c * 1e9
print('  同 T_i=100 eV 多普勒线宽：D_alpha %.6g nm、C VI %.6g nm'
      % (doppler_nm(lam_Da, 100.0, mD_u), doppler_nm(hc_eVnm / dE_Lya, 100.0, mC_u)))
print('  塞曼/多普勒（D_alpha, 100 eV）= %.6g' % (dZ_Da / doppler_nm(lam_Da, 100.0, mD_u)))
print('  塞曼/多普勒（C VI, 100 eV）  = %.6g' % (dZ_C / doppler_nm(hc_eVnm / dE_Lya, 100.0, mC_u)))
for n_e in (1.0e20, 1.0e21):
    scale = (n_e / 1.0e20) ** (2.0 / 3.0)
    print('  n_e=%.6g m^-3 -> D_alpha Stark FWHM = %.6g nm（声明系数按 (n_e)^{2/3} 外推）'
          % (n_e, Stark_Da * scale))
print('  n_e=1e20 m^-3 时 Stark/Zeeman（D_alpha）= %.6g' % (Stark_Da / dZ_Da))

# ---- [AH.6] 不透明度与自吸收 ----
print('')
print('[AH.6] 从谱线到诊断：光学薄假设与自吸收')
lam0 = lam_Da * 1e-9
dnu_D = c / lam0 * (doppler_nm(lam_Da, 100.0, mD_u) / lam_Da)
sigma_line = lam0 ** 2 / (8.0 * math.pi) * A21_Da / dnu_D
print('  D_alpha 多普勒宽度（频率）= %.6g Hz' % dnu_D)
print('  线心吸收截面 sigma_12   = %.6g m^2 = %.6g * sigma_T' % (sigma_line, sigma_line / sigma_T))
n1 = 1.0e16
for L in (0.01, 1.0, 10.0):
    print('  n_1=%.6g m^-3, L=%.6g m -> tau = %.6g' % (n1, L, n1 * sigma_line * L))
print('  轫致吸收系数（附录 X 口径）= %.6g m^-1；线吸收/轫致（同上 n_1、L=1 m）= %.6g'
      % (1.46411e-6, n1 * sigma_line / 1.46411e-6))
tau_ff_1m = 1.46411e-6 * 1.0
print('  1 m 厚轫致吸收光学深度    = %.6g' % tau_ff_1m)

# ---- [AH.7] 练习解答 ----
print('')
print('[AH.7] 练习解答')
print('  1) He II %.6g nm, T_i=200 eV -> Delta-lambda = %.6g nm；与 C VI 同温线宽之比 = %.6g'
      % (lam_HeII, doppler_nm(lam_HeII, 200.0, mHe_u),
         doppler_nm(lam_HeII, 200.0, mHe_u) / doppler_nm(hc_eVnm / dE_Lya, 200.0, mC_u)))
ne_, f = ion_frac(0.7, 1.0e21, chi_H)
print('  2) D, n=1e21 m^-3, T=0.7 eV -> 电离度 = %.6g (=%.6g%%)' % (f, 100.0 * f))
dZ_C_keV = dZ_C
dD_C_keV = doppler_nm(hc_eVnm / dE_Lya, 1000.0, mC_u)
print('  3) C VI Ly-alpha, B=%.6g T -> dlambda_Z = %.6g nm；T_i=1 keV 多普勒线宽 = %.6g nm；比 = %.6g'
      % (B0, dZ_C_keV, dD_C_keV, dZ_C_keV / dD_C_keV))
print('  4) 若 R = I(Ly-beta)/I(Ly-alpha) = 0.5 -> T = dE/ln(1/R) = %.6g eV'
      % (dE_gap / math.log(1.0 / 0.5)))
tau_target = 1.0
n1_ex = 1.0e18
print('  5) 要 tau=%.6g，n_1=%.6g m^-3 时需要 L = %.6g m（用上节的 sigma_12）'
      % (tau_target, n1_ex, tau_target / (n1_ex * sigma_line)))
print('  6) W(%.6g u) 谱线 lambda=5 nm, T_i=%.6g eV -> Delta-lambda = %.6g nm'
      % (mW_u, Ti, doppler_nm(5.0, Ti, mW_u)))
print('')
# ---- [AH.8] 派生核对：正文里用到的比值与换算 ----
print('')
print('[AH.8] 派生核对：正文里用到的比值与换算')


def ratio_int(T):
    return math.exp(-dE_gap / T)


print('  R(300)/R(500) = %.6g；R(1000)/R(500) = %.6g；R(1000)/R(300) = %.6g'
      % (ratio_int(300.0) / ratio_int(500.0), ratio_int(1000.0) / ratio_int(500.0),
         ratio_int(1000.0) / ratio_int(300.0)))
print('  T 从 500 降到 300 eV 时 R 的降幅 = %.6g%%'
      % (100.0 * (1.0 - ratio_int(300.0) / ratio_int(500.0))))
print('  主工作点 T=%.6g eV 时 dE/T = %.6g（强度比几乎不随温度变）' % (Te, dE_gap / Te))
for lam, m_u, T_i in ((hc_eVnm / dE_Lya, mC_u, 100.0), (lam_Da, mD_u, 100.0), (lam_HeII, mHe_u, 200.0)):
    dl = doppler_nm(lam, T_i, m_u)
    print('  lambda=%.6g nm, T_i=%.6g eV：dl/lambda = %.6g，速度展宽 = %.6g m/s （%.4g c）'
          % (lam, T_i, dl / lam, c * dl / lam, dl / lam))
print('  塞曼/多普勒 两条线的比值再比 = %.6g'
      % ((dZ_Da / doppler_nm(lam_Da, 100.0, mD_u)) / (dZ_C / doppler_nm(hc_eVnm / dE_Lya, 100.0, mC_u))))
print('  n_e 涨 10 倍时斯塔克线宽倍数 10^(2/3) = %.6g；Stark(1e21)/Stark(1e20) = %.6g'
      % (10.0 ** (2.0 / 3.0), (1.0e21 / 1.0e20) ** (2.0 / 3.0)))
print('  边缘 0.5 eV 时的中性份额 = %.6g（= %.6g%%）'
      % (1.0 - ion_frac(0.5, ne_edge, chi_H)[1], 100.0 * (1.0 - ion_frac(0.5, ne_edge, chi_H)[1])))
print('  C VI Ly-alpha 与边缘 T_e 之比 = %.6g' % (dE_Lya / Te_sol))
print('  碳与氢在同一温度下的电离度之比：T=0.5 eV -> %.6g'
      % (ion_frac(0.5, ne_edge, chi_C1)[1] / ion_frac(0.5, ne_edge, chi_H)[1]))
# ---- [AH.9] 练习与表格里的比值 ----
print('')
print('[AH.9] 练习与表格比值')
f_07_20 = ion_frac(0.7, 1.0e20, chi_H)[1]
f_07_21 = ion_frac(0.7, 1.0e21, chi_H)[1]
print('  T=0.7 eV：n=1e20 电离度 %.6g%% / n=1e21 电离度 %.6g%% = %.6g 倍'
      % (100.0 * f_07_20, 100.0 * f_07_21, f_07_20 / f_07_21))
L_tau1 = 1.0 / (1.0e18 * sigma_line)
print('  解答 5：tau=1 需要的长度 L = %.6g m；装置小半径 a=2 m / L = %.6g' % (L_tau1, 2.0 / L_tau1))
print('  解答 6：钨/碳质量比 = %.6g；其平方根 = %.6g；波长比 5 nm / C VI Ly-alpha = %.6g'
      % (mW_u / mC_u, (mW_u / mC_u) ** 0.5, 5.0 / (hc_eVnm / dE_Lya)))
print('  解答 1：He II 200 eV 线宽 / C VI 200 eV 线宽 = %.6g'
      % (doppler_nm(lam_HeII, 200.0, mHe_u) / doppler_nm(hc_eVnm / dE_Lya, 200.0, mC_u)))
print('  AH.6 表格：n_1=1e16、L=0.01/1/10 m 的 tau 与 1/3 之比 = %.6g / %.6g / %.6g'
      % (1e16 * sigma_line * 0.01 / (1.0 / 3.0), 1e16 * sigma_line * 1.0 / (1.0 / 3.0),
         1e16 * sigma_line * 10.0 / (1.0 / 3.0)))
print('以上为附录 AH 正文与门禁 P40 可以引用的全部数值。')
