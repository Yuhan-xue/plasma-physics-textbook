# -*- coding: utf-8 -*-
"""附录 AC（等离子体数值方法：从 PIC 到 MHD）的唯一数值来源。

工作点（DECLARED，与附录 M/S/U/V/W/X/Y/Z/AA/AB 同一 D1 参数）：
  等离子体：n_e = 1e20 m^-3、T_e = 9.194 keV、B = 5.3 T、Z_eff = 2、m_i = 2.5 m_p（D-T 混合）
  几何：R = 6.2 m、a = 2.0 m
  模拟参数（声明）：盒子边长 100 个德拜长度、网格步长 = 1 个德拜长度、每格 20 个粒子、
                    目标时长 1000 个电子等离子体周期、MHD 网格 1 cm / 256^3
所有打印值即正文与门禁可以引用的值。
"""
import math

# ---------------- 常数 ----------------
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
mp = 1.67262192369e-27
eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6
keV = 1.0e3 * e

# ---------------- DECLARED ----------------
n_e = 1.0e20                     # m^-3
T_e = 9.194 * keV                # J
B0 = 5.3                         # T
Zeff = 2.0
A_ion = 2.5                      # m_i = 2.5 m_p
R0 = 6.2                         # m
a0_minor = 2.0                   # m
L_box_lambdaD = 100.0            # 盒边长 = 100 λ_De
dx_lambdaD = 1.0                 # 网格步长 = 1 λ_De
Nppc = 20.0                      # 每格粒子数
T_target_wpe = 1000.0            # 目标时长 = 1000 / ω_pe
dt_frac_cyc = 0.2                # Δt = 0.2 / ω_ce
dx_mhd = 1.0e-2                  # m（MHD 网格）
N_mhd = 256                      # 每维格点数
CFL = 0.5
n_avg = 100.0                    # 时间平均步数（AC.4）
rate_per_core = 1.0e8            # 粒子步/秒/核（声明）
rate_gpu = 1.0e9                 # 粒子步/秒/GPU（声明）
mem_per_cell = 8 * 8             # 每格 8 个双精度量

DAO = '《等离子体物理自学教材》附录 AC 数值脚本'
print('=' * 72)
print(DAO, '—— 全部数字的唯一来源')
print('=' * 72)

print('\n[DECLARED] 输入一览')
print('  n_e = %.6g m^-3' % n_e)
print('  T_e = %.6g keV = %.6g J' % (T_e / keV, T_e))
print('  B = %.6g T' % B0)
print('  Z_eff = %.6g, m_i = %.6g m_p' % (Zeff, A_ion))
print('  几何 R = %.6g m, a = %.6g m' % (R0, a0_minor))
print('  模拟盒 = %.6g λ_De, Δx = %.6g λ_De, 每格 %.6g 粒子' % (L_box_lambdaD, dx_lambdaD, Nppc))
print('  目标时长 = %.6g / ω_pe, Δt = %.6g / ω_ce' % (T_target_wpe, dt_frac_cyc))
print('  MHD：Δx = %.6g m, %d^3 网格, CFL = %.6g' % (dx_mhd, N_mhd, CFL))

# ---------------- AC.1 分辨率 ----------------
m_i = A_ion * mp
lam_De = math.sqrt(eps0 * T_e / (n_e * e ** 2))
wpe = math.sqrt(n_e * e ** 2 / (eps0 * me))
fpe = wpe / (2 * math.pi)
wce = e * B0 / me
fce = wce / (2 * math.pi)
N_D = n_e * lam_De ** 3
v_the = math.sqrt(T_e / me)
dt_wpe = 0.2 / wpe
dt_wce = dt_frac_cyc / wce
dt_use = dt_wce
T_ce = 2 * math.pi / wce
print('\n[AC.1] 分辨率：三条判据')
print('  lambda_De = %.6g m = %.4g um' % (lam_De, lam_De * 1e6))
print('  omega_pe = %.6g rad/s = %.5g Hz' % (wpe, fpe))
print('  omega_ce = %.6g rad/s = %.5g Hz' % (wce, fce))
print('  德拜球内粒子数 N_D = n λ_De^3 = %.6g' % N_D)
print('  电子速度分布宽度 sigma_ve = %.6g m/s = sigma_ve/c = %.5g' % (v_the, v_the / c))
print('  Δt 上限（等离子体振荡 0.2/ω_pe）= %.6g s' % dt_wpe)
print('  Δt 上限（回旋 0.2/ω_ce）= %.6g s' % dt_wce)
print('  两者之比 (0.2/ω_pe)/(0.2/ω_ce) = %.6g' % (dt_wpe / dt_wce))
print('  采用 Δt = %.6g s；每步热位移 v_the Δt = %.6g m = %.5g λ_De'
      % (dt_use, v_the * dt_use, v_the * dt_use / lam_De))
print('  一个回旋周期的步数 T_ce/Δt = %.6g' % (T_ce / dt_use))
print('  目标模拟时长 = %.6g 个回旋周期' % (T_target_wpe / wpe / T_ce))
# 碰撞时标对照（Spitzer，lnΛ 声明为 17）
lnL = 17.0
tau_ei = (3.0 * (2 * math.pi) ** 1.5 * eps0 ** 2 * me ** 2 * v_the ** 3
          / (n_e * e ** 4 * Zeff * lnL))
print('  lnΛ = %.6g（声明）；电子-离子碰撞时间 tau_ei = %.6g s' % (lnL, tau_ei))
t_sim_ac1 = T_target_wpe / wpe
print('  碰撞频率 nu_ei = %.6g Hz' % (1.0 / tau_ei))
print('  目标时长 t_sim = %.6g s；t_sim/tau_ei = %.5g（模拟窗口内碰撞可忽略）'
      % (t_sim_ac1, t_sim_ac1 / tau_ei))

# ---------------- AC.2 推进粒子：相位误差 ----------------
err_step = (wce * dt_use) ** 2 / 24
err_cyc = err_step * (T_ce / dt_use)
gamma_1MeV = 1.0 + 1.0e6 * e / (me * c ** 2)
print('\n[AC.2] 蛙跳与 Boris：相位误差')
print('  无量纲步长 ω_ce Δt = %.6g' % (wce * dt_use))
print('  每步相位误差 (ω_ce Δt)^2/24 = %.6g rad' % err_step)
print('  每回旋累计相位误差 = %.6g rad = %.5g 度' % (err_cyc, math.degrees(err_cyc)))
print('  1 MeV 电子的 gamma = %.6g（相对论推进器的必要性）' % gamma_1MeV)
print('  D-T 离子回旋频率 ω_ci = %.6g rad/s' % (e * B0 / m_i))
print('  离子与电子回旋频率之比 ω_ci/ω_ce = %.6g' % ((e * B0 / m_i) / wce))
print('  离子/电子回旋周期之比 = %.6g（同一 Δt 对离子「过粗」的倍数）'
      % (wce / (e * B0 / m_i)))

# ---------------- AC.3 PIC 的成本 ----------------
L_box = L_box_lambdaD * lam_De
N_cell_1d = L_box / (dx_lambdaD * lam_De)
N_cell_3d = N_cell_1d ** 3
N_part = N_cell_3d * Nppc
t_sim = T_target_wpe / wpe
N_step = t_sim / dt_use
psteps = N_part * N_step
cpu_s = psteps / rate_per_core
gpu_s = psteps / rate_gpu
V_plasma = 832.2                  # m^3（D1 等离子体体积，与附录 M/S/U/V/W/X/Y/Z 同口径）
N_real = n_e * V_plasma
print('\n[AC.3] PIC 的成本')
print('  盒边长 L = %.6g m = %.5g cm' % (L_box, L_box * 1e2))
print('  一维格点数 = %.6g，三维格点数 = %.6g' % (N_cell_1d, N_cell_3d))
print('  总粒子数 = %.6g' % N_part)
print('  目标时长 t_sim = %.6g s（%.6g 个等离子体周期）' % (t_sim, T_target_wpe / (2 * math.pi)))
print('  时间步数 = %.6g' % N_step)
print('  粒子步总数 = %.6g' % psteps)
print('  单核估算（1e8 粒子步/秒）= %.6g s = %.5g 小时' % (cpu_s, cpu_s / 3600))
print('  GPU 估算（1e9 粒子步/秒）= %.6g s = %.5g 小时' % (gpu_s, gpu_s / 3600))
print('  真实 D1 等离子体体积 V = %.6g m^3 内粒子数 = %.6g' % (V_plasma, N_real))
print('  真实/模拟粒子数之比 = %.6g' % (N_real / N_part))

# ---------------- AC.4 形函数、噪声与时间平均 ----------------
def noise(npp):
    return 1.0 / math.sqrt(npp)

print('\n[AC.4] 粒子噪声与时间平均')
for npp in (20.0, 25.0, 100.0, 1000.0, 1.0e4):
    print('  每格 %.6g 个粒子 → 单格相对噪声 1/sqrt(N_p) = %.5g%%' % (npp, 100 * noise(npp)))
print('  时间平均 %.6g 步后噪声 = %.5g%%（降低 sqrt(M) = %.5g 倍）'
      % (n_avg, 100 * noise(Nppc) / math.sqrt(n_avg), math.sqrt(n_avg)))
for npp in (20, 25, 100, 1000, 10000):
    print('  每格 %.6g 个粒子 + 时间平均 %.6g 步 → 噪声 = %.5g%%'
          % (npp, n_avg, 100 * noise(npp) / math.sqrt(n_avg)))
print('  要达到 1%% 精度、不做时间平均 → 每格需要 %.6g 个粒子' % (1.0 / 0.01 ** 2))
print('  要达到 2%% 精度、不做时间平均 → 每格需要 %.6g 个粒子' % (1.0 / 0.02 ** 2))

# ---------------- AC.5 场求解与离散误差 ----------------
fft_ops = N_cell_3d * math.log2(N_cell_3d)
dir_ops = N_cell_3d ** 2
print('\n[AC.5] 场求解与离散误差')
print('  FFT 泊松：N log2 N = %.6g 次运算' % fft_ops)
print('  直接求和：N^2 = %.6g 次运算' % dir_ops)
print('  两者之比 = %.6g' % (dir_ops / fft_ops))
print('  二阶中心差分在 Δx = λ_De 时的相对误差量级 (Δx/λ_De)^2/12 = %.6g'
      % (dx_lambdaD ** 2 / 12.0))
print('  能量守恒验收线：总能量漂移 < 1%%；噪声主导时能量漂移 ~ %.5g%%（每格 %.6g 粒子）'
      % (100 * noise(Nppc), Nppc))
print('  收敛性：二阶格式误差 ∝ Δx²，Δx 减半 → 误差降到 %.6g 倍' % 0.5 ** 2)
for r in (2.0, 1.0, 0.5):
    print('    Δx/λ_De = %.6g → 截断误差量级 (Δx/λ_De)^2/12 = %.5g%%' % (r, 100 * r ** 2 / 12.0))

# ---------------- AC.6 流体与 MHD ----------------
rho = n_e * m_i
vA = B0 / math.sqrt(mu0 * rho)
c_s = math.sqrt(T_e / m_i)   # 离子声速（冷离子近似，与附录 AE 的 c_s 同口径）
v_ms = math.sqrt(vA ** 2 + c_s ** 2)
tau_A = R0 / vA
dt_mhd = CFL * dx_mhd / v_ms
n_step_mhd = 10 * tau_A / dt_mhd
mem = N_mhd ** 3 * mem_per_cell
print('\n[AC.6] 流体与 MHD')
print('  rho = n_e m_i = %.6g kg/m^3' % rho)
print('  v_A = %.6g m/s = %.5g km/s' % (vA, vA / 1e3))
print('  c_s = %.6g m/s = %.5g km/s' % (c_s, c_s / 1e3))
print('  快磁声速 v_ms = %.6g m/s = %.5g km/s' % (v_ms, v_ms / 1e3))
print('  阿尔芬时间 tau_A = R/v_A = %.6g s = %.4g us' % (tau_A, tau_A * 1e6))
print('  Δt（CFL=%.6g、Δx=%.6g m）= %.6g s' % (CFL, dx_mhd, dt_mhd))
print('  模拟 10 tau_A 需要的步数 = %.6g' % n_step_mhd)
print('  %d^3 网格内存 = %.6g B = %.5g GB' % (N_mhd, mem, mem / 1e9))
print('  ∇·B 保持：约束输运（CT）与投影法都是把误差压到机器精度量级')
print('  显式 CFL 限制 vs 隐式：隐式步长可放大 %.3g 倍（声明，取 100 倍量级）' % 100.0)

# ---------------- AC.7 练习解答 ----------------
print('\n[AC.7] 六道练习解答')
# 1) n=1e19, T=1 keV
n1, T1 = 1.0e19, 1.0 * keV
lam1 = math.sqrt(eps0 * T1 / (n1 * e ** 2))
wpe1 = math.sqrt(n1 * e ** 2 / (eps0 * me))
ND1 = n1 * lam1 ** 3
print('  1) n=1e19 m^-3, T=1 keV: λ_De = %.6g m = %.5g um；ω_pe = %.6g rad/s；N_D = %.6g'
      % (lam1, lam1 * 1e6, wpe1, ND1))
# 2) Δx = 2 λ_De
print('  2) Δx = 2 λ_De：网格已分辨不了德拜屏蔽（Δx/λ_De = %.6g），会出现有限网格不稳定性；'
      '判据要求 Δx <= λ_De，即比值 <= 1' % 2.0)
# 3) B = 2 T
B3 = 2.0
wce3 = e * B3 / me
dt3 = 0.2 / wce3
print('  3) B=2 T: ω_ce = %.6g rad/s；Δt <= 0.2/ω_ce = %.6g s；每回旋步数 = %.6g'
      % (wce3, dt3, 2 * math.pi / 0.2))
# 4) 盒 50 λ_De, 每格 50 粒子, 3D, 1000/ω_pe
L4 = 50.0 * lam_De
n4 = L4 / lam_De
part4 = n4 ** 3 * 50.0
step4 = (1000.0 / wpe) / dt_use
print('  4) 盒 50 λ_De / 每格 50 粒子：格点 %.6g，粒子 %.6g；步数 %.6g；粒子步 %.6g'
      % (n4 ** 3, part4, step4, part4 * step4))
print('     单核（1e8 粒子步/秒）= %.6g s = %.4g 小时' % (part4 * step4 / 1e8, part4 * step4 / 1e8 / 3600))
# 5) MHD Δx = 5 mm
dx5 = 5.0e-3
dt5 = CFL * dx5 / v_ms
print('  5) MHD Δx=5 mm：Δt = CFL Δx/v_ms = %.6g s；模拟 10 τ_A 需 %.6g 步' % (dt5, 10 * tau_A / dt5))
# 6) 噪声：每格 25 粒子，要 2%
print('  6) 每格 25 粒子 → 噪声 %.5g%%；要 2%% 需每格 %.6g 个粒子（不做时间平均）'
      % (100 * noise(25.0), 1.0 / 0.02 ** 2))

print('\n' + '=' * 72)
print('以上为附录 AC 正文与门禁 P35 可以引用的全部数值')
print('=' * 72)
