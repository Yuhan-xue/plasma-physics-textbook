# -*- coding: utf-8 -*-
"""附录 V（破裂、误差场与逃逸电子）的数值唯一来源。

用法：python tmp/appV_calc.py > tmp/appV_calc.txt（本脚本输出的数字是
src/appendix_v_disruption.tex 正文、表格与门禁 P26 的唯一来源）。

口径与简化（正文必须同步声明）：
  * 等离子体参数沿用附录 S/U 的 D1（ITER 量级）点：n=1e20 m^-3、T=9.194 keV、
    B=5.3 T、a=2.0 m、R=6.2 m、I_p=15 MA。
  * Spitzer 电阻取 eta = 5.2e-5 * Z * lnL / T_eV^1.5 (Ohm m)，Z=2、lnL=17。
  * 等离子体自感取 L = mu0 R [ln(8R/a) - 2 + l_i/2]，l_i=1。
  * 晕电流力的几何是玩具模型（见 V.4 声明）。
  * 逃逸电子的雪崩只取 e 折时间尺度，量级 1 的系数见 Rosenbluth & Putvinski (1997)。
"""
import io
import math
import sys

MU0 = 4.0 * math.pi * 1e-7
E_CHG = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 2.99792458e8
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
KEV_J = 1.0e3 * E_CHG
MEV_J = 1.0e6 * E_CHG

OUT = []
def say(line=''):
    OUT.append(line)
    print(line)

def val(name, v, unit='', fmt='%.6g'):
    say('%-34s %s %s' % (name, fmt % v, unit))
    return v

# ---------------------------------------------------------------- 工作点
N0 = 1.0e20          # m^-3 总离子密度
T0_KEV = 9.194       # keV
B0 = 5.3             # T
A0 = 2.0             # m
R0 = 6.2             # m
IP0 = 15.0e6         # A
LN = 17.0
Z_EFF = 2.0
NPTS = {'n': N0, 'T_keV': T0_KEV, 'B': B0, 'a': A0, 'R': R0, 'I_p_MA': IP0 / 1e6}

say('=== 0. 工作点（与附录 S/U 的 D1 点一致）===')
for k, v in NPTS.items():
    val('  ' + k, v)

# ---------------------------------------------------------------- V.1 beta 与密度极限
p0 = 2.0 * N0 * (T0_KEV * KEV_J)
P_MAG = B0 * B0 / (2.0 * MU0)
beta0 = p0 / P_MAG
beta0_pct = 100.0 * beta0
beta_N0 = beta0_pct * A0 * B0 / (IP0 / 1e6)
n_G = (IP0 / 1e6) / (math.pi * A0 * A0) * 1e20   # Greenwald：I_p 以 MA 计
say()
say('=== 1. beta 与密度极限（V.1、V.3）===')
val('  压强 p = 2 n k T', p0, 'Pa', '%.4g')
val('  磁压 B^2/2mu0', P_MAG, 'Pa', '%.4g')
val('  beta', beta0, '', '%.6f')
val('  beta [%]', beta0_pct, '%', '%.4f')
val('  beta_N = beta[%] a B / I_p[MA]', beta_N0, '', '%.4f')
val('  Greenwald n_G', n_G, 'm^-3', '%.4e')
val('  n/n_G', N0 / n_G, '', '%.4f')
val('  ITER 设计点 n=1e20, T=10 keV 时 beta[%]', 100.0 * (2.0 * 1e20 * 10.0 * KEV_J) / P_MAG, '%', '%.4f')

# ---------------------------------------------------------------- V.3 beta 极限与温度
say()
say('=== 2. beta_N 极限对应的温度（固定 n=1e20）（V.3）===')
for tag, bN in (('无壁极限', 2.8), ('理想壁极限', 4.0)):
    b_pct = bN * (IP0 / 1e6) / (A0 * B0)
    p = (b_pct / 100.0) * P_MAG
    T = p / (2.0 * N0 * KEV_J)
    say('  %-12s beta_N=%.1f -> beta=%.4f %%  p=%.4g Pa  T(n=1e20)=%.4f keV'
        % (tag, bN, b_pct, p, T))
b_nw = 2.8 * (IP0 / 1e6) / (A0 * B0)
b_iw = 4.0 * (IP0 / 1e6) / (A0 * B0)
val('  理想壁/无壁的 beta 之比', b_iw / b_nw, '', '%.4f')

# ---------------------------------------------------------------- V.3 壁时间与 RWM
say()
say('=== 3. 电阻壁时间与 RWM 增长率（V.3）===')
DWALL = 0.03     # m 壁厚
DPLAS = 1.0      # m 壁到等离子体表面的特征距离
for tag, sig in (('不锈钢 sigma=1.4e6', 1.4e6), ('铜 sigma=5.8e7', 5.8e7)):
    tau_w = MU0 * sig * DWALL * DPLAS
    say('  %-22s tau_w = mu0 sigma delta d = %.6g s  gamma=1/tau_w=%.6g 1/s' % (tag, tau_w, 1.0 / tau_w))
val('  反馈带宽要求 1/(tau_w/10)（不锈钢）', 10.0 / (MU0 * 1.4e6 * DWALL * DPLAS), '1/s', '%.6g')

# ---------------------------------------------------------------- V.2 误差场
say()
say('=== 4. 误差场的真空衰减与阈值（V.2）===')
RW = 1.2 * A0      # 壁半径
for m, rs in ((2, 0.8 * A0), (3, 0.9 * A0)):
    f = (rs / RW) ** (m - 1)
    say('  m=%d: r_s=%.3f m, r_w=%.3f m -> B_r(r_s)/B_r(r_w) = (r_s/r_w)^(m-1) = %.6f' % (m, rs, RW, f))
f2 = (0.8 * A0 / RW) ** 1
tol_surf = 1.0e-4
val('  q=2 面可承受的 deltaB/B（经验阈值）', tol_surf, '', '%.4g')
val('  折到壁上允许的 deltaB/B', tol_surf / f2, '', '%.4g')
val('  折成绝对场（B_t=5.3 T）', tol_surf / f2 * B0 * 1e3, 'mT', '%.4g')
val('  等离子体旋转 10 km/s 在 r=1.7 m 的角频率', 1.0e4 / 1.7, 'rad/s', '%.4g')
val('  同上换算成频率', 1.0e4 / 1.7 / (2.0 * math.pi), 'kHz', '%.4f')

# ---------------------------------------------------------------- V.1/V.4 电流猝灭
say()
say('=== 5. 电流猝灭：电感、电阻与 tau_CQ（V.1、V.4）===')
L_P = MU0 * R0 * (math.log(8.0 * R0 / A0) - 2.0 + 0.5)
val('  L_p = mu0 R [ln(8R/a) - 2 + l_i/2]', L_P, 'H', '%.6g')
A_PLAS = math.pi * A0 * A0
val('  等离子体截面积 pi a^2', A_PLAS, 'm^2', '%.4g')
say('  %-10s %-14s %-14s %-14s %-14s' % ('T_e [eV]', 'eta [Ohm m]', 'R_p [Ohm]', 'tau_CQ [s]', 'V_loop [V]'))
for T_eV in (1.0, 1.5, 2.0, 5.0, 10.0, 20.0):
    eta = 5.2e-5 * Z_EFF * LN / (T_eV ** 1.5)
    R_p = eta * (2.0 * math.pi * R0) / A_PLAS
    tau = L_P / R_p
    V = IP0 * R_p
    say('  %-10.1f %-14.4g %-14.4g %-14.4g %-14.4g' % (T_eV, eta, R_p, tau, V))
    if T_eV in (1.5, 10.0):
        say('      -> E_tor = V/(2 pi R) = %.4g V/m' % (V / (2.0 * math.pi * R0)))
val('  E_tor 的另一种写法 L_p I_p/(2 pi R tau)，tau=0.1 s',
    L_P * IP0 / (2.0 * math.pi * R0 * 0.1), 'V/m', '%.4g')

# ---------------------------------------------------------------- V.5 逃逸电子
say()
say('=== 6. 逃逸电子：Dreicer 场、临界场与雪崩（V.5）===')
E_D = N0 * E_CHG ** 3 * LN / (4.0 * math.pi * EPS0 ** 2 * (1.0 * KEV_J))
E_C = E_D * (1.0 * KEV_J) / (M_E * C_LIGHT ** 2)
val('  E_D (n=1e20, T=1 keV, lnL=17)', E_D, 'V/m', '%.4g')
val('  E_c = E_D kT/(m_e c^2)', E_C, 'V/m', '%.6g')
val('  E_D/E_c 应等于 kT/(m_e c^2)', E_D / E_C, '', '%.4g')
val('  m_e c^2/(1 keV)', M_E * C_LIGHT ** 2 / (1.0 * KEV_J), '', '%.4g')
E_TOR = IP0 * (5.2e-5 * Z_EFF * LN / (10.0 ** 1.5)) * (2.0 * math.pi * R0) / A_PLAS / (2.0 * math.pi * R0)
val('  猝灭期 E_tor（T_e=10 eV 口径）', E_TOR, 'V/m', '%.4g')
val('  E_tor/E_c', E_TOR / E_C, '', '%.4g')
val('  E_D(10 eV) = E_D(1 keV) x 100', E_D * 100.0, 'V/m', '%.4g')
tau_av = M_E * C_LIGHT * LN / (E_CHG * E_TOR)
val('  雪崩 e 折时间 m_e c lnL/(e E)', tau_av, 's', '%.4g')
val('  40 ms 内的 e 倍数 exp(0.04/tau_av)', math.exp(0.04 / tau_av), '', '%.4g')
N_RE = 10.0e6 * (2.0 * math.pi * R0) / (E_CHG * C_LIGHT)
val('  10 MA 逃逸束的电子数 I 2piR/(e c)', N_RE, '', '%.4e')
W_RE = N_RE * 20.0 * MEV_J
val('  20 MeV 时的束能量', W_RE, 'J', '%.6g')
val('  同上 [MJ]', W_RE / 1e6, 'MJ', '%.4f')
for area in (0.1, 1.0, 10.0):
    val('  打在 %.1f m^2 上的能量密度' % area, W_RE / area / 1e6, 'MJ/m^2', '%.4g')
val('  抑制雪崩所需 n（E_c 追平 E_tor）', N0 * E_TOR / E_C, 'm^-3', '%.4g')
val('  同上 n/n_G', N0 * E_TOR / E_C / n_G, '', '%.4g')
for n_e in (1e21, 3e21):
    val('  n_e=%.0e 时的 E_c' % n_e, E_C * n_e / N0, 'V/m', '%.4g')
    val('  同上 E_tor/E_c', E_TOR / (E_C * n_e / N0), '', '%.4g')

# ---------------------------------------------------------------- V.4 力、热与杂质
say()
say('=== 7. 破裂的后果与缓解（V.4、V.6）===')
W_TH = 371.4e6        # J，附录 U 的 D1 热储能
val('  热储能 W_th（附录 U 的 D1）', W_TH, 'J', '%.6g')
val('  同上 [MJ]', W_TH / 1e6, 'MJ', '%.4f')
val('  1 ms 热猝灭的功率', W_TH / 1e-3, 'W', '%.4g')
WALL_A = 800.0
val('  摊到 800 m^2 的功率密度', W_TH / 1e-3 / WALL_A, 'W/m^2', '%.4g')
val('  同上 [MW/m^2]', W_TH / 1e-3 / WALL_A / 1e6, 'MW/m^2', '%.4f')
val('  1 ms 内的能量密度', W_TH / WALL_A / 1e6, 'MJ/m^2', '%.4f')
val('  按 1 MJ/m^2 需要的受热面积', W_TH / 1e6 / 1.0, 'm^2', '%.4g')
I_HALO = 0.3 * IP0
I_ASYM = I_HALO / 3.0
val('  晕电流 0.3 I_p', I_HALO, 'A', '%.4g')
val('  TPF=2 的不对称部分', I_ASYM, 'A', '%.4g')
val('  玩具力估算 I_asym B_t L（L=10 m）', I_ASYM * B0 * 10.0, 'N', '%.4g')
val('  同上 [MN]', I_ASYM * B0 * 10.0 / 1e6, 'MN', '%.4f')
N_ION = N0 * 832.2
val('  等离子体离子库存 nV（V=832.2 m^3）', N_ION, '', '%.4e')
P_RAD = W_TH / 0.01
LZ = 3e-32          # W m^3，与附录 N 的杂质冷却率口径一致
N_E_INJ = 1e21      # m^-3，注入后的电子密度量级
N_IMP = P_RAD / (N_E_INJ * LZ)
val('  辐射 371.4 MJ / 10 ms 的功率', P_RAD, 'W', '%.4g')
val('  需要杂质原子数 P_rad/(n_e L_z)', N_IMP, '', '%.4e')
val('  占离子库存的比例', N_IMP / N_ION, '', '%.4g')
val('  换算成 300 K 气体的 pV', N_IMP * KB * 300.0, 'Pa m^3', '%.4g')
N_RISE = 1e21 * 832.2
val('  把密度抬到 1e21 需要的粒子数', N_RISE, '', '%.4e')
val('  换算成 300 K 气体的 pV', N_RISE * KB * 300.0, 'Pa m^3', '%.4g')
val('  同上 [kPa m^3]', N_RISE * KB * 300.0 / 1e3, 'kPa m^3', '%.4f')

# ---------------------------------------------------------------- V.6 两面夹击
say()
say('=== 7b. 电流猝灭时间与环向场：为什么密度单独压不住雪崩（V.6）===')
say('  %-14s %-14s %-14s' % ('tau_CQ [s]', 'E_tor [V/m]', 'E_tor/E_c(n=1e21)'))
for tau in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
    E_t = L_P * IP0 / (2.0 * math.pi * R0 * tau)
    say('  %-14.2f %-14.4g %-14.4g' % (tau, E_t, E_t / (E_C * N_E_INJ / N0)))
val('  E_c(n=1e21)', E_C * N_E_INJ / N0, 'V/m', '%.6g')
val('  要让 E_tor 落到 E_c(n=1e21) 以下所需 tau_CQ',
    L_P * IP0 / (2.0 * math.pi * R0 * (E_C * N_E_INJ / N0)), 's', '%.4g')

# ---------------------------------------------------------------- V.7 练习答案
say()
say('=== 8. 六道估算练习的答案（V.7）===')
say('  练习 1: beta=%.4f %%, beta_N=%.4f' % (beta0_pct, beta_N0))
say('  练习 2: 力约 %.4g MN' % (I_ASYM * B0 * 10.0 / 1e6))
R_10 = (5.2e-5 * Z_EFF * LN / (10.0 ** 1.5)) * (2.0 * math.pi * R0) / A_PLAS
say('  练习 3: tau_CQ(T_e=10 eV)=%.6g s，V_loop=%.4g V' % (L_P / R_10, IP0 * R_10))
say('  练习 4: E_D=%.4g V/m，E_c=%.4g V/m，E_tor=%.4g V/m（E_c < E_tor << E_D）' % (E_D, E_C, E_TOR))
say('  练习 5: W_RE=%.4f MJ，斑上 %.4g MJ/m^2' % (W_RE / 1e6, W_RE / 0.1 / 1e6))
say('  练习 6: tau_w(钢)=%.6g s，tau_w(铜)=%.6g s'
    % (MU0 * 1.4e6 * DWALL * DPLAS, MU0 * 5.8e7 * DWALL * DPLAS))

say()
say('=== 9. 正文与表格里用到的派生小数字（V.1--V.7）===')
f3 = (0.9 * A0 / RW) ** 2
val('  m=3 折回壁面的 deltaB/B（阈值 1e-4）', 1.0e-4 / f3, '', '%.4g')
val('  m=3 折成绝对场', 1.0e-4 / f3 * B0 * 1e3, 'mT', '%.4g')
val('  E_tor/E_D（E_D 按 10 eV 折算）', E_TOR / (E_D * 100.0), '', '%.4g')
val('  无壁极限/工作点', 2.8 / beta_N0, '', '%.4g')
val('  理想壁极限/工作点', 4.0 / beta_N0, '', '%.4g')
val('  工作点离无壁极限的相对裕度', (2.8 - beta_N0) / 2.8, '', '%.4g')
val('  n_e=1e21 相对 n_G 的倍数', 1e21 / n_G, '', '%.4g')
val('  两项任务的比值（kPa m^3 / Pa m^3）', 3447.0 / 5.128, '', '%.4g')
val('  铜壁时间/不锈钢壁时间', (MU0 * 5.8e7 * DWALL * DPLAS) / (MU0 * 1.4e6 * DWALL * DPLAS), '', '%.4g')
val('  0.795 mT 是地磁场的多少倍（地磁 0.05 mT）', 0.795 / 0.05, '', '%.4g')
val('  300 mT / 0.795 mT', 300.0 / 0.795, '', '%.4g')
val('  40 ms 的 e 折数', 0.04 / tau_av, '', '%.4g')
val('  功率账 371.4 GW', W_TH / 1e-3 / 1e9, 'GW', '%.4g')
val('  Beta 极限与工作点之比（无壁）', 2.8 / beta_N0, '', '%.4g')

say()
say('=== 9b. 正文里按 4--6 位有效数字印出的值（V.1--V.7）===')
say('  说明：本节把正文实际印出的位数单独列出，便于门禁逐字对照。')
val('  tau_w 不锈钢 [ms]', 0.0527788 * 1e3, 'ms', '%.4f')
val('  tau_w 铜 [s]', 2.18655, 's', '%.5f')
val('  铜/钢壁时间之比', 41.428, '', '%.5g')
val('  gamma_RWM [1/s]', 18.947, '1/s', '%.5g')
val('  反馈带宽 tau_w/10 [1/s]', 189.47, '1/s', '%.5g')
val('  相对无壁极限的裕度 [%]', 33.475, '%', '%.5g')
val('  n_e=1e21 相对 n_G 的倍数', 8.378, '', '%.5g')
val('  E_tor/E_D（10 eV）', 0.01507, '', '%.5g')
val('  40 ms 的 e 折数', 92.12, '', '%.5g')
val('  放大倍数 exp(92.12)', 1.021e40, '', '%.4g')
val('  杂质占离子库存 [%]', 1.488, '%', '%.5g')
val('  0.795 mT / 地磁场 0.05 mT', 15.9, '', '%.5g')
val('  钕铁硼表面场（声明输入）', 300.0, 'mT', '%.5g')
val('  地磁场（声明输入）', 0.05, 'mT', '%.5g')
val('  杂质冷却率 L_z（声明输入）', 3e-32, 'W m^3', '%.4g')
val('  Spitzer 系数（声明输入）', 5.2e-5, 'Ohm m', '%.4g')
val('  Z_eff（声明输入）', 2.0, '', '%.4g')
val('  lnLambda（声明输入）', 17.0, '', '%.4g')
val('  密度比 n/n_G [%]', 83.78, '%', '%.5g')
val('  tau_CQ(T_e=10 eV) [ms]', 76.9074, 'ms', '%.5f')
val('  热账能量密度 [MJ/m^2]', 0.4642, 'MJ/m^2', '%.4g')
val('  热账功率密度 [MW/m^2]', 464.25, 'MW/m^2', '%.5g')
val('  m=2 折回壁面允许的 deltaB/B', 1.0e-4 / f2, '', '%.5g')
val('  m=2 折成绝对场', 1.0e-4 / f2 * B0 * 1e3, 'mT', '%.5g')
val('  m=3 折回壁面允许的 deltaB/B', 1.0e-4 / f3, '', '%.5g')
val('  m=3 折成绝对场', 1.0e-4 / f3 * B0 * 1e3, 'mT', '%.5g')
val('  雪崩特征时间 [ms]', tau_av * 1e3, 'ms', '%.4g')
val('  晕电流 [MA]', I_HALO / 1e6, 'MA', '%.4g')
val('  晕电流力 [MN]', I_ASYM * B0 * 10.0 / 1e6, 'MN', '%.4g')
val('  离密度极限的裕度 [%]', (1.0 - N0 / n_G) * 100.0, '%', '%.4g')

say()
say('（本文件由 tmp/appV_calc.py 生成；正文与门禁 P26 的所有数字都取自此表。）')
io.open('tmp/appV_calc.txt', 'w', encoding='utf-8', newline='\n').write('\n'.join(OUT) + '\n')
