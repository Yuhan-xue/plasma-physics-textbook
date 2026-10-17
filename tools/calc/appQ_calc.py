# -*- coding: utf-8 -*-
"""附录 Q（等离子体诊断）全部数值的唯一来源。
用法：python tmp/appQ_calc.py > tmp/appQ_calc.txt
常量与 tmp/appP_calc.py 保持同源；正文与门禁都只能引用本脚本的输出。
"""
import io
import math
import os

E = 1.602176634e-19
ME = 9.1093837015e-31
U = 1.66053906660e-27
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
MU0 = 1.25663706212e-6
C = 2.99792458e8
H = 6.62607015e-34
M_AR = 39.948 * U
SIGMA_T = 6.6524587321e-29
K_EV = E / KB          # 1 eV = 11604.518 K

OUT = []


def say(s=''):
    OUT.append(s)


def debye(te, ne):
    """德拜长度，te 以 eV，ne 以 m^-3。"""
    return math.sqrt(EPS0 * te * E / (ne * E * E))


def cs(te, mi):
    """玻姆速度（离子声速），te 以 eV。"""
    return math.sqrt(te * E / mi)


def vth_e(te):
    """电子平均热速度 sqrt(8kT/pi m)。"""
    return math.sqrt(8.0 * te * E / (math.pi * ME))


def doppler_pm(lam_nm, t_k, m_amu):
    """Doppler FWHM（pm）。"""
    return 7.16e-7 * lam_nm * math.sqrt(t_k / m_amu) * 1000.0


def n_crit(f):
    """截止密度（m^-3），f 以 Hz。"""
    return EPS0 * ME * (2.0 * math.pi * f) ** 2 / (E * E)


say('=== 附录 Q 数值（tmp/appQ_calc.py 生成）===')
say('常量：E=%.9e ME=%.10e U=%.11e EPS0=%.10e KB=%.7e MU0=%.11e C=%.8e H=%.8e'
    % (E, ME, U, EPS0, KB, MU0, C, H))
say('常数：M_AR=%.6e  SIGMA_T=%.10e  1 eV=%.6f K' % (M_AR, SIGMA_T, K_EV))

# ---------------------------------------------------------------- Q.1 参数地图
say()
say('--- Q.1 参数地图（附录 P 的五个工况，供诊断手段对照）---')
regimes = [
    ('集成电路刻蚀(CCP/ICP)', 3.0, 1.0e17),
    ('溅射与PECVD', 2.0, 1.0e16),
    ('等离子体医学(射流)', 1.5, 1.0e19),
    ('大气压DBD微放电', 2.0, 1.0e21),
    ('聚变边界(对照,附录N)', 15.0, 2.0e19),
]
for name, te, ne in regimes:
    ld = debye(te, ne)
    say('  %-24s Te=%.1f eV  ne=%.1e m^-3  lambda_D=%.4e m = %.2f um'
        % (name, te, ne, ld, ld * 1e6))

# ------------------------------------------------- Q.2 探针几何与 OML 判据
say()
say('--- Q.2 探针半径与德拜长度之比（判据）---')
rp = 0.5e-3
say('  圆柱探针半径 r_p = %.1f mm = %.0f um' % (rp * 1e3, rp * 1e6))
for ne, te in [(1.0e17, 3.0), (1.0e16, 3.0), (1.0e15, 3.0)]:
    ld = debye(te, ne)
    say('  ne=%.0e m^-3 Te=%.1f eV -> lambda_D=%.2f um, r_p/lambda_D=%.2f'
        % (ne, te, ld * 1e6, rp / ld))
ld_etch = debye(3.0, 1.0e17)
say('  OML 经验上限 r_p/lambda_D ~ 3 时允许的最大半径（刻蚀工况）= %.1f um' % (3.0 * ld_etch * 1e6))
say('  刻蚀工况下 OML 要求 r_p <= %.2f um（比 0.5 mm 细 %.0f 倍）'
    % (3.0 * ld_etch * 1e6, rp / (3.0 * ld_etch)))

# ------------------------------------------------ Q.3 由 I-V 曲线斜率读 T_e
say()
say('--- Q.3 单探针 I-V 曲线：对数斜率法读 T_e ---')
TE3 = 3.0
IES0 = 10.0e-3
IIS = 1.0e-3
say('  合成数据：Te=%.1f eV, 电子饱和电流 I_es0=%.1f mA, 离子饱和电流 I_is=%.1f mA'
    % (TE3, IES0 * 1e3, IIS * 1e3))
say('  理论：I_e = I_es0*exp((V-V_p)/Te)  =>  d(ln I_e)/dV = 1/Te')
say('  %6s %12s %12s' % ('V-Vp/V', 'I_e/mA', 'I_net/mA'))
for dv in [0.0, -1.0, -2.0, -3.0, -4.0, -6.0, -8.0]:
    ie = IES0 * math.exp(dv / TE3)
    say('  %6.1f %12.4f %12.4f' % (dv, ie * 1e3, (ie - IIS) * 1e3))
ie_a = IES0 * math.exp(-3.0 / TE3)
ie_b = IES0 * math.exp(-6.0 / TE3)
slope_raw = (math.log(ie_a - IIS) - math.log(ie_b - IIS)) / 3.0
slope_ok = (math.log(ie_a) - math.log(ie_b)) / 3.0
say('  取 V-Vp = -3 V 与 -6 V 两点：')
say('    正确（先减 I_is）：slope=%.6f /V -> Te=%.4f eV' % (slope_ok, 1.0 / slope_ok))
say('    错误（不减 I_is）：slope=%.6f /V -> Te=%.4f eV，相对误差 %+.2f%%'
    % (slope_raw, 1.0 / slope_raw, (1.0 / slope_raw - TE3) / TE3 * 100.0))
ie_c = IES0 * math.exp(-1.0 / TE3)
ie_d = IES0 * math.exp(-2.0 / TE3)
s2 = (math.log(ie_c - IIS) - math.log(ie_d - IIS)) / 1.0
say('  取更靠近 V_p 的 -1 V 与 -2 V 两点（仍不减 I_is）：Te=%.4f eV，相对误差 %+.2f%%'
    % (1.0 / s2, (1.0 / s2 - TE3) / TE3 * 100.0))
say('  结论：离子电流占净电流的比例越小（越靠近 V_p），漏减造成的偏差越小。')

# ------------------------------------------- Q.4 由离子饱和电流读 n_e
say()
say('--- Q.4 离子饱和电流 -> 密度 ---')
r_p4, l_p4 = 0.5e-3, 10.0e-3
A_p4 = 2.0 * math.pi * r_p4 * l_p4
cs_ar3 = cs(3.0, M_AR)
ji_ar3 = E * 1.0e17 * cs_ar3
say('  圆柱探针 r_p=%.1f mm, l_p=%.1f mm -> 收集面积 A_p=%.6e m^2'
    % (r_p4 * 1e3, l_p4 * 1e3, A_p4))
say('  Te=3.0 eV 的氩离子玻姆速度 c_s=%.2f m/s' % cs_ar3)
say('  j_i = e*n_e*c_s = %.4f A/m^2 = %.4f mA/cm^2（与附录 P 的刻蚀基准一致）'
    % (ji_ar3, ji_ar3 * 0.1))
iis_thin = ji_ar3 * A_p4
iis_bohm = 0.6 * ji_ar3 * A_p4
say('  A_p 上的离子饱和电流：裸玻姆 %.6f mA，含 0.6 鞘层扩张因子 %.6f mA'
    % (iis_thin * 1e3, iis_bohm * 1e3))
I_meas = 1.0e-3
ne_ok = I_meas / (0.6 * E * cs_ar3 * A_p4)
ne_bad = I_meas / (E * cs_ar3 * A_p4)
say('  实测 I_is=%.3f mA -> n_e=%.4e m^-3（含 0.6）；漏掉 0.6 会得到 %.4e m^-3（低估 %.1f%%）'
    % (I_meas * 1e3, ne_ok, ne_bad, (1.0 - ne_bad / ne_ok) * 100.0))
say('  误差传递：d n_e/n_e = d I_is/I_is + d A_p/A_p + 0.5*d T_e/T_e')

# ------------------------------------------------- Q.5 悬浮电位
say()
say('--- Q.5 悬浮电位与等离子体电位 ---')
ve = vth_e(3.0)
ratio = ve / (4.0 * cs_ar3)
ratio06 = ratio / 0.6
say('  Te=3.0 eV：电子平均热速度 v_e=%.4e m/s' % ve)
say('  理想比（随机通量电子 : 玻姆离子）v_e/(4 c_s) = %.4f -> V_p-V_f = %.4f * Te = %.4f V'
    % (ratio, math.log(ratio), math.log(ratio) * 3.0))
say('  含 0.6 因子：比 = %.4f -> V_p-V_f = %.4f * Te = %.4f V'
    % (ratio06, math.log(ratio06), math.log(ratio06) * 3.0))
say('  两式之差 = %.4f V = Te*ln(1/0.6) = %.4f V'
    % (math.log(ratio06) * 3.0 - math.log(ratio) * 3.0, 3.0 * math.log(1.0 / 0.6)))
say('  与附录 P 的 e*phi_f/kTe = %.4f 一致' % math.log(ratio))

# ------------------------------------------------- Q.6 Rogowski 线圈
say()
say('--- Q.6 磁诊断：Rogowski 线圈与 RC 积分 ---')
N_r, A_c, r_mean = 1000, 1.0e-4, 0.05
M = MU0 * N_r * A_c / (2.0 * math.pi * r_mean)
say('  N=%d 匝, 单匝面积 %.1f cm^2, 平均半径 %.1f cm -> 互感 M=%.4e H = %.1f nH'
    % (N_r, A_c * 1e4, r_mean * 1e2, M, M * 1e9))
I_pk, t_rise = 100.0, 1.0e-6
dIdt = I_pk / t_rise
V_coil = M * dIdt
say('  放电 I=%.0f A 在 %.0f us 内上升 -> dI/dt=%.3e A/s -> 线圈输出 %.2f V' % (I_pk, t_rise * 1e6, dIdt, V_coil))
R_i, C_i = 50.0, 1.0e-6
RC = R_i * C_i
V_out = M * I_pk / RC
say('  RC 积分：R=%.0f ohm, C=%.1f uF -> 时间常数 RC=%.1f us（是脉宽的 %.0f 倍）'
    % (R_i, C_i * 1e6, RC * 1e6, RC / t_rise))
say('  积分后输出电压 V_out = M*I/(RC) = %.4f V；灵敏度 S=M/(RC)=%.4e V/A'
    % (V_out, M / RC))
say('  校验：V_out = V_coil*脉宽/RC = %.4f V' % (V_coil * t_rise / RC))

# ------------------------------------------------- Q.7 光学诊断
say()
say('--- Q.7 光学诊断：Doppler 宽度与谱线比 ---')
for t_k, tag in [(300.0, '室温 300 K'), (K_EV * 0.5, 'Ti=0.5 eV'), (K_EV * 3.0, 'Ti=3 eV')]:
    say('  Ar I 696.5 nm, %-12s -> Doppler FWHM = %.4f pm' % (tag, doppler_pm(696.5, t_k, 39.948)))
say('  常见光谱仪仪器宽度 50 pm（0.05 nm）-> 比室温 Doppler 宽度大 %.0f 倍，测不出离子温度'
    % (50.0 / doppler_pm(696.5, 300.0, 39.948)))
dE_oes = 0.40
say('  OES 双线比（阈值能差 %.2f eV）：I_1/I_2 = K*exp(-%.2f/Te)，取 K=1' % (dE_oes, dE_oes))
say('  %8s %10s' % ('强度比', 'Te/eV'))
for r in [0.95, 0.90, 0.875, 0.85, 0.75, 0.5]:
    say('  %8.3f %10.4f' % (r, dE_oes / (-math.log(r)) if r < 1.0 else float('inf')))
t_lo = dE_oes / (-math.log(0.850))
t_hi = dE_oes / (-math.log(0.900))
say('  比值 ±2.9%% 时 Te 从 %.4f 变到 %.4f eV（-%.1f%% / +%.1f%%），误差被放大近 10 倍'
    % (t_lo, t_hi, (1.0 - t_lo / (dE_oes / (-math.log(0.875)))) * 100.0,
       (t_hi / (dE_oes / (-math.log(0.875))) - 1.0) * 100.0))
P_oes = 1.0e-3
lam_oes = 750.0e-9
eff = 0.01
rate = P_oes * lam_oes / (H * C) * eff
say('  1 mW 的 750 nm 光进入 1%% 通量的光谱仪 -> %.4e 光子/s；曝光 1 ms 得 %.3e 个光子'
    % (rate, rate * 1e-3))

# ------------------------------------------------- Q.8 汤姆孙散射与微波
say()
say('--- Q.8 汤姆孙散射与微波干涉仪 ---')
lam_i = 532.0e-9
theta = math.pi / 2.0
k_s = 2.0 * (2.0 * math.pi / lam_i) * math.sin(theta / 2.0)
say('  入射 532 nm、散射角 90 deg -> k = 2 k_i sin(theta/2) = %.4e 1/m' % k_s)
for name, ne, te in [('刻蚀', 1.0e17, 3.0), ('大气压 DBD', 1.0e21, 2.0), ('聚变边界', 2.0e19, 15.0)]:
    ld = debye(te, ne)
    say('  %-12s lambda_D=%.4e m -> k*lambda_D=%.2f -> Salpeter alpha=1/(k*lambda_D)=%.4e'
        % (name, ld, k_s * ld, 1.0 / (k_s * ld)))
E_L, t_l = 1.0, 10.0e-9
E_ph = H * C / lam_i
N_ph = E_L / E_ph
n_e_t, L_t, dOm = 1.0e17, 0.01, 0.01
N_s = n_e_t * L_t * SIGMA_T * dOm * N_ph
say('  激光 %.1f J / %.0f ns：单发 %.4e 个 532 nm 光子（每个 %.4e J）' % (E_L, t_l * 1e9, N_ph, E_ph))
say('  n_e=%.1e m^-3, L=%.0f mm, 收集立体角 %.2f sr -> 散射光子 %.1f 个，统计误差 %.2f%%'
    % (n_e_t, L_t * 1e3, dOm, N_s, 100.0 / math.sqrt(N_s)))
f_mw, L_mw = 10.0e9, 0.3
nc10 = n_crit(f_mw)
say('  微波 %.0f GHz -> 截止密度 n_c=%.4e m^-3' % (f_mw * 1e-9, nc10))
say('  n_e=1.0e17 m^-3 时 n_e/n_c=%.6f；L=%.1f m 的相移 dphi=(omega/c)*(n_e/(2 n_c))*L=%.4f rad = %.4f 条纹'
    % (1.0e17 / nc10, L_mw, (2.0 * math.pi * f_mw / C) * (1.0e17 / (2.0 * nc10)) * L_mw,
       (2.0 * math.pi * f_mw / C) * (1.0e17 / (2.0 * nc10)) * L_mw / (2.0 * math.pi)))
ne_high = 1.0e19
say('  要测 n_e=%.1e m^-3（附录 P 的医学射流），截止密度必须更高：f > %.4e Hz = %.2f GHz'
    % (ne_high, math.sqrt(ne_high / n_crit(1.0e9)) * 1.0e9, math.sqrt(ne_high / n_crit(1.0e9))))

# ------------------------------------------------- Q.9 六道练习
say()
say('--- Q.9 六道练习的答案 ---')
te_ex1, ies_ex1, iis_ex1 = 2.0, 5.0e-3, 0.5e-3
ie1 = ies_ex1 * math.exp(-2.0 / te_ex1)
ie2 = ies_ex1 * math.exp(-4.0 / te_ex1)
sl_ok = (math.log(ie1) - math.log(ie2)) / 2.0
sl_bad = (math.log(ie1 - iis_ex1) - math.log(ie2 - iis_ex1)) / 2.0
say('  练习1: Te=%.1f eV, I_es0=%.1f mA, I_is=%.1f mA；两点 V-Vp=-2/-4 V'
    % (te_ex1, ies_ex1 * 1e3, iis_ex1 * 1e3))
say('     I_e=%.4f / %.4f mA -> 正确 Te=%.4f eV；漏减 I_is 得 Te=%.4f eV（%+.1f%%）'
    % (ie1 * 1e3, ie2 * 1e3, 1.0 / sl_ok, 1.0 / sl_bad, (1.0 / sl_bad - te_ex1) / te_ex1 * 100.0))
r_p5, l_p5, te5 = 0.4e-3, 8.0e-3, 2.5
A_p5 = 2.0 * math.pi * r_p5 * l_p5
cs5 = cs(te5, M_AR)
ne5 = 0.60e-3 / (0.6 * E * cs5 * A_p5)
say('  练习2: r_p=%.1f mm, l_p=%.1f mm, Te=%.1f eV -> A_p=%.6e m^2, c_s=%.2f m/s, I_is=0.60 mA -> n_e=%.4e m^-3'
    % (r_p5 * 1e3, l_p5 * 1e3, te5, A_p5, cs5, ne5))
ld5 = debye(4.0, 5.0e16)
say('  练习3: ne=5e16 m^-3, Te=4.0 eV -> lambda_D=%.4e m = %.2f um, r_p=1 mm 时 r_p/lambda_D=%.2f'
    % (ld5, ld5 * 1e6, 1.0e-3 / ld5))
say('     要回到 OML 适用区（<=3）需 r_p <= %.1f um' % (3.0 * ld5 * 1e6))
N_r2, A_c2, r_mean2 = 500, 2.0e-4, 0.04
M2 = MU0 * N_r2 * A_c2 / (2.0 * math.pi * r_mean2)
I2, t2 = 200.0, 2.0e-6
RC2 = 100.0 * 1.0e-6
say('  练习4: N=500, A_c=2 cm^2, r=4 cm -> M=%.4e H；I=200 A / 2 us -> V_coil=%.2f V；RC=%.0f us -> V_out=%.4f V, S=%.4e V/A'
    % (M2, M2 * I2 / t2, RC2 * 1e6, M2 * I2 / RC2, M2 / RC2))
ne5b, L5b, dOm5b, E5b = 5.0e18, 0.01, 0.02, 0.5
N_ph5 = E5b / E_ph
N_s5 = ne5b * L5b * SIGMA_T * dOm5b * N_ph5
say('  练习5: E_L=%.1f J, ne=%.1e m^-3, L=1 cm, dOmega=%.2f sr -> %.0f 个散射光子, 统计误差 %.3f%%'
    % (E5b, ne5b, dOm5b, N_s5, 100.0 / math.sqrt(N_s5)))
ne6 = 3.0e18
f6 = math.sqrt(ne6 / n_crit(1.0e9)) * 1.0e9
nc6 = n_crit(20.0e9)
phi6 = (2.0 * math.pi * 20.0e9 / C) * (ne6 / (2.0 * nc6)) * 0.2
say('  练习6: 要测 ne=%.1e m^-3 需 f > %.4e Hz = %.2f GHz；取 20 GHz 时 n_c=%.4e, ne/nc=%.4f,'
    % (ne6, f6, f6 * 1e-9, nc6, ne6 / nc6))
say('     L=0.2 m 的相移=%.4f rad = %.4f 条纹' % (phi6, phi6 / (2.0 * math.pi)))

say()
say('以上数字是正文与门禁的唯一来源。')

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'appQ_calc.txt')
io.open(path, 'w', encoding='utf-8', newline='\n').write('\n'.join(OUT) + '\n')
print('\n'.join(OUT))
print()
print('written:', path)
