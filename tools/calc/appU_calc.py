# -*- coding: utf-8 -*-
"""附录 U（聚变增益 Q 与点火）与门禁 P23 的数值唯一来源。

用法：python tmp/appU_calc.py > tmp/appU_calc.txt

口径（正文必须与之一致）：
  * 装置沿用附录 S/tmp/appT_calc.py 的三个自洽工作点（同一组输入，不另立一套）。
  * D-T 反应率：n_D = n_T = n/2，反应率 R = n^2 <sigma v>/4；
      <sigma v>_DT = 1.1e-24 T_keV^2  [m^3/s]（第 10 章卡片公式的米制口径）
      E_fus = 17.6 MeV = 2.819e-12 J，P_alpha = P_fus/5，P_neutron = (4/5) P_fus。
    该 T^2 近似适用于约 8--20 keV；低温端（T < 5 keV）它高估真实 Bosch-Hale 值约 2 倍，
    因此只在 U.5「Q 随 P_aux 下降」的趋势表里用它，并在正文明确标注。
  * 增益定义：Q = P_fus/P_aux；alpha 加热份额 f_alpha = P_alpha/(P_alpha+P_aux) = Q/(Q+5)。
  * 燃烧闭环（与附录 S 的口径差别，正文必须写清）：附录 S 把 P 当作外部加热功率；
    本附录的不动点用 IPB98(y,2) 的正规口径 P = P_tot = P_aux + P_alpha，
    且 W = tau_E P_tot = 3 n T V 反解 T。因 tau_E ∝ P_tot^-0.69 ⇒ T ∝ P_tot^0.31，
    P_alpha ∝ T^2 ∝ P_tot^0.62 < 1，故不动点稳定（自稳负反馈，不是热失控）。
  * 点火：P_alpha = P_loss = W/tau_E ⇒ n tau_E = 60 T_J/(<sigma v> E_fus)
    ⇒ n T tau_E > 3.10e21 keV s /m^3（在 T^2 模型内与 T 无关；真实曲线在 ~15 keV 取极小）。
  * 工程增益（教学简化，只计加热功率）：Q_eng = eta_th eta_aux (Q+1)，eta_th = eta_aux = 0.4；
    由此 Q_eng = 1 需要 Q = 1/(eta_th eta_aux) - 1 = 5.25。
"""
import math

e = 1.602176634e-19
u = 1.66053906660e-27
m_e = 9.1093837015e-31
m_D = 2.013553212745 * u
mu0 = 4.0 * math.pi * 1e-7

E_FUS = 17.6e6 * e            # J，D-T 每次反应释放的能量
ETA_TH = 0.4                  # 热-电转换效率（教学取值）
ETA_AUX = 0.4                 # 加热系统插头效率（教学取值）

# 名称, I[MA], B[T], n19, P[MW], R[m], a[m], kappa, M, q0, q95(典型/设计)
DEV = [
    ('D1 ITER 量级', 15.0, 5.3, 10.1, 150.0, 6.20, 2.00, 1.7, 2.5, 1.05, 3.0),
    ('D2 JET 量级', 3.0, 3.0, 4.0, 15.0, 2.96, 1.25, 1.4, 2.0, 1.00, 3.9),
    ('D3 DIII-D 量级', 1.5, 2.1, 5.0, 10.0, 1.67, 0.67, 1.8, 2.0, 0.85, 4.0),
]


def ipb98y2(I, B, n19, P, R, a, kap, M, H=1.0):
    eps = a / R
    return (0.0562 * H * I ** 0.93 * B ** 0.15 * n19 ** 0.41 * P ** -0.69
            * R ** 1.97 * eps ** 0.58 * kap ** 0.78 * M ** 0.19)


def volume(R, a, kap):
    return 2.0 * math.pi ** 2 * R * a ** 2 * kap


def sigma_dt(T_keV):
    """第 10 章卡片的 D-T 反应率近似，返回 m^3/s（T 单位 keV）。"""
    return 1.1e-24 * T_keV ** 2


def t_from_W(W_J, n, V):
    """由总储能反解温度，返回 keV。"""
    return W_J / (3.0 * n * V) / 1000.0 / e


def pfus(n, V, T_keV):
    """总聚变功率 [MW]。"""
    rate = n ** 2 / 4.0 * sigma_dt(T_keV)          # 1/(m^3 s)
    return rate * E_FUS * V / 1e6


def point(I, B, n19, P, R, a, kap, M):
    """返回 (tau_E, V, T[keV], W[MJ], chi_perp[m^2/s])——与附录 S/T 完全同口径。"""
    tau = ipb98y2(I, B, n19, P, R, a, kap, M)
    V = volume(R, a, kap)
    n = n19 * 1e19
    T = t_from_W(tau * P * 1e6, n, V)
    W = 3.0 * n * T * 1000.0 * e * V
    chi = 3.0 * a ** 2 / (8.0 * tau)
    return tau, V, T, W / 1e6, chi


print('=== 常数 ===')
print('e = %.9e C, m_D = %.9e kg, mu0 = %.9e H/m' % (e, m_D, mu0))
print('E_fus = 17.6 MeV = %.6e J, eta_th = %.2f, eta_aux = %.2f'
      % (E_FUS, ETA_TH, ETA_AUX))

print()
print('=== U.1 三个装置在同一外部加热功率下的增益 ===')
print('（P_aux 取附录 S 的工作点功率；tau_E、T、W 与附录 S 表 S.2 完全一致）')
print('%-16s %8s %8s %8s %9s %10s %7s'
      % ('装置', 'P_aux/MW', 'tau_E/s', 'T/keV', 'W/MJ', 'P_fus/MW', 'Q'))
q1 = {}
for name, I, B, n19, P, R, a, kap, M, q0, q95 in DEV:
    tau, V, T, W, chi = point(I, B, n19, P, R, a, kap, M)
    n = n19 * 1e19
    Pf = pfus(n, V, T)
    Q = Pf / P
    q1[name] = dict(tau=tau, V=V, T=T, W=W, Pf=Pf, Q=Q, n=n, P=P)
    print('%-16s %8.1f %8.4f %8.3f %9.1f %10.1f %7.3f'
          % (name, P, tau, T, W, Pf, Q))

print()
print('=== U.1b alpha 功率与中子功率（D1 工作点） ===')
d = q1['D1 ITER 量级']
print('P_fus = %.1f MW, P_alpha = P_fus/5 = %.1f MW, P_neutron = 0.8 P_fus = %.1f MW'
      % (d['Pf'], d['Pf'] / 5.0, 0.8 * d['Pf']))
print('alpha 功率密度 = %.3e W/m^3（V = %.1f m^3）' % (d['Pf'] / 5.0 * 1e6 / d['V'], d['V']))

print()
print('=== U.2 加热份额 f_alpha = Q/(Q+5) = 1/(1+5/Q) ===')
print('%8s %12s %12s' % ('Q', 'f_alpha', 'P_alpha/P_aux'))
for Q in (1.0, 2.0, 3.710, 5.0, 10.0, 20.0, 50.0):
    print('%8.3f %12.4f %12.4f' % (Q, Q / (Q + 5.0), 1.0 / (1.0 + 5.0 / Q)))
print('%8s %12.4f %12s' % ('inf', 1.0, '-'))

print()
print('=== U.3 自加热闭环的不动点迭代 ===')
print('（P_tot = P_aux + P_alpha，tau_E = IPB98(P=P_tot)，W = tau_E P_tot = 3 n T V）')


def fixed_point(I, B, n19, P_aux, R, a, kap, M, P0=None, steps=12):
    """返回迭代轨迹 [(P_tot, tau_E, T, P_fus, P_alpha, Q)] 与收敛值。"""
    n = n19 * 1e19
    V = volume(R, a, kap)
    P_tot = P_aux if P0 is None else P0
    traj = []
    for _ in range(steps):
        tau = ipb98y2(I, B, n19, P_tot, R, a, kap, M)
        T = t_from_W(tau * P_tot * 1e6, n, V)
        Pf = pfus(n, V, T)
        Pa = Pf / 5.0
        traj.append((P_tot, tau, T, Pf, Pa, Pf / P_aux))
        P_tot = P_aux + Pa
    return traj


I1, B1, n1, P1, R1, a1, k1, M1 = 15.0, 5.3, 10.1, 150.0, 6.20, 2.00, 1.7, 2.5
traj = fixed_point(I1, B1, n1, P1, R1, a1, k1, M1)
print('%-4s %10s %9s %9s %10s %10s %8s'
      % ('步', 'P_tot/MW', 'tau_E/s', 'T/keV', 'P_fus/MW', 'P_alpha/MW', 'Q'))
for i, (Pt, tau, T, Pf, Pa, Q) in enumerate(traj, 1):
    print('%-4d %10.2f %9.4f %9.3f %10.2f %10.2f %8.3f' % (i, Pt, tau, T, Pf, Pa, Q))
burn = traj[-1]
print('收敛：P_tot = %.1f MW, tau_E = %.4f s, T = %.3f keV, P_fus = %.1f MW, '
      'P_alpha = %.1f MW, Q = %.3f'
      % (burn[0], burn[1], burn[2], burn[3], burn[4], burn[5]))
print('注：迭代为「代入 -> 求新 P_tot」，逐步收敛；相邻两步差 = %.3f MW'
      % abs(traj[-1][0] - traj[-2][0]))

print()
print('=== U.3b 指数为何是 0.31 / 0.62（自稳的根据） ===')
for Pa in (100.0, 150.0, 225.0):
    tau = ipb98y2(I1, B1, n1, Pa, R1, a1, k1, M1)
    V1 = volume(R1, a1, k1)
    T = t_from_W(tau * Pa * 1e6, n1 * 1e19, V1)
    print('P_tot=%6.1f MW -> tau_E=%7.4f s, T=%6.3f keV, T/P_tot^0.31=%8.5f'
          % (Pa, tau, T, T / Pa ** 0.31))
slope = math.log(traj[3][2] / traj[1][2]) / math.log(traj[3][0] / traj[1][0])
print('由迭代第 2、4 步反解温度指数 dlnT/dlnP_tot = %.4f（理论 0.31）' % slope)
print('=> P_alpha ∝ T^2 ∝ P_tot^%.3f < 1，故闭环增益小于 1，不动点稳定'
      % (2.0 * slope))

print()
print('=== U.4 点火：三重积阈值 ===')
print('由 P_alpha = P_loss = W/tau_E 得 n tau_E = 60 T_J/(<sigma v> E_fus)')
for Tk in (5.0, 8.0, 10.0, 15.0, 20.0):
    TJ = Tk * 1000.0 * e
    ntau = 60.0 * TJ / (sigma_dt(Tk) * E_FUS)
    print('T=%5.1f keV -> n tau_E = %.4e m^-3 s, n T tau_E = %.4e keV s /m^3'
          % (Tk, ntau, ntau * Tk))
Tk = d['T']
TJ = Tk * 1000.0 * e
ntau_ign = 60.0 * TJ / (sigma_dt(Tk) * E_FUS)
ntt_ign = ntau_ign * Tk
print('取 D1 温度 T = %.3f keV：n tau_E 阈值 = %.4e m^-3 s，n T tau_E 阈值 = %.4e'
      % (Tk, ntau_ign, ntt_ign))
print()
print('%-16s %8s %8s %8s %12s %10s' % ('装置', 'n19', 'T/keV', 'tau_E/s', 'n T tau_E', '占阈值'))
for name, I, B, n19, P, R, a, kap, M, q0, q95 in DEV:
    rec = q1[name]
    ntt = rec['n'] * rec['T'] * rec['tau']
    print('%-16s %8.1f %8.3f %8.4f %12.3e %9.1f%%'
          % (name, n19, rec['T'], rec['tau'], ntt, 100.0 * ntt / ntt_ign))
need = ntt_ign / (d['T'] * d['tau'])
print()
print('D1 若要在同一 T、tau_E 下点火：n 需 %.3f -> %.3f (x%.4f)'
      % (d['n'] / 1e19, need / 1e19, need / d['n']))
P_need = 5.0 * d['P'] * need / d['n']    # 需要的 P_alpha -> P_fus
print('等价地说：P_alpha 要从 %.1f MW 提到 P_loss = %.1f MW，即 P_fus 从 %.1f 提到 %.1f MW'
      % (d['Pf'] / 5.0, d['P'], d['Pf'], 5.0 * d['P']))
print('用 P_fus ∝ n^2 的粗估：n x sqrt(%.1f/%.1f) = x%.4f -> n19 = %.3f'
      % (5.0 * d['P'], d['Pf'], math.sqrt(5.0 * d['P'] / d['Pf']),
         n1 * math.sqrt(5.0 * d['P'] / d['Pf'])))
n_G = I1 / (math.pi * a1 ** 2)           # 1e20 m^-3
print('Greenwald: n_G = I_p/(pi a^2) = %.4f x 1e20 m^-3' % n_G)
print('D1 工作点 n/n_G = %.1f%%；点火所需 n/n_G = %.1f%%'
      % (100.0 * n1 / 1e1 / n_G, 100.0 * n1 * math.sqrt(5.0 * d['P'] / d['Pf']) / n_G))

print()
print('=== U.5 Q 随外部加热功率下降（固定 D1 的 n、B、I_p） ===')
print('%-10s %9s %9s %9s %10s %8s' % ('P_aux/MW', 'tau_E/s', 'T/keV', 'W/MJ', 'P_fus/MW', 'Q'))
V1 = volume(R1, a1, k1)
nD1 = n1 * 1e19
for Pa in (150.0, 100.0, 50.0, 20.0, 10.0, 5.0):
    tau = ipb98y2(I1, B1, n1, Pa, R1, a1, k1, M1)
    T = t_from_W(tau * Pa * 1e6, nD1, V1)
    Pf = pfus(nD1, V1, T)
    print('%-10.1f %9.4f %9.3f %9.1f %10.1f %8.3f'
          % (Pa, tau, T, 3.0 * nD1 * T * 1000.0 * e * V1 / 1e6, Pf, Pf / Pa))
lo, hi = 1.0, 150.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    tau = ipb98y2(I1, B1, n1, mid, R1, a1, k1, M1)
    T = t_from_W(tau * mid * 1e6, nD1, V1)
    if pfus(nD1, V1, T) / mid > 10.0:
        lo = mid
    else:
        hi = mid
tau10 = ipb98y2(I1, B1, n1, lo, R1, a1, k1, M1)
T10 = t_from_W(tau10 * lo * 1e6, nD1, V1)
Pf10 = pfus(nD1, V1, T10)
print('Q = 10 的解：P_aux = %.3f MW, tau_E = %.4f s, T = %.3f keV, P_fus = %.1f MW'
      % (lo, tau10, T10, Pf10))
print('注：T < 5 keV 时 T^2 近似高估 <sigma v> 约 2 倍，故低温行的绝对值不可当真，只看趋势。')

print()
print('=== U.5b 工程增益 Q_eng = eta_th eta_aux (Q+1)（教学简化） ===')
for Q in (3.710, 6.072, 10.0, 5.25):
    print('Q = %6.3f -> Q_eng = %.3f' % (Q, ETA_TH * ETA_AUX * (Q + 1.0)))
print('Q_eng = 1 的阈值：Q = 1/(eta_th eta_aux) - 1 = %.3f'
      % (1.0 / (ETA_TH * ETA_AUX) - 1.0))
Qb = burn[5]
Pf_it, Pa_it = 500.0, 50.0
print('ITER 设计点（P_fus = %.0f MW, P_aux = %.0f MW, Q = %.0f）：' % (Pf_it, Pa_it, Pf_it / Pa_it))
print('  热功率 P_th = %.0f MW -> P_e = %.0f MW；插头功率 P_grid = %.0f MW；Q_eng = %.3f'
      % (Pf_it + Pa_it, ETA_TH * (Pf_it + Pa_it), Pa_it / ETA_AUX,
         ETA_TH * ETA_AUX * (Pf_it / Pa_it + 1.0)))

print()
print('=== U.6 练习数值 ===')
print('练习 1：Q = 3.710 -> f_alpha = %.4f；Q = 10 -> f_alpha = %.4f'
      % (d['Q'] / (d['Q'] + 5.0), 10.0 / 15.0))
print('练习 2：P_neutron = %.1f MW；第一壁取 800 m^2 -> 平均中子壁负载 %.4f MW/m^2'
      % (0.8 * d['Pf'], 0.8 * d['Pf'] / 800.0))
print('练习 3：三步迭代 P_tot = %.1f / %.1f / %.1f MW，Q = %.3f / %.3f / %.3f'
      % (traj[0][0], traj[1][0], traj[2][0], traj[0][5], traj[1][5], traj[2][5]))
print('练习 4：n T tau_E 阈值 = %.4e，D1 实有 %.4e（%.1f%%）'
      % (ntt_ign, d['n'] * d['T'] * d['tau'], 100.0 * d['n'] * d['T'] * d['tau'] / ntt_ign))
print('练习 5：固定 n、B、I_p 时 Q = 10 只需 P_aux = %.2f MW，但 P_fus 仅 %.1f MW（T = %.2f keV）'
      % (lo, Pf10, T10))
print('练习 6：Q = 3.710 -> Q_eng = %.3f（<1）；Q = 10 -> Q_eng = %.3f（>1）；阈值 Q = %.2f'
      % (ETA_TH * ETA_AUX * (d['Q'] + 1.0), ETA_TH * ETA_AUX * 11.0,
         1.0 / (ETA_TH * ETA_AUX) - 1.0))

print()
print('=== U.5c 三个代表点的完整工程账（教学简化：只计加热功率） ===')
print('%-14s %9s %9s %9s %9s %10s %8s'
      % ('工作点', 'P_fus/MW', 'P_aux/MW', 'P_th/MW', 'P_e/MW', 'P_grid/MW', 'Q_eng'))
for name, Pf, Pa in (('D1 外部加热点', d['Pf'], d['P']),
                     ('自加热闭环点', burn[3], d['P']),
                     ('ITER 设计点', 500.0, 50.0)):
    Pth = Pf + Pa
    print('%-14s %9.1f %9.1f %9.1f %9.1f %10.1f %8.3f'
          % (name, Pf, Pa, Pth, ETA_TH * Pth, Pa / ETA_AUX,
             ETA_TH * ETA_AUX * (Pf / Pa + 1.0)))

print()
print('以上数字是附录 U 正文与门禁 P23 的唯一来源。')
