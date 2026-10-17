# -*- coding: utf-8 -*-
"""附录 T（磁面、磁岛与磁通量管）与门禁 P22 的数值唯一来源。

用法：python tmp/appT_calc.py > tmp/appT_calc.txt

口径（正文必须与之一致）：
  * 装置沿用附录 S 的三个自洽工作点（同一组输入），本附录只增加 q 剖面与拓扑量：
      q(rho) = q0 + (q95 - q0) rho^2 ，rho = r/a          （教学剖面，q 面上凸）
      q95 取"装置典型/设计值"；同时用第 10 章 10.8 节的简化公式
      q95 = 5 a^2 B (1+kappa^2)/(2 R I_p) 复核，差值单独列出。
  * 有理面：q(r)=m/n；极向磁场 B_theta(r) = B r/(q R)；剪切 s = r q'/q，q' = dq/dr。
  * 岛宽（柱位形、单螺旋扰动，标准结果）：w = 4 sqrt( Psi_tilde /(B_theta s) )
      Psi_tilde = 螺旋磁通扰动（Wb/rad）。w 为岛的全宽。
  * 平行热导（Braginskii，Z=1）：chi_par = 3.16 T tau_e / m_e；电子碰撞时间
      tau_e = 3.44e5 T_eV^{3/2} /(n_cm3 lnLambda)，lnLambda = 24 - ln(sqrt(n_cm3)/T_eV)（NRL，T_e>10 eV）
  * 岛内等温时间 tau_flat = w^2/chi_par；横越磁面时间尺度 tau_perp = w^2/chi_perp（chi_perp 取附录 S）
  * 岛重叠（Chirikov 量级判据）：相邻有理面间距 dr 被两个岛链的半宽之和吃掉 ->
      w/2 + w/2 > dr，即 w > dr；据此反解所需 Psi_tilde。
  * 修正 Rutherford 方程的自举项（教学模型）：Delta'_bs = 5.6 mu0 j_bs L_q/(B_theta w)，
      临界岛宽 w_marg = 5.6 mu0 j_bs L_q/(B_theta |Delta'_0|)，j_bs = f_bs I_p/(pi a^2 kappa)，L_q = q/q'
      （|Delta'_0| 与 f_bs 是教学取值，做敏感性表而不给单一"答案"）
  * 锯齿：q0 < 1 时反转半径 r_inv = a sqrt((1-q0)/(q95-q0))，否则无 q=1 面
"""
import math

e = 1.602176634e-19
u = 1.66053906660e-27
m_e = 9.1093837015e-31
m_D = 2.013553212745 * u
mu0 = 4.0 * math.pi * 1e-7

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


def point(I, B, n19, P, R, a, kap, M):
    """返回 (tau_E, V, T[keV], W[MJ], chi_perp[m^2/s])——与附录 S 完全同口径。"""
    tau = ipb98y2(I, B, n19, P, R, a, kap, M)
    V = 2.0 * math.pi ** 2 * R * a ** 2 * kap
    n = n19 * 1e19
    T = tau * P * 1e6 / (3.0 * n * V)          # J
    W = 3.0 * n * T * V                        # J
    chi = 3.0 * a ** 2 / (8.0 * tau)
    return tau, V, T / 1000.0 / e, W / 1e6, chi


print('=== 常数 ===')
print('e = %.9e C, m_e = %.9e kg, m_D = %.9e kg, mu0 = %.9e H/m'
      % (e, m_e, m_D, mu0))

print()
print('=== T.1 三个装置的 q 剖面参数与 q95 复核 ===')
print('%-16s %5s %6s %8s %9s %9s %8s' % ('装置', 'q0', 'q95设计', 'q95公式', '差值', 'a/m', 'B/T'))
qrows = []
for name, I, B, n19, P, R, a, kap, M, q0, q95 in DEV:
    q95f = 5.0 * a ** 2 * B * (1.0 + kap ** 2) / (2.0 * R * I)
    print('%-16s %5.2f %6.2f %8.2f %9.2f %9.2f %8.2f'
          % (name, q0, q95, q95f, q95 - q95f, a, B))
    qrows.append((name, I, B, n19, P, R, a, kap, M, q0, q95, q95f))

print()
print('=== T.2 有理面半径与剪切 ===')
print('%-16s %7s %8s %8s %8s %8s %8s %8s %8s'
      % ('装置', 'q面', 'rho', 'r/m', 'B_theta', 'q\'', 's', 'L_q/m', 'k_par*'))
qsurf = {}
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    qsurf[name] = {}
    for qt in (1.0, 1.5, 2.0, 3.0):
        if qt <= q0 or qt >= q95:
            print('%-16s %7.2f   (该面不在 q0<q<q95 区间内)' % (name, qt))
            continue
        rho = math.sqrt((qt - q0) / (q95 - q0))
        r = rho * a
        qp = 2.0 * (q95 - q0) * r / a ** 2
        s = r * qp / qt
        Bt = B * r / (qt * R)
        Lq = qt / qp
        kpar = (1.0 - 1.0) / (qt * R)   # 在有理面上 k_par = 0（m/n = q）
        print('%-16s %7.2f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f'
              % (name, qt, rho, r, Bt, qp, s, Lq, kpar))
        qsurf[name][qt] = dict(rho=rho, r=r, qp=qp, s=s, Bt=Bt, Lq=Lq)

print()
print('=== T.2b 非共振处的平行波数 k_par = (m - n q)/(q R)（D1, R=6.20）===')
R1 = 6.20
for (m, n) in ((2, 1), (3, 2)):
    for frac in (0.9, 0.95, 1.0, 1.05, 1.1):
        qt = m / n * frac
        kpar = (m - n * qt) / (qt * R1)
        print('  m/n=%d/%d  q=%.3f  k_par=%+9.3e 1/m   平行波长 2pi/k_par=%s m'
              % (m, n, qt, kpar, ('%.1f' % (2 * math.pi / kpar)) if kpar else 'inf'))

print()
print('=== T.3 平行与垂直输运：岛内为什么等温 ===')
print('%-16s %7s %7s %8s %10s %10s %11s %11s'
      % ('装置', 'T/keV', 'n19', 'lnLam', 'tau_e/s', 'chi_par', 'chi_perp', 'chi_par/chi_p'))
trans = {}
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    tau, V, TkeV, W, chi = point(I, B, n19, P, R, a, kap, M)
    n_cm3 = n19 * 1e19 / 1e6
    lnL = 24.0 - math.log(math.sqrt(n_cm3) / (TkeV * 1000.0))
    tau_e = 3.44e5 * (TkeV * 1000.0) ** 1.5 / (n_cm3 * lnL)
    Tj = TkeV * 1000.0 * e
    chi_par = 3.16 * Tj * tau_e / m_e
    print('%-16s %7.3f %7.1f %8.2f %10.3e %10.3e %11.4f %11.3e'
          % (name, TkeV, n19, lnL, tau_e, chi_par, chi, chi_par / chi))
    trans[name] = dict(tau=tau, T=TkeV, chi=chi, chi_par=chi_par, tau_e=tau_e,
                       lnL=lnL, V=V, W=W)

print()
print('=== T.3b 同一个 1 cm 管子的两条路（D1）===')
nm, I, B, n19, P, R, a, kap, M, q0, q95, q95f = qrows[0]
cp = trans[nm]['chi_par']
cperp = trans[nm]['chi']
L_par = 2.0 * math.pi * R
t_par = L_par ** 2 / cp
t_perp = (0.01) ** 2 / cperp
print('  绕环一周 L_par = %.2f m；沿磁力线扩散时间 L^2/chi_par = %.3e s' % (L_par, t_par))
print('  横越 1 cm：  L^2/chi_perp = %.3e s' % t_perp)
print('  比值 t_perp/t_par = %.3e  -> 平行方向快 %d 个量级'
      % (t_perp / t_par, round(math.log10(t_perp / t_par))))

print()
print('=== T.4 磁岛宽度 w = 4 sqrt(Psi/(B_theta s))，取 Psi = 1e-4 Wb/rad ===')
print('%-16s %7s %9s %9s %9s %10s %9s'
      % ('装置', 'q面', 'r/m', 'B_theta/T', 's', 'w/m', 'w/a'))
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    for qt in (1.5, 2.0):
        d = qsurf[name].get(qt)
        if not d:
            continue
        Psi = 1e-4
        w = 4.0 * math.sqrt(Psi / (d['Bt'] * d['s']))
        print('%-16s %7.2f %9.4f %9.4f %9.4f %10.4f %9.4f'
              % (name, qt, d['r'], d['Bt'], d['s'], w, w / a))
        qsurf[name].setdefault('w%d' % int(qt * 10), w)

print()
print('=== T.4b 岛内等温时间 tau_flat = w^2/chi_par 与 tau_E 之比 ===')
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    w = qsurf[name]['w15']
    cp = trans[name]['chi_par']
    tf = w ** 2 / cp
    print('%-16s w(1.5)=%8.4f m  tau_flat=%10.3e s  tau_E=%8.4f s  比值=%.3e'
          % (name, w, tf, trans[name]['tau'], tf / trans[name]['tau']))

print()
print('=== T.5 岛链重叠：需要多大的螺旋扰动 ===')
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    d15 = qsurf[name].get(1.5)
    d20 = qsurf[name].get(2.0)
    if not (d15 and d20):
        print('%-16s 缺 q=1.5 或 q=2 面，跳过' % name)
        continue
    dr = d20['r'] - d15['r']
    Psi_overlap = d20['Bt'] * d20['s'] * (dr / 4.0) ** 2
    print('%-16s r(1.5)=%7.4f  r(2)=%7.4f  dr=%7.4f m' % (name, d15['r'], d20['r'], dr))
    print('%-16s 重叠要求 w > dr -> Psi > %.4e Wb/rad  (= %6.1f 倍教学值 1e-4)'
          % ('', Psi_overlap, Psi_overlap / 1e-4))
    print('%-16s 教学岛宽 w(1e-4)=%7.4f m（半宽和 %7.4f m vs 间距 %7.4f m）'
          % ('', d15['s'] and qsurf[name]['w15'], qsurf[name]['w15'], dr))

print()
print('=== T.6 新经典撕裂模：临界岛宽 w_marg = 5.6 mu0 j_bs L_q/(B_theta |Delta0|) ===')
f_bs = 0.25
print('教学取值：自举电流份额 f_bs = %.2f；|Delta0| 做敏感性表' % f_bs)
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    d = qsurf[name].get(2.0)
    A = math.pi * a ** 2 * kap
    j_phi = I * 1e6 / A
    j_bs = f_bs * j_phi
    print('%-16s j_phi=%.3e A/m^2  j_bs=%.3e A/m^2  L_q=%.3f m  B_theta=%.4f T'
          % (name, j_phi, j_bs, d['Lq'], d['Bt']))
    for d0 in (20.0, 50.0, 100.0):
        wm = 5.6 * mu0 * j_bs * d['Lq'] / (d['Bt'] * d0)
        print('%-16s   |Delta0|=%6.1f 1/m -> w_marg=%8.4f m = %6.2f cm (w/a=%.3f)'
              % ('', d0, wm, wm * 100, wm / a))

print()
print('=== T.7 锯齿：q0<1 时的反转半径 ===')
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    if q0 < 1.0:
        rho = math.sqrt((1.0 - q0) / (q95 - q0))
        print('%-16s q0=%.2f < 1 -> r_inv = %.4f m (rho=%.4f)，混合区体积份额 ~ %.4f'
              % (name, q0, rho * a, rho, rho ** 2))
    else:
        print('%-16s q0=%.2f >= 1 -> 无 q=1 面，不出现 m/n=1/1 锯齿（需靠 q0 剖面控制）'
              % (name, q0))

print()
print('=== T.4c q=2 面的岛内等温时间（表 T.4 第二、四、六行）===')
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    w = qsurf[name]['w20']
    tf = w ** 2 / trans[name]['chi_par']
    print('%-16s w(2.0)=%8.4f m  tau_flat=%10.3e s  tau_flat/tau_E=%.3e'
          % (name, w, tf, tf / trans[name]['tau']))

print()
print('=== T.2c 剪切恒等式核对 s = 2(q - q0)/q ===')
for (name, I, B, n19, P, R, a, kap, M, q0, q95, q95f) in qrows:
    for qt in (1.0, 1.5, 2.0, 3.0):
        d = qsurf[name].get(qt)
        if not d:
            continue
        print('%-16s q=%.2f  s(表)=%.4f  s(恒等式)=%.4f  B_theta/B=%.4f'
              % (name, qt, d['s'], 2.0 * (qt - q0) / qt, d['Bt'] / B))

print()
print('=== T.8 练习答案 ===')
nm, I, B, n19, P, R, a, kap, M, q0, q95, q95f = qrows[0]
s15 = qsurf[nm][1.5]
s20 = qsurf[nm][2.0]
w15 = qsurf[nm]['w15']
w20 = qsurf[nm]['w20']
print('① D1 的 q=1.5 面：r = %.4f m，rho = %.4f' % (s15['r'], s15['rho']))
print('② q=1.5 面的剪切 s = %.4f；B_theta = %.4f T' % (s15['s'], s15['Bt']))
Psi2 = 4e-4
w2 = 4.0 * math.sqrt(Psi2 / (s15['Bt'] * s15['s']))
print('③ Psi 由 1e-4 增到 4e-4：w = %.4f m（是 1e-4 时的 %.2f 倍；w ∝ sqrt(Psi)）'
      % (w2, w2 / w15))
print('④ 岛内等温时间 tau_flat = %.3e s（chi_par=%.3e m^2/s）'
      % (w15 ** 2 / trans[nm]['chi_par'], trans[nm]['chi_par']))
print('⑤ 重叠所需 Psi = %.4e Wb/rad（D1 的 2/1 与 3/2 链，dr=%.4f m）'
      % (s20['Bt'] * s20['s'] * ((s20['r'] - s15['r']) / 4.0) ** 2, s20['r'] - s15['r']))
print('⑥ D1 的 q=2 面教学岛宽 w = %.4f m（w/a=%.4f），比 q=1.5 面略窄'
      % (w20, w20 / a))
print()
print('以上数字是附录 T 正文与门禁 P22 的唯一来源。')
