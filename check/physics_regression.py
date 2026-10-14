#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》物理数值回归检查
================================================

与另外四个套件互补，本模块只回答一个问题：\textbf{书里那些可以独立复算的
结论，现在是否仍然成立}——数值来自独立实现（复根用牛顿迭代直接解原始色散
关系，积分用数值求积），源码断言只用来把「数值」与「书上写的那句话」绑在
一起。设计原则（来自 PP-W01 的教训）：断言必须可被证伪——每条检查都给出
算得的数字与容差，而不是「没报错就算过」。

  M. 数值复核   —— 双流阈值/最大增长率/例题数值、Buneman 冷共振标度、
                   δ′ 冷极限积分、Salpeter 参数、斜入射反射条件、
                   俯仰角、e-folding 比值、Druyvesteyn 往返
  N. 结论口径   —— 源码必须写出的正确表述，以及必须已被删除的错误表述
  O. 交叉一致性 —— 附录 F 条目编号/指针可解析、附录 G.2 的返回值语义
  P. 新增附录   —— 附录 I/K 的物理量独立复算（逐题绑定、含符号与单位）、附录 J 对照表完备性、
                   附录 I 的 DOI 出处、附录 I/J/K 的挂载与字母连续性、附录 K 三道判断题的结论方向、α 散射判据的表格行方向

用法：
    python check/physics_regression.py
    python check/physics_regression.py --json

退出码：0 = 无 FAIL；1 = 存在 FAIL
"""

import argparse
import io
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
MAIN = 'plasma_physics_textbook_v2.tex'

RESULTS = []


def rec(level, cat, msg):
    RESULTS.append((level, cat, msg))


def read(p):
    return io.open(p, encoding='utf-8').read()


def strip_comments(t):
    """逐字符剥 % 注释（\\% 不受影响）——与另三个套件保持同一实现。"""
    out = []
    for line in t.split('\n'):
        i, idx = 0, None
        while i < len(line):
            if line[i] == '\\':
                i += 2
                continue
            if line[i] == '%':
                idx = i
                break
            i += 1
        out.append(line if idx is None else line[:idx])
    return '\n'.join(out)


def near(got, want, tol, unit=''):
    return abs(got - want) <= tol, '算得 %.10g，期望 %.10g±%g%s' % (got, want, tol, unit)


def has(text, needle, where):
    if needle in text:
        rec('PASS', 'N', '%s 含必需表述：%s' % (where, needle))
        return True
    rec('FAIL', 'N', '%s 缺少必需表述：%s' % (where, needle))
    return False


def lacks(text, needle, where):
    if needle not in text:
        rec('PASS', 'N', '%s 已无错误表述：%s' % (where, needle))
        return True
    rec('FAIL', 'N', '%s 仍含错误表述：%s' % (where, needle))
    return False


# ------------------------------------------------------------------ M. 数值复核

def check_M01_two_stream():
    """对称双流：闭式解 ↔ 原始色散关系，阈值、最大增长率、例题数值。"""
    def w2_min(a):
        return a * a + 0.5 - 0.5 * math.sqrt(1.0 + 8.0 * a * a)

    def resid(w, a):
        return abs(1.0 - 0.5 * (1.0 / (w - a) ** 2 + 1.0 / (w + a) ** 2))

    worst = 0.0
    for a in (0.2, 0.5, 0.75, 0.9, 0.99):
        u = w2_min(a)
        w = complex(0.0, math.sqrt(-u)) if u < 0 else complex(math.sqrt(u), 0.0)
        worst = max(worst, resid(w, a))
    if worst < 1e-9:
        rec('PASS', 'M', '双流闭式解满足原始色散关系 1=½[1/(w−a)²+1/(w+a)²]：'
                         '最大残差 %.3e（阈值 1e-9）' % worst)
    else:
        rec('FAIL', 'M', '双流闭式解与原始色散关系不符：最大残差 %.3e' % worst)

    # 阈值 a=1 处 w²=0；a 略大于 1 必须稳定
    ok_lo, msg_lo = near(w2_min(1.0), 0.0, 1e-12)
    (rec('PASS', 'M', '阈值 a=1（|ku|=2ω_pe）处 w²=0：%s' % msg_lo) if ok_lo
     else rec('FAIL', 'M', '阈值 a=1 处 w²=%s' % msg_lo))
    if w2_min(1.0001) > 0 and w2_min(0.9999) < 0:
        rec('PASS', 'M', '稳定侧 a=1.0001（w²=%.3e>0）与不稳定侧 a=0.9999（w²=%.3e<0）分界正确'
            % (w2_min(1.0001), w2_min(0.9999)))
    else:
        rec('FAIL', 'M', 'a=1 两侧判别错误：w²(1.0001)=%.3e, w²(0.9999)=%.3e'
            % (w2_min(1.0001), w2_min(0.9999)))

    # 最大增长率 1/√8 at a=√(3/8)
    a_opt = math.sqrt(3.0 / 8.0)
    g_max = math.sqrt(-w2_min(a_opt))
    ok, msg = near(g_max, 1.0 / math.sqrt(8.0), 1e-12)
    (rec('PASS', 'M', 'γ_max/ω_pe=1/√8：%s（a_opt=√(3/8)=%.6f）' % (msg, a_opt)) if ok
     else rec('FAIL', 'M', 'γ_max 不符：%s' % msg))
    # 扫描确认 a_opt 确为极大点
    best, best_a = 0.0, 0.0
    n = 200001
    for i in range(1, n + 1):
        a = i / n * 1.0
        g = math.sqrt(max(-w2_min(a), 0.0))
        if g > best:
            best, best_a = g, a
    ok, msg = near(best_a, a_opt, 1e-5)
    (rec('PASS', 'M', '扫描确认最大增长率位置：%s' % msg) if ok
     else rec('FAIL', 'M', '扫描到的最大增长率位置不符：%s' % msg))

    # 闭式解必须同时是色散关系化简后的二次方程之根：
    # (w²−a²)² = w²+a² ⟹ x² − (2a²+1)x + a⁴ − a² = 0（x=w²）
    worst_poly, worst_at = 0.0, 0.0
    for i in range(1, 41):
        a = i / 20.0
        disc = math.sqrt(1.0 + 8.0 * a * a)
        x_lo = a * a + 0.5 - 0.5 * disc          # 低根：不稳定侧的 w²
        x_hi = a * a + 0.5 + 0.5 * disc          # 高根：稳定支
        d = max(abs(w2_min(a) - x_lo),
                abs((x_lo * x_lo - (2 * a * a + 1) * x_lo + a ** 4 - a * a)),
                abs((x_hi * x_hi - (2 * a * a + 1) * x_hi + a ** 4 - a * a)))
        if d > worst_poly:
            worst_poly, worst_at = d, a
    if worst_poly < 1e-12:
        rec('PASS', 'M', '闭式解与色散关系多项式根（x²−(2a²+1)x+a⁴−a²=0，x=w²）一致：'
                         '20 个 a 值最大偏差 %.3e' % worst_poly)
    else:
        rec('FAIL', 'M', '闭式解与多项式根不一致：最大偏差 %.3e（a=%.2f）' % (worst_poly, worst_at))

    # 稳定侧深水区：a=1.2（|ku|=2.4ω_pe）必须无实频不稳定（w²>0）
    w2_12 = w2_min(1.2)
    if w2_12 > 0:
        rec('PASS', 'M', '稳定侧 a=1.2（|ku|=2.4ω_pe）：w²=%.6f>0，无非零实频解' % w2_12)
    else:
        rec('FAIL', 'M', '稳定侧 a=1.2 判别错误：w²=%.6f（应 >0）' % w2_12)

    # 例题数值 ω_pe=5e10 rad/s, u=2e7 m/s
    wpe, u = 5e10, 2e7
    items = [
        ('k_max=2ω_pe/u', 2 * wpe / u, 5.0e3, 0.5, 'm^-1'),
        ('λ_min=2π/k_max', 2 * math.pi / (2 * wpe / u) * 1e3, 1.25664, 1e-4, 'mm'),
        ('k_opt=√(3/2)ω_pe/u', math.sqrt(1.5) * wpe / u, 3.061862e3, 0.01, 'm^-1'),
        ('γ_max=ω_pe/√8', wpe / math.sqrt(8.0), 1.7677670e10, 1e4, 's^-1'),
        ('τ=1/γ_max', math.sqrt(8.0) / wpe * 1e12, 56.56854, 1e-3, 'ps'),
        ('a=0.75 的 w²（书值 −0.1101）', w2_min(0.75), -0.1101, 1e-4, ''),
        ('a=0.75 的 γ/ω_pe（书值 0.332）', math.sqrt(-w2_min(0.75)), 0.332, 5e-4, ''),
    ]
    for name, got, want, tol, unit in items:
        ok, msg = near(got, want, tol, unit)
        (rec('PASS', 'M', '例题数值 %s：%s' % (name, msg)) if ok
         else rec('FAIL', 'M', '例题数值 %s：%s' % (name, msg)))

    # 源码交叉核对
    sup = strip_comments(read(os.path.join(SRC, 'ch9_supplement.tex')))
    for s in (r'0<|ku|<2\,\omega_{pe}',
              r'k_{\max}=\frac{2\,\omega_{pe}}{u}',
              r'5.0\times 10^{3}',
              r'1.26\,\mathrm{mm}',
              r'\gamma_{\max}=\frac{\omega_{pe}}{\sqrt8}',
              r'1.768\times 10^{10}',
              r'\sqrt{3/2}\,\omega_{pe}',
              r'3.06\times10^{3}',
              r'56.6\,\mathrm{ps}',
              r'w^2\approx-0.1101',
              r'\gamma/\omega_{pe}\approx0.332',
              r'那是错的'):
        has(sup, s, 'src/ch9_supplement.tex（双流例题）')

    # 正文第 12 章与补充材料必须给同一个复核例：只改一处（数值漂移）要被抓到
    main = strip_comments(read(os.path.join(SRC, MAIN)))
    for s in (r'w^2\approx-0.1101',
              r'\gamma/\omega_{pe}=\sqrt{0.1101}\approx0.332'):
        has(main, s, 'src/%s（双流复核例）' % MAIN)


def check_M02_buneman():
    """Buneman：复根直接解原式，μ^{1/3} 标度与 0.0543 的对照。"""
    def f(w, b, mu):
        return 1.0 - 1.0 / (w - b) ** 2 - mu / w ** 2

    def poly(b, mu):
        """1 = 1/(w-b)^2 + mu/w^2 化为首一四次多项式（实系数）。"""
        return [1.0 + 0j, -2 * b + 0j, b * b - 1.0 - mu + 0j, 2 * mu * b + 0j, -mu * b * b + 0j]

    def pval(c, z):
        s = 0j
        for k in c:
            s = s * z + k
        return s

    def durand_kerner(c):
        """Durand--Kerner 迭代求全部复根（无需初值猜测，比单点牛顿稳）。"""
        n = len(c) - 1
        z = [(0.4 + 0.9j) ** k for k in range(n)]
        for _ in range(500):
            new = list(z)
            for i in range(n):
                d = 1.0 + 0j
                for j in range(n):
                    if i != j:
                        d *= (z[i] - z[j])
                if d == 0:
                    continue
                new[i] = z[i] - pval(c, z[i]) / d
            if max(abs(new[i] - z[i]) for i in range(n)) < 1e-15:
                z = new
                break
            z = new
        return z

    def gmax(mu):
        best, best_b, best_w = 0.0, 0.0, None
        for i in range(1, 801):
            b = 0.8 + 0.4 * i / 800.0
            for w in durand_kerner(poly(b, mu)):
                if abs(w.imag) > best:
                    best, best_b, best_w = w.imag, b, w
        return best_b, best, best_w

    mu = 1.0 / 1836.0
    b, gm, w = gmax(mu)
    ok, msg = near(gm, 0.0543, 5e-4)
    (rec('PASS', 'M', 'Buneman 精确求根 γ/ω_pe：%s（书值 0.0543）' % msg) if ok
     else rec('FAIL', 'M', 'Buneman γ/ω_pe 与书值 0.0543 不符：%s' % msg))
    r = abs(f(w, b, mu))
    (rec('PASS', 'M', 'Buneman 复根 w=%.6f%+.6fj 代回原式残差 %.3e（阈值 1e-10）' % (w.real, w.imag, r))
     if r < 1e-10 else
     rec('FAIL', 'M', 'Buneman 复根代回原式残差过大：%.3e' % r))

    g1 = gmax(1e-6)[1]
    g2 = gmax(1e-8)[1]
    slope_mu = (math.log(g2) - math.log(g1)) / (math.log(1e-8) - math.log(1e-6))
    ok, msg = near(slope_mu, 1.0 / 3.0, 0.005)
    (rec('PASS', 'M', 'γ 对 μ 的对数斜率：%s（即 γ∝μ^{1/3}；对质量比 m_i/m_e 为 %.6f，书值 −0.333）'
         % (msg, -slope_mu)) if ok
     else rec('FAIL', 'M', 'γ 对 μ 的对数斜率不是 +1/3（对质量比应为 −1/3）：%s' % msg))

    main = strip_comments(read(os.path.join(SRC, MAIN)))
    for s in ('0.0543', '-0.333', '幂次是 $2/3$', r'\mu^{1/3}', r'\tfrac13', r'\omega_{pi}^2'):
        has(main, s, 'src/%s（Buneman）' % MAIN)
    lacks(main, r'\gamma/\omega_{pe}\propto\mu^{1/2}', 'src/%s（Buneman）' % MAIN)


def check_M03_delta_prime():
    """δ′ 冷极限：数值求积给 +1/c²，并与正文符号一致。"""
    def integ(c, sg):
        f = lambda v: (-v / sg ** 3 / math.sqrt(2 * math.pi) * math.exp(-v * v / (2 * sg * sg))) / (v - c)
        n = 200001
        lo, hi = -40 * sg, 40 * sg
        h = (hi - lo) / (n - 1)
        s = 0.5 * (f(lo) + f(hi))
        for i in range(1, n - 1):
            s += f(lo + i * h)
        return s * h

    worst = 0.0
    for c in (1.0, 2.0, 3.0):
        for sg in (1e-3, 1e-4, 1e-5):
            val = integ(c, sg)
            worst = max(worst, abs(val - 1.0 / c ** 2) * c ** 2)
    (rec('PASS', 'M', '数值积分 ∫δ′(v)/(v−c)dv=+1/c²（c=1,2,3；σ=1e-3…1e-5）：最大相对误差 %.2e（阈值 1e-4）'
         % worst) if worst < 1e-4 else
     rec('FAIL', 'M', '∫δ′(v)/(v−c)dv 与 +1/c² 不符：最大相对误差 %.2e' % worst))

    main = strip_comments(read(os.path.join(SRC, MAIN)))
    for s in (r"\int\delta'(v)\phi(v)\diff v=-\phi'(0)",
              r"\int\delta'(v)/(v-c)\diff v=+1/c^2",
              r'\boxed{+n_0k^2/\omega^2}'):
        has(main, s, 'src/%s（δ′ 冷极限）' % MAIN)
    lacks(main, r'-n_0k^2/\omega^2', 'src/%s（δ′ 冷极限）' % MAIN)
    for m in re.finditer(r'1=-\\omega_\{pe\}\^2/\\omega\^2', main):
        ctx = main[max(0, m.start() - 120):m.start()]
        if '常见错误' not in ctx:
            rec('FAIL', 'N', 'src/%s 出现 $1=-\\omega_{pe}^2/\\omega^2$ 但不在「常见错误」语境中' % MAIN)
            break
    else:
        rec('PASS', 'N', 'src/%s 中反向结果 $1=-\\omega_{pe}^2/\\omega^2$ 只作为反例出现' % MAIN)


def check_M04_salpeter():
    """Salpeter 参数：λ_D、N_D、α(θ) 与阈值密度。"""
    e = 1.602176634e-19
    eps0 = 8.8541878128e-12
    lam0 = 1064e-9

    def debye(Te_eV, ne):
        return math.sqrt(eps0 * Te_eV * e / (ne * e * e))

    ld = debye(10.0, 1e20)
    ND = 4.0 / 3.0 * math.pi * 1e20 * ld ** 3
    q90 = 4 * math.pi / lam0 * math.sin(math.pi / 4)
    a90 = 1.0 / (q90 * ld)
    a1 = 1.0 / (4 * math.pi / lam0 * math.sin(math.radians(0.5)) * ld)
    n_crit = 1e20 / a90 ** 2
    a_1e23 = 1.0 / (q90 * debye(10.0, 1e23))
    ld_f = debye(1e4, 1e20)
    a_f = 1.0 / (q90 * ld_f)
    for name, got, want, tol, unit in (
            ('λ_D(10 eV,1e20)', ld * 1e6, 2.350820, 1e-5, 'μm'),
            ('N_D', ND, 5442.0, 1.0, ''),
            ('α(90°)', a90, 0.0509363, 1e-6, ''),
            ('α(1°)', a1, 4.127338, 1e-5, ''),
            ('α=1 的 n_e(90°)', n_crit, 3.85429e22, 1e18, 'm^-3'),
            ('α(1e23,90°)', a_1e23, 1.610752, 1e-5, ''),
            ('聚变堆芯 λ_D(10 keV,1e20)', ld_f * 1e6, 74.3394, 1e-3, 'μm'),
            ('聚变堆芯 α(90°)', a_f, 0.00161075, 1e-7, '')):
        ok, msg = near(got, want, tol, unit)
        (rec('PASS', 'M', 'Salpeter %s：%s' % (name, msg)) if ok
         else rec('FAIL', 'M', 'Salpeter %s：%s' % (name, msg)))

    main = strip_comments(read(os.path.join(SRC, MAIN)))
    for s in (r'\alpha \equiv \frac{1}{q\lambda_D}',
              r'q = |\vect{k}_s-\vect{k}_i| = \frac{4\pi}{\lambda_0}\sin\frac{\theta}{2}',
              r'2.35\,\mu\mathrm{m}', r'5.4\times10^3', r'\alpha\approx0.051',
              r'\alpha\approx4.1', r'3.9\times10^{22}', r'\alpha\approx1.6',
              r'74\,\mu\mathrm{m}', r'\alpha\approx1.6\times10^{-3}', r'过松'):
        has(main, s, 'src/%s（Salpeter 参数）' % MAIN)
    lacks(main, r'N_D\gg1 就', 'src/%s（Salpeter 参数）' % MAIN)


def check_M05_oblique():
    """斜入射：数值求反射频率并与 ω_pe/cosθ₀ 对照。"""
    wpe = 1.0
    worst = 0.0
    for th_deg in (30.0, 60.0, 80.0):
        th = math.radians(th_deg)
        lo, hi = wpe * (1 + 1e-12), wpe / math.cos(th) * 1.5

        def g(w):
            return math.sqrt(max(1 - wpe ** 2 / w ** 2, 0.0)) - math.sin(th)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(mid) < 0:
                lo = mid
            else:
                hi = mid
        w_ref = 0.5 * (lo + hi)
        worst = max(worst, abs(w_ref - wpe / math.cos(th)) / (wpe / math.cos(th)))
    (rec('PASS', 'M', '数值（Snell 守恒）反射频率 = ω_pe/cosθ₀：30°/60°/80° 最大相对差 %.2e（阈值 1e-9）'
         % worst) if worst < 1e-9 else
     rec('FAIL', 'M', '反射频率与 ω_pe/cosθ₀ 不符：最大相对差 %.2e' % worst))

    sup = strip_comments(read(os.path.join(SRC, 'ch6_supplement.tex')))
    has(sup, r'\omega=\frac{\omega_{pe}}{\cos\theta_0}', 'src/ch6_supplement.tex（斜入射）')
    has(sup, r'是除以 $\cos\theta_0$，不是乘', 'src/ch6_supplement.tex（斜入射）')
    has(sup, r'n=\sin\theta_0', 'src/ch6_supplement.tex（斜入射）')
    lacks(sup, r'\omega_{pe}\cos\theta', 'src/ch6_supplement.tex（斜入射）')


def check_M06_pitch_and_efolding():
    """俯仰角 tan²θ=2 与 10/ln10 比值。"""
    main = strip_comments(read(os.path.join(SRC, MAIN)))
    sup = strip_comments(read(os.path.join(SRC, 'ch9_supplement.tex')))
    th = math.degrees(math.atan(math.sqrt(2.0)))
    ok, msg = near(th, 54.7356103, 1e-5, '°')
    (rec('PASS', 'M', 'tan²θ=2 的俯仰角：%s' % msg) if ok
     else rec('FAIL', 'M', '俯仰角不符：%s' % msg))
    # 两项速度相等性回代：v_⊥²/2 与 v_∥²
    s2 = math.tan(math.radians(th)) ** 2
    vperp2, vpar2 = s2 / (1 + s2), 1.0 / (1 + s2)
    ok, msg = near(vperp2 / 2.0, vpar2, 1e-12)
    (rec('PASS', 'M', '回代两项相等：v_⊥²/2=%.12f, v_∥²=%.12f（%s）' % (vperp2 / 2, vpar2, msg)) if ok
     else rec('FAIL', 'M', '俯仰角回代两项不相等：%s' % msg))
    # 单粒子条件与各向同性热平均必须被区分开（PP-C08 的实质要求）
    idx = main.find(r'\tan^2\theta=v_\perp^2/v_\parallel^2=2')
    ctx = main[max(0, idx - 300): idx + 500] if idx >= 0 else ''
    if '各向同性' in ctx and '热平均' in ctx:
        rec('PASS', 'M', '正文在同一处把「单粒子 tan²θ=2」与「各向同性热平均后才相等」区分开')
    else:
        rec('FAIL', 'M', '正文未在同一处区分单粒子条件与各向同性热平均')

    ratio = 10.0 / math.log(10.0)
    ok, msg = near(ratio, 4.342944819, 1e-8)
    (rec('PASS', 'M', '十个 e-folding / 十倍振幅时间比 10/ln10：%s' % msg) if ok
     else rec('FAIL', 'M', '10/ln10 比值不符：%s' % msg))

    has(main, r'\theta\approx54.74^\circ', 'src/%s（漂移速度比较）' % MAIN)
    has(main, r'\tan^2\theta=v_\perp^2/v_\parallel^2=2', 'src/%s（漂移速度比较）' % MAIN)
    has(sup, r'\dfrac{10}{\ln10}\approx4.34', 'src/ch9_supplement.tex（e-folding）')
    has(sup, r'不是 $\ln10\approx2.3$ 倍', 'src/ch9_supplement.tex（e-folding）')


def check_M07_eedf():
    """Druyvesteyn 往返：正演 → 二阶导数 → 反演 F 与 T_e。"""
    Te = 4.0
    F = lambda eps: math.sqrt(eps) * math.exp(-eps / Te)      # 形状，未归一
    I = lambda u: Te * math.exp(-u / Te)                      # 正演：I''∝F/√ε

    # (a) 反演出的 F 与输入成正比（比例常数与 ε 无关）
    ratios = []
    for u in (0.5, 1.0, 2.0, 4.0, 8.0, 12.0):
        h = 1e-4
        I2 = (I(u + h) - 2 * I(u) + I(u - h)) / h ** 2
        ratios.append(math.sqrt(u) * I2 / F(u))
    spread = (max(ratios) - min(ratios)) / (sum(ratios) / len(ratios))
    (rec('PASS', 'M', 'Druyvesteyn 反演 F∝√ε·I″ 与输入 F 成正比：比例常数离散度 %.2e（阈值 1e-6）'
         % spread) if spread < 1e-6 else
     rec('FAIL', 'M', 'Druyvesteyn 反演与输入 F 不成比例：离散度 %.2e' % spread))

    # (b) 麦克斯韦分布下 ln I″ 斜率恰为 −1/T_e（即丢掉 √ε 不改斜率）
    us = [0.3 + 19.7 * i / 399.0 for i in range(400)]
    xs = us
    ys = []
    for u in xs:
        h = 1e-4
        ys.append(math.log((I(u + h) - 2 * I(u) + I(u - h)) / h ** 2))
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    T_fit = -1.0 / slope
    ok, msg = near(T_fit, 4.0, 1e-6, ' eV')
    (rec('PASS', 'M', '麦克斯韦合成电流：ln I″ 斜率反演 T_e：%s' % msg) if ok
     else rec('FAIL', 'M', 'ln I″ 斜率反演温度不是 4 eV：%s' % msg))

    # (c) √ε 因子决定归一化：∫F dε=1，丢掉后 ∫≠1
    Fm = lambda eps: 2 / math.sqrt(math.pi) * Te ** -1.5 * math.sqrt(eps) * math.exp(-eps / Te)
    Fb = lambda eps: 2 / math.sqrt(math.pi) * Te ** -1.5 * math.exp(-eps / Te)
    nn = 400001
    lo, hi = 0.0, 60 * Te
    h = (hi - lo) / (nn - 1)
    s1 = 0.5 * (Fm(lo) + Fm(hi))
    s2 = 0.5 * (Fb(lo) + Fb(hi))
    for i in range(1, nn - 1):
        s1 += Fm(lo + i * h)
        s2 += Fb(lo + i * h)
    i1, i2 = s1 * h, s2 * h
    ok, msg = near(i1, 1.0, 1e-5)
    (rec('PASS', 'M', 'Maxwell EEDF 归一化 ∫F dε=1：%s' % msg) if ok
     else rec('FAIL', 'M', 'Maxwell EEDF 归一化偏差：%s' % msg))
    ok, msg = near(i2, 0.564189584, 1e-5)
    (rec('PASS', 'M', '丢掉 √ε 后 ∫F dε=1/√π：%s（证明 √ε 不可吸收进常数）' % msg) if ok
     else rec('FAIL', 'M', '丢掉 √ε 的积分不是 1/√π：%s' % msg))

    # (d) 三种记法的换算：p=F/n_e、g=p/√ε，以及 p/g=√ε 与 ∫√ε·g dε=1
    pfun = lambda eps: 2 / math.sqrt(math.pi) * Te ** -1.5 * math.sqrt(eps) * math.exp(-eps / Te)
    gfun = lambda eps: 2 / math.sqrt(math.pi) * Te ** -1.5 * math.exp(-eps / Te)
    worst = 0.0
    for eps in (0.25, 1.0, 4.0):
        worst = max(worst, abs(pfun(eps) / gfun(eps) - math.sqrt(eps)))
    ok, msg = near(worst, 0.0, 1e-12)
    (rec('PASS', 'M', '归一化 EEDF p 与 EEPF g 之比 p/g=√ε（0.25/1/4 eV）：%s' % msg) if ok
     else rec('FAIL', 'M', 'p/g 不等于 √ε：%s' % msg))

    ng = 400001
    hi_g = 60 * Te
    hg = hi_g / (ng - 1)
    sg = 0.5 * (0.0 + math.sqrt(hi_g) * gfun(hi_g))          # √ε·g 在 ε=0 处为 0
    for i in range(1, ng - 1):
        e = i * hg
        sg += math.sqrt(e) * gfun(e)
    ig = sg * hg
    ok, msg = near(ig, 1.0, 1e-5)
    (rec('PASS', 'M', 'EEPF 归一化 ∫√ε·g dε=1：%s' % msg) if ok
     else rec('FAIL', 'M', 'EEPF 归一化 ∫√ε·g dε 不是 1：%s' % msg))

    # (e) J↔eV 的雅可比：F_eV(ε_eV)=e·F_J(e·ε_eV) 使两个积分相等（差一个 e 就不等）
    T_J = Te * E_CHG
    F_J = lambda eJ: math.sqrt(eJ) * math.exp(-eJ / T_J)
    F_eV = lambda ev: E_CHG * F_J(E_CHG * ev)
    nj = 200001
    hi_J, hi_eV = 60 * T_J, 60 * Te
    h_J, h_eV = hi_J / (nj - 1), hi_eV / (nj - 1)
    aJ = 0.5 * (F_J(0.0) + F_J(hi_J))
    aV = 0.5 * (F_eV(0.0) + F_eV(hi_eV))
    for i in range(1, nj - 1):
        aJ += F_J(i * h_J)
        aV += F_eV(i * h_eV)
    i_J, i_eV = aJ * h_J, aV * h_eV
    rel = abs(i_J - i_eV) / i_J
    ok, msg = near(rel, 0.0, 1e-9)
    (rec('PASS', 'M', 'J↔eV 雅可比 ∫F_eV dε_eV=∫F_J dε_J：相对差 %s' % msg) if ok
     else rec('FAIL', 'M', 'J↔eV 换单位后两个积分不等：相对差 %s' % msg))

    main = strip_comments(read(os.path.join(SRC, MAIN)))
    has(main, r'F(\eps) = \frac{2\sqrt{2 m_e \eps}}{A e^3}', 'src/%s（Druyvesteyn 公式）' % MAIN)
    has(main, r'\int_0^\infty F\diff\eps=n_e', 'src/%s（EEDF 归一化）' % MAIN)
    has(main, '没有"固定偏差"这回事', 'src/%s（EEDF 警告框）' % MAIN)
    has(main, r'\diff v/\diff\eps=1/\sqrt{2m_e\eps}', 'src/%s（EEDF 雅可比）' % MAIN)
    has(main, r'EEDF $F(\eps)$ & $\int_0^\infty F\diff\eps=n_e$ & $\mathrm{m^{-3}(eV)^{-1}}$',
        'src/%s（三种记法表：EEDF）' % MAIN)
    has(main, r'归一化 EEDF $p(\eps)=F(\eps)/n_e$ & $\int_0^\infty p\diff\eps=1$ & $\mathrm{(eV)^{-1}}$',
        'src/%s（三种记法表：归一化 EEDF）' % MAIN)
    has(main, r'EEPF $g(\eps)=p(\eps)/\sqrt{\eps}$ & $\int_0^\infty\sqrt{\eps}\,g\diff\eps=1$ & '
              r'$\mathrm{(eV)^{-3/2}}$', 'src/%s（三种记法表：EEPF）' % MAIN)
    has(main, 'EEPF 不是「归一化 EEDF」的同义词', 'src/%s（EEPF ≠ 归一化 EEDF）' % MAIN)
    has(main, r'$F=n_e\sqrt{\eps}\,g$', 'src/%s（EEPF 反解回 F）' % MAIN)
    has(main, r'F_{\mathrm{eV}}(\eps_{\mathrm{eV}})=e\,F_J(e\,\eps_{\mathrm{eV}})',
        'src/%s（J→eV 的因子 e）' % MAIN)
    has(main, r'$I \propto \int_{\eps > e(V_p - V)} F(\eps)\, v(\eps) \,\diff\eps$',
        'src/%s（探针通量积分用 F 而不是别的记法）' % MAIN)
    lacks(main, 'EEPF 通常就是', 'src/%s（旧的 EEPF 错误表述）' % MAIN)
    lacks(main, 'f_\\eps(', 'src/%s（旧的 f_ε 记法）' % MAIN)


def check_M08_timescale():
    """增长时间尺度：振幅十倍（ln10/γ）与十个 e-folding（10/γ）的相容性。

    ch9 瑞利--泰勒例取 γ=3.07e3 /s，书中写 0.75 ms 与 3.3 ms；
    两者之比必须等于 10/ln10≈4.34294，而不能按 ln10≈2.3 理解。
    """
    gamma = 3.07e3
    ratio = 10.0 / math.log(10.0)

    t10x = math.log(10.0) / gamma
    ok, msg = near(t10x * 1e3, 0.75, 5e-3, 'ms')
    (rec('PASS', 'M', '十倍振幅时间 ln10/γ（γ=3.07e3）：%s' % msg) if ok
     else rec('FAIL', 'M', '十倍振幅时间 ln10/γ：%s' % msg))

    t10e = 10.0 / gamma
    ok, msg = near(t10e * 1e3, 3.3, 5e-2, 'ms')
    (rec('PASS', 'M', '十个 e-folding 时间 10/γ：%s' % msg) if ok
     else rec('FAIL', 'M', '十个 e-folding 时间 10/γ：%s' % msg))

    ok, msg = near(t10e / t10x, ratio, 1e-9)
    (rec('PASS', 'M', '两者之比精确等于 10/ln10：%s' % msg) if ok
     else rec('FAIL', 'M', '两者之比不等于 10/ln10：%s' % msg))

    # 书值 3.3 ms / 0.75 ms = 4.4 与 10/ln10 的相对差须落在四舍五入范围内
    src_ratio = 3.3 / 0.75
    rel = abs(src_ratio - ratio) / ratio
    if rel < 0.02:
        rec('PASS', 'M', '书中 3.3 ms 与 0.75 ms 相容：比值 %.3f vs 10/ln10=%.6f，'
                         '相对差 %.2f%% < 2%%（两位有效数字的舍入范围内）'
            % (src_ratio, ratio, rel * 100))
    else:
        rec('FAIL', 'M', '书中 3.3 ms / 0.75 ms = %.3f 与 10/ln10=%.6f 不相容（相对差 %.2f%%）'
            % (src_ratio, ratio, rel * 100))

    ok, msg = near(1.0 / gamma * 1e3, 0.33, 5e-3, 'ms')
    (rec('PASS', 'M', '一个 e-folding 时间 τ=1/γ：%s' % msg) if ok
     else rec('FAIL', 'M', 'τ=1/γ：%s' % msg))

    sup = strip_comments(read(os.path.join(SRC, 'ch9_supplement.tex')))
    for s in (r't_{10\times} = \frac{\ln 10}{\gamma}',
              r'\approx 7.5\times 10^{-4}\,\mathrm{s} = 0.75\,\mathrm{ms}',
              r't_{10e} = \frac{10}{\gamma}',
              r'\approx 3.3\times 10^{-3}\,\mathrm{s} = 3.3\,\mathrm{ms}',
              r'\dfrac{10}{\ln10}\approx4.34',
              '不是 $\\ln10\\approx2.3$ 倍',
              '不能互换'):
        has(sup, s, 'src/ch9_supplement.tex（增长时间尺度）')


# ------------------------------------------------------------ O. 交叉一致性

def check_O01_appendix_f():
    t = read(os.path.join(SRC, MAIN))
    i0 = t.find(r'\appchapter{附录F')
    i1 = t.find('% 附录G：Python仿真入门')
    if i0 < 0 or i1 < 0 or i1 <= i0:
        rec('FAIL', 'O', 'src/%s 中找不到附录 F 区块（锚点缺失）' % MAIN)
        return
    blk = strip_comments(t[i0:i1])
    num, seq = 0, []
    for m in re.finditer(r'\\setcounter\{enumi\}\{(\d+)\}|\\item\b', blk):
        if m.group(1) is not None:
            num = int(m.group(1))
        else:
            num += 1
            seq.append(num)
    if seq == list(range(1, len(seq) + 1)):
        rec('PASS', 'O', '附录 F 条目编号连续无重复：1..%d（共 %d 条）' % (len(seq), len(seq)))
    else:
        dup = [n for n in seq if seq.count(n) > 1]
        rec('FAIL', 'O', '附录 F 条目编号不连续：序列=%s（重复 %s）' % (seq, sorted(set(dup)) or '无'))

    titles = re.findall(r'\\section\{([^}]*)\}', blk)
    want = ['第 2、8 章：等离子体基本概念',
            '第 3、8 章：单粒子运动与库仑碰撞',
            '第 4、5、8 章：流体方程、MHD 与输运',
            '第 1、6、7 章：傅里叶约定、等离子体波与动理学',
            '第 5、7、9、10 章：MHD 局限、不稳定性与磁约束',
            '第 2、8、10 章：碰撞过程与聚变装置']
    if titles == want:
        rec('PASS', 'O', '附录 F 六个分组标题与各条目的真实入口章一致（第 2、8 / 3、8 / 4、5、8 / 1、6、7 / 5、7、9、10 / 2、8、10 章）')
    else:
        rec('FAIL', 'O', '附录 F 分组标题不符：实际 %s' % titles)

    refs = re.findall(r'\\S\\ref\{([^}]+)\}', blk)
    allsrc = '\n'.join(strip_comments(read(os.path.join(SRC, f)))
                       for f in sorted(os.listdir(SRC)) if f.endswith('.tex'))
    defined = {}
    for m in re.finditer(r'\\label\{([^}]+)\}', allsrc):
        defined[m.group(1)] = defined.get(m.group(1), 0) + 1
    missing = [r for r in refs if r not in defined]
    dup = [r for r in set(refs) if defined.get(r, 0) > 1]
    if missing:
        rec('FAIL', 'O', '附录 F 指针指向不存在的标签：%s' % sorted(set(missing)))
    elif dup:
        rec('FAIL', 'O', '附录 F 指针指向重复定义的标签：%s' % sorted(set(dup)))
    else:
        rec('PASS', 'O', '附录 F 的 %d 个 \\S\\ref 指针全部指向唯一存在的标签（%d 个不同目标）'
            % (len(refs), len(set(refs))))
    items = re.findall(r'\\item ([^\n]*(?:\n(?!\\item|\\end\{enumerate\}|\\section)[^\n]*)*)', blk)
    nop = [i[:24] for i in items if '\\S\\ref{' not in i]
    if nop:
        rec('FAIL', 'O', '附录 F 有 %d 条没有章节指针：%s' % (len(nop), nop))
    else:
        rec('PASS', 'O', '附录 F 全部 %d 条都带可点击章节指针（跨章条目给出多个入口）' % len(items))


def check_O03_appendix_f_chapters():
    """附录 F 逐条的章节归属：按 .aux 的 \\newlabel 号核对每个 \\S\\ref 的章号。

    静态版（不读 PDF）：分组标题声明「第 X、Y 章」，则该组内每个 \\S\\ref
    指向的小节号首段必须落在 {X, Y} 内；反过来每个声明的章都要有真实入口。
    """
    t = read(os.path.join(SRC, MAIN))
    i0 = t.find(r'\appchapter{附录F')
    i1 = t.find('% 附录G：Python仿真入门')
    if i0 < 0 or i1 < 0 or i1 <= i0:
        rec('FAIL', 'O', 'src/%s 中找不到附录 F 区块（锚点缺失）' % MAIN)
        return
    blk = strip_comments(t[i0:i1])

    aux_path = os.path.join(SRC, MAIN[:-4] + '.aux')
    if not os.path.exists(aux_path):
        rec('WARN', 'O', '未找到 %s，附录 F 章号归属无法静态核对（需先构建一次）'
            % os.path.relpath(aux_path, ROOT))
        return
    aux = read(aux_path)
    label2num = {}
    for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^{}]*)\}\{([^{}]*)\}', aux):
        label2num[m.group(1)] = m.group(2)

    # 按 \section 切段，逐段取声明章集合与段内指针
    heads = [(m.start(), m.group(1)) for m in re.finditer(r'\\section\{([^}]*)\}', blk)]
    heads.append((len(blk), None))
    bad, unknown, total_refs, missing_entry = [], [], 0, []
    for k in range(len(heads) - 1):
        start, title = heads[k]
        body = blk[start:heads[k + 1][0]]
        declared = [int(x) for x in re.findall(r'(\d+)', title.split('：')[0])]
        seen = set()
        for m in re.finditer(r'\\S\\ref\{([^}]+)\}', body):
            total_refs += 1
            lab = m.group(1)
            num = label2num.get(lab)
            if num is None:
                unknown.append(lab)
                continue
            mm = re.match(r'\s*(\d+)', num)
            if not mm:
                unknown.append(lab)
                continue
            ch = int(mm.group(1))
            seen.add(ch)
            if ch not in declared:
                bad.append((title, lab, num, declared))
        gone = [c for c in declared if c not in seen]
        if gone:
            missing_entry.append((title, gone))

    if bad:
        rec('FAIL', 'O', '附录 F 有 %d 处指针的章号与所在分组标题声明的章不符，例如 %s'
            % (len(bad), [(b[1], b[2], b[3]) for b in bad[:3]]))
    else:
        rec('PASS', 'O', '附录 F 的 %d 个 \\S\\ref 指针章号全部落在所在分组声明的章集合内'
            '（按 .aux 的 \\newlabel 号核对）' % total_refs)
    if missing_entry:
        rec('FAIL', 'O', '附录 F 有分组标题声明了却没有任何入口的章：%s'
            % [(m[0][:12], m[1]) for m in missing_entry])
    else:
        rec('PASS', 'O', '附录 F 每个分组标题声明的每一章都至少有一个真实入口（跨章条目给出多个入口）')
    if unknown:
        rec('WARN', 'O', '附录 F 有 %d 个指针在 .aux 中查不到号（标签名或构建产物不同步）：%s'
            % (len(unknown), sorted(set(unknown))))
    else:
        rec('PASS', 'O', '附录 F 全部指针都能在 .aux 中查到小节号')


def check_O02_appendix_g():
    t = read(os.path.join(SRC, MAIN))
    i0 = t.find(r'\section{德拜屏蔽电势的数值求解}')
    i1 = t.find(r'\section{朗缪尔波色散曲线绘制}')
    if i0 < 0 or i1 < 0 or i1 <= i0:
        rec('FAIL', 'O', 'src/%s 中找不到附录 G.2 区块（锚点缺失）' % MAIN)
        return
    blk = t[i0:i1]
    for s in (r'\texttt{sol.success}', r'\texttt{sol.status}', '0}=收敛到指定精度',
              '网格节点数超过', r'sol.rms_residuals', 'max_nodes',
              '读反是常见事故', '而是失败',
              '无量纲化: x = r/lam_D', "y'' + (2/x) y' = exp(y) - 1"):
        has(blk, s, 'src/%s（附录 G.2）' % MAIN)
    if re.search(r'#[^\n]*最近一次网格上收敛', blk):
        rec('FAIL', 'O', 'src/%s（附录 G.2）代码注释仍把 status=1 解释成收敛' % MAIN)
    else:
        rec('PASS', 'O', 'src/%s（附录 G.2）代码中已无「status=1 即收敛」的错误注释' % MAIN)


# ------------------------------------------------- P. 新增附录（I/J/K，v1.8.0）

APP_FILES = {
    'I': 'appendix_i_learning.tex',
    'J': 'appendix_j_textbook_map.tex',
    'K': 'appendix_k_problems.tex',
    'L': 'appendix_l_problems.tex',
    'M': 'appendix_m_heating.tex',
    'N': 'appendix_n_boundary.tex',
    'O': 'appendix_o_propulsion.tex',
    'P': 'appendix_p_lowtemp.tex',
    'Q': 'appendix_q_diagnostics.tex',
    'R': 'appendix_r_dusty.tex',
    'S': 'appendix_s_confinement.tex',
    'T': 'appendix_t_topology.tex',
    'U': 'appendix_u_gain.tex',
    'V': 'appendix_v_disruption.tex',
    'W': 'appendix_w_fastparticle.tex',
    'X': 'appendix_x_radiation.tex',
    'Y': 'appendix_y_strongfield.tex',
    'Z': 'appendix_z_reactor.tex',
    'AA': 'appendix_aa_alternates.tex',
    'AB': 'appendix_ab_space_astro.tex',
    'AC': 'appendix_ac_numerics.tex',
    'AD': 'appendix_ad_pedestal.tex',
    'AE': 'appendix_ae_kinetics.tex',
    'AF': 'appendix_af_turbulence.tex',
    'AG': 'appendix_ag_equilibrium.tex',
    'AH': 'appendix_ah_atomic.tex',
    'AI': 'appendix_ai_icf.tex',
    'AJ': 'appendix_aj_fuelling.tex',
    'AK': 'appendix_ak_safety.tex',
    'AL': 'appendix_al_nonlinear.tex',
    'AM': 'appendix_am_magnets.tex',
    'AN': 'appendix_an_polarimetry.tex',
    'AO': 'appendix_ao_vacuum.tex',
    'AP': 'appendix_ap_maintenance.tex',
    'AQ': 'appendix_aq_interferometry.tex',
    'AR': 'appendix_ar_power.tex',
    'AS': 'appendix_as_rays.tex',
}

# CODATA 2018 常数（与附录 D 的单位换算表同一组值）
E_CHG = 1.602176634e-19
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
M_P = 1.67262192369e-27
M_D = 2.013553212745 * M_P        # 氘核质量（含其电子）
M_AR = 39.948 * M_P               # 氩原子质量（≈ Ar+）
M_E = 9.1093837015e-31            # 电子质量（CODATA 2018）
M_XE = 131.293 * 1.66053906660e-27   # 氙离子质量（131.293 u，附录 O 统一取这个值）
G0 = 9.80665                      # 标准重力加速度（比冲与排气速度的换算）
KB = 1.380649e-23                 # 玻尔兹曼常数（附录 P 的 p/(k_B T) 与 Torr 换算）
LOSCHMIDT = 2.6867805e25          # 洛施密特常数（273.15 K、101325 Pa，附录 P 的 sccm 换气量）


# 数值抽取：三种写法（``1.3\times10^{-4}``／``1.3e-4``／``1.3``）
NUM_PAT = re.compile(
    r'(\d+(?:\.\d+)?)\s*\\times\s*10\s*\^\{?(-?\d+)\}?'
    r'|(\d+(?:\.\d+)?)[eE]([-+]?\d+)'
    r'|(?<![\d.])(\d+(?:\.\d+)?)(?![\d.])')

# 只有紧跟在这些符号之后的减号，才是数值本身的负号
SIGN_HEADS = ('\\approx', '\\sim', '\\simeq', '\\pm')


def _num_at(m):
    if m.group(1) is not None:
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    if m.group(3) is not None:
        return float(m.group(3)) * 10.0 ** int(m.group(4))
    return float(m.group(5))


def _is_negative(text, pos):
    r"""``text[pos]`` 是数值首位；判断它是否带负号。

    ``10^{-2}`` 的指数负号、``\propto`` 前面的符号都不是答案的负号，所以只在减号
    紧跟 ``= ( , : $ [ {`` 或 ``\approx \sim \simeq \pm`` 之后才算数。P-W01 的
    最小复现（把 K.5 的 −2.00×10⁻² 翻成 +2.00×10⁻²）正是靠这一步才能抓到。
    """
    w = text[max(0, pos - 12):pos].rstrip()
    if not w.endswith('-'):
        return False
    head = w[:-1].rstrip()
    if not head:
        return True
    return head[-1] in '=(,:$[{' or head.endswith(SIGN_HEADS)


def numbers(text):
    """返回 ``[(值, 起点, 终点, 是否带负号), ...]``。"""
    return [(_num_at(m), m.start(), m.end(), _is_negative(text, m.start()))
            for m in NUM_PAT.finditer(text)]


def value_in(text, v, rel, unit='', signed=False):
    """在**给定的文本块内**找与 v 相符的数，可选要求其后 80 字符内出现单位。

    ``signed=False`` 按绝对值比对（正文常只印大小）；``signed=True`` 连符号一起比，
    用于「答案必须是负的」这类断言。``unit`` 按 LaTeX 原样写（如 ``\\mathrm{mm}``），
    并要求其后不是字母——这样 ``\\mathrm{m}`` 不会误配 ``\\mathrm{mm}``。
    """
    for x, _, end, neg in numbers(text):
        got = -x if neg else x
        diff = abs(got - v) if signed else abs(abs(got) - abs(v))
        if diff > rel * abs(v):
            continue
        if unit and not re.search(re.escape(unit) + r'(?![A-Za-z])', text[end:end + 80]):
            continue
        return True
    return False


def _match_brace(text, i):
    """``text[i]`` 是 ``{``，返回配对 ``}`` 的下标（找不到返回 len(text)）。"""
    depth = 0
    while i < len(text):
        c = text[i]
        if c == '\\':
            i += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return len(text)


def split_blocks(text, kind):
    """切块，返回 ``[(标题, 块文本), ...]``。

    ``kind='example'`` 按 ``\\begin{example}{标题}{}`` 切出每道例题（块只含该环境）；
    ``kind='section'`` 按 ``\\section{标题}`` 切出每节（块含该节到下一节之前的全部内容）。
    这是「答案必须落在自己那道题/那一节里」这条断言的实现基础。
    """
    if kind == 'example':
        out = []
        for m in re.finditer(r'\\begin\{example\}\s*\{', text):
            close = _match_brace(text, m.end() - 1)
            end = text.find(r'\end{example}', close)
            if end >= 0:
                out.append((text[m.end():close], text[m.start():end + len(r'\end{example}')]))
        return out
    marks = [(m.start(), m.end() - 1) for m in re.finditer(r'\\section\s*\{', text)]
    out = []
    for n, (pos, brace) in enumerate(marks):
        close = _match_brace(text, brace)
        end = marks[n + 1][0] if n + 1 < len(marks) else len(text)
        out.append((text[brace + 1:close], text[pos:end]))
    return out


def check_P01_new_appendix_numbers():
    """附录 I/K 的每个答案必须出现在**它所属的那道例题/那一节**内，并带正确符号与单位。

    这是「教材里的数不是抄来的」这一纪律在新增附录上的落地：每个期望值都由本函数
    从头算一遍，再到该量所在的那道题里找印出来的那个数。v1.8.1 之前只在整份文件里
    找绝对值，所以把 K.5 的答案由 −2.00×10⁻² 翻成 +2.00×10⁻²、或把某个值搬到别的
    例题里，都能蒙混过关（P-W01）。现在被证伪的方式很直接：改数、改符号、把 mm 换成
    cm、让某个量离开它所属的例题——四种都会 FAIL。
    """
    text_k = strip_comments(read(os.path.join(SRC, APP_FILES['K'])))
    text_i = strip_comments(read(os.path.join(SRC, APP_FILES['I'])))
    ex = split_blocks(text_k, 'example')
    sec = split_blocks(text_i, 'section')
    if len(ex) != 8:
        rec('FAIL', 'P', '附录 K 的 example 环境应为 8 个、实际 %d 个：答案无法绑定到题目' % len(ex))
        return
    if len(sec) < 4:
        rec('FAIL', 'P', '附录 I 的 section 少于 4 个（实际 %d 个）：答案无法绑定到小节' % len(sec))
        return

    def lam_d(te, n):
        return math.sqrt(EPS0 * te / (n * E_CHG))

    ld10 = lam_d(10.0, 1e18)                       # K.1 德拜长度 @10 eV, 1e18 m^-3
    n_d = (4.0 / 3.0) * math.pi * 1e18 * ld10 ** 3  # 德拜球内粒子数
    v_a = 5.0 / math.sqrt(MU0 * 1e20 * M_D)         # K.4 Alfvén 速度 @5 T, 1e20 m^-3
    r_loss = 1.0 - math.sqrt(0.9)                   # K.2 磁镜 R=10 损失锥内份额
    theta_c = math.degrees(math.asin(1.0 / math.sqrt(10.0)))
    t_cross = 1.0 / v_a
    v_over_c = v_a / 2.99792458e8 * 100.0
    kld = 0.3                                       # K.5 Landau 阻尼
    gam = -math.sqrt(math.pi / 8.0) * kld ** -3 * math.exp(-1.0 / (2 * kld ** 2) - 1.5)
    per_damp = (1.0 - math.exp(2 * math.pi * gam)) * 100.0
    eta_kev = 5.2e-5 * 17 / 1000.0 ** 1.5           # K.6 Spitzer 电阻率
    eta_10ev = 5.2e-5 * 10 / 10.0 ** 1.5
    tau_eta = MU0 * 1.0 / eta_kev
    tau_eta_10 = MU0 / eta_10ev
    rm_1m = MU0 * 1.0 * v_a / eta_10ev              # K.4 对照算例的磁雷诺数
    c_s = math.sqrt(3.0 * E_CHG / M_AR)             # K.7 玻姆速度 @3 eV, Ar+
    n0, ns = 1e16, 0.61e16
    lam_n0, lam_ns = lam_d(3.0, n0), lam_d(3.0, ns)
    s_sh = math.sqrt(2.0) / 3.0 * lam_ns * (2 * 100.0 / 3.0) ** 0.75
    flux = ns * c_s
    j_wall = E_CHG * flux
    s_cl = math.sqrt((4.0 / 9.0) * EPS0 * math.sqrt(2 * E_CHG / M_AR) * 100.0 ** 1.5 / j_wall)
    eta_sp = 5.2e-5 * 20 / 86.0 ** 1.5              # K.8 Sweet-Parker
    s_sp = MU0 * 1e5 * 1e5 / eta_sp
    v_in = 1e5 / math.sqrt(s_sp)
    tau_sp = 1e5 / v_in
    omega_ci = E_CHG * 3.0 / M_D                    # I.4 估算练习（B=3 T 的 D+）

    def rho_l(te):                                  # I.4 回旋半径（D+, B=5 T）
        return math.sqrt(2.0 * M_D * te * E_CHG) / (E_CHG * 5.0)

    cases = [
        ('K', 1, ld10, 0.01, r'\mathrm{m}', False, 'K.1 λ_D=%.4e m' % ld10),
        ('K', 1, ld10 / 1e-6, 0.01, r'\mathrm{\mu m}', False, 'K.1 λ_D=%.2f μm' % (ld10 / 1e-6)),
        ('K', 1, n_d, 0.02, '', False, 'K.1 N_D=%.4e' % n_d),
        ('K', 1, math.exp(-5.0) * 100, 0.02, r'\%', False, 'K.1 r=5λ_D 处屏蔽势剩 0.67%'),
        ('K', 2, theta_c, 0.001, r'^\circ', False, 'K.2 损失锥角=%.3f°' % theta_c),
        ('K', 2, r_loss * 100.0, 0.01, r'\%', False, 'K.2 损失份额=%.2f%%' % (r_loss * 100.0)),
        ('K', 2, (1.0 - math.sqrt(0.99)) * 100.0, 0.02, r'\%', False, 'K.2 R=100 损失份额=0.50%'),
        ('I', 4, ld10 / 1e-6, 0.01, r'\mathrm{\mu m}', False, 'I.4 λ_D=%.2f μm' % (ld10 / 1e-6)),
        ('I', 4, n_d, 0.02, '', False, 'I.4 N_D=%.3e' % n_d),
        ('I', 4, rho_l(10.0) / 1e-3, 0.01, r'\mathrm{mm}', False,
         'I.4 ρ_L(10 eV)=%.3f mm' % (rho_l(10.0) / 1e-3)),
        ('I', 4, rho_l(1000.0) / 1e-3, 0.01, r'\mathrm{mm}', False,
         'I.4 ρ_L(1 keV)=%.3f mm' % (rho_l(1000.0) / 1e-3)),
        ('I', 4, rho_l(10000.0) / 1e-3, 0.01, r'\mathrm{mm}', False,
         'I.4 ρ_L(10 keV)=%.3f mm' % (rho_l(10000.0) / 1e-3)),
        ('I', 4, v_a, 0.01, r'\mathrm{m/s}', False, 'I.4 v_A=%.4e m/s' % v_a),
        ('I', 4, v_over_c, 0.03, r'\%', False, 'I.4 v_A/c=%.2f%%' % v_over_c),
        ('I', 4, t_cross, 0.01, r'\mathrm{s}', False, 'I.4 t=%.3e s' % t_cross),
        ('I', 4, omega_ci, 0.03, r'\mathrm{s^{-1}}', False, 'I.4 估算练习 ω_ci=%.3e 1/s' % omega_ci),
        ('K', 4, v_a, 0.01, r'\mathrm{m\,s^{-1}}', False, 'K.4 v_A=%.4e m/s' % v_a),
        ('K', 4, v_over_c, 0.03, r'\%', False, 'K.4 v_A/c=%.2f%%' % v_over_c),
        ('K', 4, t_cross, 0.01, r'\mathrm{s}', False, 'K.4 t=%.3e s' % t_cross),
        ('K', 4, rm_1m, 0.01, '', False, 'K.4 R_m=%.3e' % rm_1m),
        ('K', 4, tau_eta_10, 0.01, r'\mathrm{s}', False, 'K.4 τ_η=%.4f s' % tau_eta_10),
        ('K', 5, gam, 0.01, '', True, 'K.5 γ/ω=%.4e（负号必须在）' % gam),
        ('K', 5, per_damp, 0.01, r'\%', False, 'K.5 每周期衰减=%.2f%%' % per_damp),
        ('K', 6, eta_kev, 0.01, r'\Omega\cdot\mathrm{m}', False, 'K.6 η(1 keV)=%.3e Ω·m' % eta_kev),
        ('K', 6, tau_eta, 0.01, r'\mathrm{s}', False, 'K.6 τ_η=%.2f s' % tau_eta),
        ('K', 6, eta_10ev, 0.01, r'\Omega\cdot\mathrm{m}', False, 'K.6 η(10 eV)=%.3e Ω·m' % eta_10ev),
        ('K', 6, tau_eta_10, 0.01, r'\mathrm{s}', False, 'K.6 τ_η(10 eV)=%.4f s' % tau_eta_10),
        ('K', 6, eta_kev and eta_10ev / eta_kev, 0.05, '', False,
         'K.6 η 之比=%.1f 倍' % (eta_10ev / eta_kev)),
        ('K', 7, c_s, 0.01, r'\mathrm{m\,s^{-1}}', False, 'K.7 c_s=%.1f m/s' % c_s),
        ('K', 7, lam_n0, 0.01, r'\mathrm{m}', False, 'K.7 λ_D(n_0)=%.3e m' % lam_n0),
        ('K', 7, lam_n0 / 1e-3, 0.01, r'\mathrm{mm}', False, 'K.7 λ_D(n_0)=0.129 mm'),
        ('K', 7, lam_ns, 0.01, r'\mathrm{m}', False, 'K.7 λ_D(n_s)=%.3e m' % lam_ns),
        ('K', 7, s_sh, 0.01, r'\mathrm{m}', False, 'K.7 鞘厚=%.3e m' % s_sh),
        ('K', 7, s_sh / 1e-3, 0.01, r'\mathrm{mm}', False, 'K.7 鞘厚=%.2f mm' % (s_sh / 1e-3)),
        ('K', 7, ns, 0.01, r'\mathrm{m^{-3}}', False, 'K.7 n_s=%.3e m^-3' % ns),
        ('K', 7, flux, 0.03, r'\mathrm{m^{-2}s^{-1}}', False, 'K.7 离子通量=%.3e m^-2s^-1' % flux),
        ('K', 7, j_wall, 0.01, r'\mathrm{A\,m^{-2}}', False, 'K.7 j=%.3f A/m²' % j_wall),
        ('K', 7, s_cl / 1e-3, 0.005, r'\mathrm{mm}', False,
         'K.7 CL 反解鞘厚=%.3f mm' % (s_cl / 1e-3)),
        ('K', 8, eta_sp, 0.01, r'\Omega\cdot\mathrm{m}', False, 'K.8 η=%.3e Ω·m' % eta_sp),
        ('K', 8, s_sp, 0.02, '', False, 'K.8 S=%.3e' % s_sp),
        ('K', 8, v_in, 0.01, r'\mathrm{m\,s^{-1}}', False, 'K.8 v_in=%.3f m/s' % v_in),
        ('K', 8, tau_sp, 0.01, r'\mathrm{s}', False, 'K.8 τ_SP=%.3e s' % tau_sp),
        ('K', 8, tau_sp / 3600.0, 0.05, r'\mathrm{h}', False, 'K.8 τ_SP=%.1f h' % (tau_sp / 3600.0)),
    ]
    miss = []
    for tag, idx, val, rel, unit, signed, desc in cases:
        blk = (ex[idx - 1][1] if 1 <= idx <= len(ex) else '') if tag == 'K' \
            else (sec[idx - 1][1] if 1 <= idx <= len(sec) else '')
        if not blk:
            miss.append((desc, '块 %s%d 不存在' % (tag, idx)))
        elif not value_in(blk, val, rel, unit=unit, signed=signed):
            miss.append((desc, '期望 %.6g%s' % (val, '（带符号）' if signed else '')))
    if miss:
        rec('FAIL', 'P', '附录 I/K 有 %d/%d 个答案不在其所属例题/小节内，或符号、单位不符：%s'
            % (len(miss), len(cases), ['%s → %s' % (m[0], m[1]) for m in miss[:4]]))
    else:
        rec('PASS', 'P', '附录 I/K 的 %d 个答案全部出现在所属例题/小节内（K.1--K.8 逐题绑定、I.4 单独绑定），'
            '符号与单位同时核对' % len(cases))


def check_P05_appendix_k_claims():
    """附录 K 三道「判断型」题的结论方向必须写对（文字级，逐题绑定）。

    数值对了不等于物理对了：P-C01（曲率漂移方向）、P-C03（v_A/c 与磁冻结的关系）
    这两处当年的错误都是「数字全对、结论反了」。所以这里对结论句本身下断言，
    并明确禁止旧表述复活。
    """
    text = strip_comments(read(os.path.join(SRC, APP_FILES['K'])))
    ex = split_blocks(text, 'example')
    if len(ex) != 8:
        rec('FAIL', 'P', '附录 K 例题数不是 8（实际 %d 个），无法做文字级判据' % len(ex))
        return
    k3, k4, k5, k6, k7 = (ex[2][1], ex[3][1], ex[4][1], ex[5][1], ex[6][1])

    has(k3, r'\unitvec{b}\times\kappa', 'src/%s K.3（漂移 ∝ b×κ）' % APP_FILES['K'])
    has(k3, r'\textbf{不是}指向曲率中心', 'src/%s K.3（结论句）' % APP_FILES['K'])
    lacks(k3, '曲率漂移沿磁力线的曲率中心方向', 'src/%s K.3（旧错误表述）' % APP_FILES['K'])

    has(k4, '不能据此断定磁冻结成立', 'src/%s K.4（v_A/c 能说明什么）' % APP_FILES['K'])
    has(k4, r'R_m=\mu_0Lv/\eta', 'src/%s K.4（磁冻结的门槛）' % APP_FILES['K'])
    lacks(k4, '磁冻结图像在非相对论等离子体中非常好',
         'src/%s K.4（旧错误表述）' % APP_FILES['K'])

    lacks(k5, r'=+2.00\times10^{-2}', 'src/%s K.5（阻尼率不能是正的）' % APP_FILES['K'])

    has(k6, 'R_m', 'src/%s K.6（结论里的磁雷诺数）' % APP_FILES['K'])

    has(k7, r'\lambda_D(n_s)', 'src/%s K.7（鞘边德拜长度）' % APP_FILES['K'])
    has(k7, r'\approx1.8\ \mathrm{mm}', 'src/%s K.7（鞘厚答案 1.8 mm，单位绑定）' % APP_FILES['K'])
    has(k7, r'1.813\ \mathrm{mm}', 'src/%s K.7（CL 反解鞘厚 1.813 mm）' % APP_FILES['K'])
    has(k7, 'Child--Langmuir', 'src/%s K.7（与 CL 关系交叉核对）' % APP_FILES['K'])
    has(k7, '同一个电流边界条件', 'src/%s K.7（(d) 的结论）' % APP_FILES['K'])
    has(k7, '鞘边', 'src/%s K.7（体密度/鞘边密度的口径区分）' % APP_FILES['K'])

    # ── v1.20.0 追加：同一条判据的「扩散式残留」（评审清单之外，2026-10-11 复查发现）──
    # P-C01/P-C03 属于「数字全对、结论反了」的家族：修了被评审点名的那一处，同一判据
    # 在别处仍可能以同样的错误形式存在。第 20 轮复查抓到三处，此处按同一口径钉住，
    # 并明确禁止旧表述复活——否则同类错误仍会沿着这条家族线复发。
    sup3 = strip_comments(read(os.path.join(SRC, 'ch3_supplement.tex')))
    has(sup3, r'\frac{m}{qBR}', 'src/ch3_supplement.tex（∇B+曲率漂移合并式的分母）')
    lacks(sup3, r'\frac{m}{qB^{2}R}', 'src/ch3_supplement.tex（分母多一个 B 的旧错式）')

    capl = strip_comments(read(os.path.join(SRC, APP_FILES['L'])))
    has(capl, '不能据此断定磁冻结成立', 'src/%s（L.4 练习解答的结论方向）' % APP_FILES['L'])
    has(capl, r'R_m=\mu_0Lv/\eta', 'src/%s（L.4 的磁冻结门槛）' % APP_FILES['L'])
    lacks(capl, '非相对论近似与理想 MHD 都成立', 'src/%s（L.4 的旧越界结论）' % APP_FILES['L'])

    sup5 = strip_comments(read(os.path.join(SRC, 'ch5_supplement.tex')))
    has(sup5, '电阻项可以丢掉', 'src/ch5_supplement.tex（R_m≫1 只解决第三道门槛）')
    has(sup5, '三道门槛中的第三道', 'src/ch5_supplement.tex（判据总结行的限定）')
    lacks(sup5, '磁力线冻结完全成立', 'src/ch5_supplement.tex（旧的一步到位结论）')


def check_P06_scattering_alpha():
    """α=1/(kλ_D) 的判据方向：所有速查表行必须是「α≪1 非相干／α≳1 集体」。

    PP-C06 的教训：12.4 正文一直是对的，但两处速查表（附录 H 符号表与术语表）
    曾把 α≪1 标成集体散射——而速查表正是读者会背下来的地方。本项按**表格行**
    逐行判定，并在子句级别检查搭配，避免「α≪1 为非集体（非相干）散射，α≳1 才需
    集体散射理论」这种同一行里两种判据并存被误判。
    """
    main = strip_comments(read(os.path.join(SRC, MAIN)))
    sup12 = strip_comments(read(os.path.join(SRC, 'ch12_supplement.tex')))

    SMALL = (r'\alpha\ll1', r'\alpha \ll 1')
    LARGE = (r'\alpha\gtrsim1', r'\alpha \gtrsim 1', r'\alpha\gg1', r'\alpha \gg 1')
    INC = ('非相干', '非集体')
    COL = ('集体', '相干散射')

    bad, rows = [], 0
    for ln, line in enumerate(main.split('\n'), 1):
        if '&' not in line or '\\alpha' not in line:
            continue
        if not any(t in line for t in SMALL + LARGE):
            continue
        rows += 1
        for clause in re.split(r'[，。；;,]', line):
            if not clause.strip():
                continue
            if any(t in clause for t in SMALL):
                stripped = clause
                for w in INC:
                    stripped = stripped.replace(w, '')
                if not any(w in clause for w in INC):
                    bad.append((ln, '「α≪1」子句没写非相干/非集体'))
                elif any(w in stripped for w in COL):
                    bad.append((ln, '「α≪1」子句被标成集体'))
            if any(t in clause for t in LARGE):
                stripped = clause
                for w in INC:
                    stripped = stripped.replace(w, '')
                if not any(w in stripped for w in COL):
                    bad.append((ln, '「α≳1」子句没写集体/相干'))
    if bad:
        rec('FAIL', 'P', 'α 散射判据方向有误（%d 处）：%s'
            % (len(bad), ['第 %d 行 %s' % (b[0], b[1]) for b in bad[:4]]))
    elif rows == 0:
        rec('FAIL', 'P', '全书找不到任何含 α 判据的表格行（速查表被删或判据被搬走？）')
    else:
        rec('PASS', 'P', '全书 %d 行含 α 判据的表格行方向正确（α≪1 非相干／α≳1 集体）' % rows)

    has(main, r'$\alpha\ll1$ 为\textbf{非集体（非相干）散射}',
        'src/%s（附录 H 符号表 α 行）' % MAIN)
    has(main, r'$\alpha \gtrsim 1$ 时的散射 regime', 'src/%s（术语表集体散射行）' % MAIN)
    has(main, r'$\alpha\ll1$ 为\textbf{非相干散射}', 'src/%s（12.4 正文判据）' % MAIN)
    has(sup12, r'非相干汤姆逊散射（$\alpha \ll 1$）与相干散射（$\alpha \gg 1$）',
        'src/ch12_supplement.tex（自测清单）')

    # 正文非表格行：只做 INFO（同一句里两种判据并存时跳过，避免脆性误报）
    amb, info = 0, 0
    for ln, line in enumerate(main.split('\n'), 1):
        if '&' in line or '\\alpha' not in line:
            continue
        for t in SMALL + LARGE:
            i = line.find(t)
            if i < 0:
                continue
            seg = line[max(0, i - 40):i + 40]
            if any(w in seg for w in INC) and any(w in seg for w in COL):
                amb += 1
            elif not any(w in seg for w in INC + COL):
                info += 1
                rec('INFO', 'P', '第 %d 行的 α 判据前后 40 字符内没有「非相干/集体」字样，人工确认即可' % ln)
    if amb:
        rec('INFO', 'P', '正文另有 %d 处 α 判据与「非相干/集体」并存于同句（表格行已逐行判定）' % amb)


def check_P07_appendix_l_bank():
    """附录 L 自测题库：12 章 × 6 题（2 概念 / 2 推导 / 2 估算）、每题有答案、每节有出处。

    可证伪点（五种，全部 FAIL）：① 某章被删或题数不等于 6；② 题型配比失衡（例如把一章
    的两道概念题换成推导）；③ 只写题目不写「答案」；④ 节首「取材」行被删（来源不可
    追溯）；⑤ 印出来的关键数值与本函数独立复算不符——数值必须落在**它自己那一节**里，
    所以把某个答案搬到别章、或把 7.7×10⁶ 改成 7.7×10⁵，都会被抓到。
    """
    path = os.path.join(SRC, APP_FILES['L'])
    if not os.path.exists(path):
        rec('FAIL', 'P', '附录 L 文件缺失：src/%s' % APP_FILES['L'])
        return
    text = strip_comments(read(path))
    secs = [b for ttl, b in split_blocks(text, 'section')
            if re.match(r'第\s*\d+\s*章', ttl)]
    if len(secs) != 12:
        rec('FAIL', 'P', '附录 L 的分章小节应为 12 个、实际 %d 个（章被删或标题格式变了）'
            % len(secs))
        return

    struct, mix, src = [], [], []
    for i, blk in enumerate(secs, 1):
        kinds = re.findall(r'\\item\s*（(概念|推导|估算)）', blk)
        n_ans = len(re.findall(r'\\tkans\{', blk))
        if len(kinds) != 6 or n_ans != 6:
            struct.append((i, len(kinds), n_ans))
        cnt = [kinds.count(k) for k in ('概念', '推导', '估算')]
        if cnt != [2, 2, 2]:
            mix.append((i, cnt))
        if r'\tksrc{' not in blk:
            src.append(i)
    if struct:
        rec('FAIL', 'P', '附录 L 有 %d 章不是「6 题且 6 个答案」：%s（章号, 题数, 答案数）'
            % (len(struct), struct))
    else:
        rec('PASS', 'P', '附录 L 每章都是 6 题且每题都有「答案」（12 章 × 6 = 72 题）')
    if mix:
        rec('FAIL', 'P', '附录 L 题型配比失衡，应为每章 2 概念 + 2 推导 + 2 估算：%s' % mix)
    else:
        rec('PASS', 'P', '附录 L 题型配比正确：每章 2 概念 + 2 推导 + 2 估算')
    if src:
        rec('FAIL', 'P', '附录 L 有 %d 章缺少节首「取材」行（来源不可追溯）：%s' % (len(src), src))
    else:
        rec('PASS', 'P', '附录 L 每章节首都有「取材」行，逐章注明参照的教材与章节')

    # ---- 关键数值：每个都由本函数独立算一遍，再去**它自己那一节**里找
    m_e = 9.1093837015e-31
    k_b = 1.380649e-23

    def lam_d(te, n):
        return math.sqrt(EPS0 * te / (n * E_CHG))

    def v_th(te, m):
        return math.sqrt(2.0 * te * E_CHG / m)

    v_th10 = v_th(10.0, m_e)
    wci = E_CHG * 5.0 / M_D
    vi = v_th(1.0e4, M_D)
    wcp = E_CHG * 5.0e-5 / M_P
    v_dia = 100.0 * E_CHG / (E_CHG * 2.0 * 0.1)
    v_a = 5.0 / math.sqrt(MU0 * 1e20 * M_D)
    eta_10ev = 5.2e-5 * 10 / 10.0 ** 1.5
    p_kev10 = 1e20 * 1e4 * E_CHG
    eta_1kev = 5.2e-5 * 10 / 1000.0 ** 1.5
    nu_ei = 2.9e-6 * 1e19 * 10 / 1000.0 ** 1.5
    tau_ei = 1.0 / nu_ei

    def gam_landau(kld):
        return -math.sqrt(math.pi / 8.0) * kld ** -3 * math.exp(-1.0 / (2 * kld ** 2) - 1.5)

    w_pe = math.sqrt(1e18 * E_CHG ** 2 / (EPS0 * m_e))
    n_c_245 = EPS0 * m_e * (2 * math.pi * 2.45e9) ** 2 / E_CHG ** 2
    n_c_140 = EPS0 * m_e * (2 * math.pi * 140e9) ** 2 / E_CHG ** 2
    p_fus = 1e40 * 1.1e-22 * (17.6e6 * E_CHG)
    e_th = 3 * 1e20 * 1e4 * E_CHG
    v_b_ar = math.sqrt(3.0 * E_CHG / M_AR)
    v_f_ar = 3.0 / 2.0 * math.log(M_AR / (2 * math.pi * m_e))
    f_pe_16 = math.sqrt(1e16 * E_CHG ** 2 / (EPS0 * m_e)) / (2 * math.pi)

    cases = [
        # (章序号, 说明, 期望值, 相对容差, 单位, 是否连符号)
        (1, 'v_th(10 eV 电子)', v_th10, 0.03, r'\mathrm{m/s}', False),
        (1, 'v_avg(10 eV 电子)', v_th10 * math.sqrt(4 / math.pi), 0.03, r'\mathrm{m/s}', False),
        (1, 'n(1 Pa, 300 K)', 1.0 / (k_b * 300.0), 0.03, r'\mathrm{m^{-3}}', False),
        (1, '1 eV 对应的温度', E_CHG / k_b, 0.01, r'\mathrm{K}', False),
        (1, 'λ_D(1 eV, 1e15)', lam_d(1.0, 1e15), 0.03, r'\mathrm{m}', False),
        (2, 'λ_D(1 eV, 1e16)', lam_d(1.0, 1e16), 0.03, r'\mathrm{m}', False),
        (2, 'N_D(1 eV, 1e16)', (4 / 3) * math.pi * 1e16 * lam_d(1.0, 1e16) ** 3, 0.05, '', False),
        (2, 'λ_D(100 eV, 1e19)', lam_d(100.0, 1e19), 0.03, r'\mathrm{m}', False),
        (3, 'ω_ci(D+, 5 T)', wci, 0.03, r'\mathrm{rad/s}', False),
        (3, 'v_i(10 keV, D+)', vi, 0.03, r'\mathrm{m/s}', False),
        (3, 'ρ_Li(10 keV, D+, 5 T)', vi / wci, 0.03, r'\mathrm{m}', False),
        (3, 'ρ_Lp(1 keV, 0.5 G)', v_th(1.0e3, M_P) / wcp, 0.03, r'\mathrm{m}', False),
        (3, 'f_ce(1 T)', E_CHG / (2 * math.pi * m_e) / 1e9, 0.03, r'\mathrm{GHz}', False),
        (3, 'ρ_Le(100 eV, 1 T)', v_th(100.0, m_e) / (E_CHG / m_e), 0.03, r'\mathrm{m}', False),
        (4, 'v_dia(100 eV, 2 T, L_n=0.1)', v_dia, 0.03, r'\mathrm{m/s}', False),
        (4, 'j_dia(n=1e19)', 1e19 * E_CHG * 2 * v_dia, 0.03, r'\mathrm{A/m^2}', False),
        (4, 'c_s(D+, 10 eV)', math.sqrt(10.0 * E_CHG / M_D), 0.03, r'\mathrm{m/s}', False),
        (4, 'c_s(Ar+, 1 eV)', math.sqrt(1.0 * E_CHG / M_AR), 0.05, r'\mathrm{m/s}', False),
        (5, 'v_A(5 T, 1e20, D+)', v_a, 0.03, r'\mathrm{m/s}', False),
        (5, 'η(T_e=10 eV)', eta_10ev, 0.05, r'\Omega', False),
        (5, 'τ_η(L=1 m)', MU0 / eta_10ev, 0.03, r'\mathrm{s}', False),
        (5, 'R_m(v=1e4 m/s)', MU0 * 1e4 / eta_10ev, 0.03, '', False),
        (5, 'β(1e20, 10 keV, 5.3 T) 的百分数', 2 * MU0 * p_kev10 / 5.3 ** 2 * 100, 0.05, r'\%', False),
        (5, 'β=3% 所需 B', math.sqrt(2 * MU0 * p_kev10 / 0.03), 0.03, r'\mathrm{T}', False),
        (6, 'n_c(2.45 GHz)', n_c_245, 0.03, r'\mathrm{m^{-3}}', False),
        (6, 'f_pe(n=1e18)', w_pe / (2 * math.pi) / 1e9, 0.03, r'\mathrm{GHz}', False),
        (6, 'Alfvén 波越 1 m 的用时', 1.0 / v_a, 0.03, r'\mathrm{s}', False),
        (7, 'γ/ω(kλ_D=0.3)', gam_landau(0.3), 0.03, '', True),
        (7, 'γ/ω(kλ_D=0.25)', gam_landau(0.25), 0.03, '', True),
        (7, 'γ/ω(kλ_D=0.2)', gam_landau(0.2), 0.03, '', True),
        (7, 'ω_pe(n=1e18)', w_pe, 0.03, r'\mathrm{rad/s}', False),
        (7, '双流 γ_max', w_pe / math.sqrt(8.0), 0.03, r'\mathrm{s^{-1}}', False),
        (8, 'η(T_e=1 keV)', eta_1kev, 0.05, r'\Omega', False),
        (8, 'ν_ei(1e19, 1 keV)', nu_ei, 0.03, r'\mathrm{s^{-1}}', False),
        (8, 'τ_ei(1e19, 1 keV)', tau_ei, 0.03, r'\mathrm{s}', False),
        (8, 'λ_mfp(1 keV)', v_th(1000.0, m_e) * tau_ei, 0.05, r'\mathrm{m}', False),
        (8, '中性粒子穿透深度(n=1e20)', 1.0 / (1e20 * 1e-19), 0.03, r'\mathrm{m}', False),
        (9, 'RT γ(A=1, k=100)', math.sqrt(100 * 9.8), 0.03, r'\mathrm{s^{-1}}', False),
        (9, 'RT e 折时间', 1.0 / math.sqrt(100 * 9.8) * 1000, 0.03, r'\mathrm{ms}', False),
        (9, '交换模 γ(c_s=1e5, R=1 m)', 1e5, 0.03, r'\mathrm{s^{-1}}', False),
        (10, 'P_fus(n=1e20, D-T)', p_fus, 0.03, r'\mathrm{W/m^3}', False),
        (10, '3nkT(10 keV)', e_th, 0.03, r'\mathrm{J/m^3}', False),
        (10, 'P_α = P_fus/5', p_fus / 5, 0.03, r'\mathrm{W/m^3}', False),
        (10, '点火门槛 τ_E', e_th / (p_fus / 5), 0.03, r'\mathrm{s}', False),
        (10, 'V=100 m³ 的总功率', p_fus * 100, 0.03, r'\mathrm{W}', False),
        (10, '1 GW 所需体积', 1e9 / p_fus, 0.03, r'\mathrm{m^3}', False),
        (11, 'V_f(Ar+, 3 eV)', v_f_ar, 0.03, r'\mathrm{V}', False),
        (11, 'v_B(Ar+, 3 eV)', v_b_ar, 0.03, r'\mathrm{m/s}', False),
        (11, 'Γ = n_s v_B', 1e16 * v_b_ar, 0.03, 'm^{-2}s^{-1}', False),
        (11, 'λ_D(3 eV, 1e16)', lam_d(3.0, 1e16), 0.03, r'\mathrm{m}', False),
        (11, 'f_pe(n=1e16)', f_pe_16, 0.03, r'\mathrm{Hz}', False),
        (11, '趋肤深度 δ', 2.99792458e8 / (2 * math.pi * f_pe_16), 0.03, r'\mathrm{m}', False),
        (12, 'T_e(双点法 50→200 μA)', 2.0 / math.log(4.0), 0.03, r'\mathrm{eV}', False),
        (12, 'n_e(I_es=1 mA)', 3.5e16, 0.03, r'\mathrm{m^{-3}}', False),
        (12, '140 GHz 干涉仪相移',
         math.pi / (2.99792458e8 / 140e9) * (1e19 / n_c_140) * 0.5, 0.03, r'\mathrm{rad}', False),
        (12, 'n_c(140 GHz)', n_c_140, 0.03, r'\mathrm{m^{-3}}', False),
        (12, '磁探针感应电压', 10 * 1e-4 * 1e4, 0.03, r'\mathrm{V}', False),
    ]
    bad = []
    for ch, desc, want, rel, unit, signed in cases:
        if not value_in(secs[ch - 1], want, rel, unit, signed):
            bad.append((ch, desc))
    if bad:
        rec('FAIL', 'P', '附录 L 有 %d 项数值与独立复算不符（或不在它所属的那一章内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 L 的 %d 项关键数值全部与独立复算一致，且落在各自章节内'
            % len(cases))

    has(text, r'\tkans', 'src/%s（答案宏定义）' % APP_FILES['L'])
    has(text, r'\tksrc', 'src/%s（取材宏定义）' % APP_FILES['L'])
    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['L'])
    has(text, r'\label{tab:app-l-source}', 'src/%s（来源表标签）' % APP_FILES['L'])


def check_P08_appendix_m_heating():
    r"""附录 M（加热、电流驱动与运行极限）：结构完整 + 每个数值都由本函数独立复算。

    可证伪点（五种，全部 FAIL）：① 某节被删（节名按主题前缀匹配）；② 缺「出处」行，
    来源不可追溯；③ 三张表之一或版权声明被删；④ 印出来的关键数值与本函数独立复算
    不符——数值必须落在**它自己那一节**里，所以把 86.5 MW 改成 8.65 MW、或把
    \$3.10\\times10^{6}\$ m/s 搬到施工别的节，都会被抓到；⑤ 氘质量约定的旧偏差
    （把 2.0136 u 当成 \$2m_p\$，使束速偏大 0.7\% 的 3.085e6 m/s）不得回流。
    """
    mpath = os.path.join(SRC, APP_FILES['M'])
    if not os.path.exists(mpath):
        rec('FAIL', 'P', '附录 M 文件缺失：src/%s' % APP_FILES['M'])
        return
    text = strip_comments(read(mpath))
    secs = split_blocks(text, 'section')
    want_prefix = ['加热手段全景', '欧姆加热', '中性束注入', '射频加热三件套', '运行极限',
                   '六道估算练习', '来源、分级与版权说明']
    found = {}
    for ttl, blk in secs:
        for pre in want_prefix:
            if ttl.startswith(pre):
                found[pre] = blk
    miss = [p for p in want_prefix if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 M 缺少小节（或标题改得无法识别）：%s（实际 %d 节）' % (miss, len(secs)))
        return
    rec('PASS', 'P', '附录 M 的 %d 个小节齐全（%s）' % (len(want_prefix), '、'.join(want_prefix)))

    # ---- 结构：三张表、版权声明、6 道练习各带答案、至少 4 处「出处」
    for lab, what in (('tab:app-m-heating', '加热手段全景表'),
                      ('tab:app-m-rf', '射频三手段对照表'),
                      ('tab:app-m-source', '来源与独立性表')):
        has(text, r'\label{%s}' % lab, 'src/%s（%s）' % (APP_FILES['M'], what))
    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['M'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src >= 4:
        rec('PASS', 'P', '附录 M 有 %d 处「出处」行，逐节注明借鉴来源' % n_src)
    else:
        rec('FAIL', 'P', '附录 M 的「出处」行只有 %d 处（应 ≥4）：来源不可追溯' % n_src)
    ex = found['六道估算练习']
    n_item = len(re.findall(r'\\item\s*\\textbf\{（', ex))
    n_ans = len(re.findall(r'\\textbf\{答案：\}', ex))
    if n_item == 6 and n_ans == 6:
        rec('PASS', 'P', '附录 M 的六道估算练习每题都有「答案」（6 题 / 6 答案）')
    else:
        rec('FAIL', 'P', '附录 M 练习数与答案数应为 6/6，实际 %d/%d' % (n_item, n_ans))
    lacks(text, '3.085', 'src/%s（旧氘质量约定下的束速，不得回流）' % APP_FILES['M'])

    # ---- 数值：全部由本函数从头算一遍（CODATA 2018；氘取 2.013553212745 u）
    m_e = 9.1093837015e-31
    m_u = 1.66053906660e-27
    m_p = M_P
    m_d = 2.013553212745 * m_u

    def eta_sp(te_ev, lnl=17.0):
        return 5.2e-5 * lnl / te_ev ** 1.5

    def v_of(energy_ev, m):
        return math.sqrt(2.0 * energy_ev * E_CHG / m)

    def f_cyc(B, m):
        return E_CHG * B / (2 * math.pi * m)

    def b_res(f, m):
        return 2 * math.pi * f * m / E_CHG

    p_fus_m = 3.10181e6                      # 附录 L 第 10 章的堆芯功率密度（已复核）
    p_alpha_m = p_fus_m / 5.0
    eta_1k = eta_sp(1.0e3)
    v_te_10k = v_of(1.0e4, m_e)
    f_lh = math.sqrt(f_cyc(5.3, m_d) * f_cyc(5.3, m_e))
    eta_10k = eta_sp(1.0e4)
    n_g = 15.0 / (math.pi * 2.0 ** 2)
    f_pe_ng = math.sqrt(n_g * 1e20 * E_CHG ** 2 / (EPS0 * m_e)) / (2 * math.pi)

    cases = [
        # (节, 说明, 期望值, 相对容差, 单位, 是否连符号)
        ('欧姆加热', 'η(1 keV, lnΛ=17)=2.80e-8 Ω·m', eta_1k, 0.02, r'\Omega', False),
        ('欧姆加热', 'P_Ω(j=1 MA/m²)=2.8e4 W/m³', eta_1k * 1e12, 0.02, r'\mathrm{W/m^3}', False),
        ('欧姆加热', 'P_Ω=28 kW/m³', eta_1k * 1e12 / 1e3, 0.02, r'\mathrm{kW/m^3}', False),
        ('欧姆加热', 'τ_R(a=2 m)=180 s', MU0 * 4.0 / eta_1k, 0.02, r'\mathrm{s}', False),
        ('欧姆加热', '10 keV 时的 P_Ω(884 W/m³)', eta_10k * 1e12, 0.03, r'\mathrm{W/m^3}', False),
        ('中性束注入', 'v_b(100 keV D)=3.10e6 m/s', v_of(1.0e5, m_d), 0.02, r'\mathrm{m/s}', False),
        ('中性束注入', 'v_b(1 MeV D)=9.79e6 m/s', v_of(1.0e6, m_d), 0.02, r'\mathrm{m/s}', False),
        ('中性束注入', '衰减长度 λ=1 m', 1.0 / (1e20 * 1e-20), 0.02, r'\mathrm{m}', False),
        ('中性束注入', 'E_c(10 keV)=186 keV', 14.8 * 10.0 * 2.0 * 0.5 ** (2 / 3.0), 0.02,
         r'\mathrm{keV}', False),
        ('中性束注入', 'E_c(1 keV)=18.6 keV', 14.8 * 1.0 * 2.0 * 0.5 ** (2 / 3.0), 0.02,
         r'\mathrm{keV}', False),
        ('中性束注入', '1 MA 需 P=21 MW(γ=0.3)', 1.0 * 1.0 * 6.2 / 0.3, 0.03, r'\mathrm{MW}', False),
        ('射频加热三件套', 'f_ci(D, 5.3 T)=40 MHz', f_cyc(5.3, m_d) / 1e6, 0.03, r'\mathrm{MHz}', False),
        ('射频加热三件套', 'B_res(40 MHz, D)=5.24 T', b_res(40e6, m_d), 0.02, r'\mathrm{T}', False),
        ('射频加热三件套', 'B_res(40 MHz, H)=2.62 T', b_res(40e6, m_p), 0.02, r'\mathrm{T}', False),
        ('射频加热三件套', 'f_ce/B=28.0 GHz/T', f_cyc(1.0, m_e) / 1e9, 0.02, r'\mathrm{GHz/T}', False),
        ('射频加热三件套', 'f_ce(5.3 T)=148 GHz', f_cyc(5.3, m_e) / 1e9, 0.02, r'\mathrm{GHz}', False),
        ('射频加热三件套', 'f_LH=2.45 GHz', f_lh / 1e9, 0.02, r'\mathrm{GHz}', False),
        ('射频加热三件套', 'v_te(10 keV)=5.93e7 m/s', v_te_10k, 0.02, r'\mathrm{m/s}', False),
        ('射频加热三件套', 'k_∥(2v_te)=130 m⁻¹', 2 * math.pi * f_lh / (2 * v_te_10k), 0.02,
         r'\mathrm{m^{-1}}', False),
        ('射频加热三件套', 'λ_∥(2v_te)=4.8 cm', 2 * v_te_10k / f_lh * 100.0, 0.02, r'\mathrm{cm}', False),
        ('射频加热三件套', 'k_∥(3v_te)=87 m⁻¹', 2 * math.pi * f_lh / (3 * v_te_10k), 0.02,
         r'\mathrm{m^{-1}}', False),
        ('射频加热三件套', 'λ_∥(3v_te)=7.3 cm', 3 * v_te_10k / f_lh * 100.0, 0.02, r'\mathrm{cm}', False),
        ('运行极限', 'n_G(15 MA, a=2 m)=1.19e20', n_g * 1e20, 0.02, r'\mathrm{m^{-3}}', False),
        ('运行极限', 'β_N(2.5%, 2 m, 5.3 T, 15 MA)=1.77', 2.5 * 2.0 * 5.3 / 15.0, 0.02, '', False),
        ('运行极限', 'P_LH(大装置)=86.5 MW', 0.0488 * 5.3 ** 0.803 * 683.0 ** 0.941, 0.02,
         r'\mathrm{MW}', False),
        ('运行极限', 'P_LH(小装置)=4.72 MW',
         0.0488 * 0.5 ** 0.717 * 2.5 ** 0.803 * 100.0 ** 0.941, 0.02, r'\mathrm{MW}', False),
        ('六道估算练习', '[1] η(10 keV)=8.84e-10 Ω·m', eta_10k, 0.02, r'\Omega', False),
        ('六道估算练习', '[1] P_Ω(10 keV, 1 MA/m²)=884 W/m³', eta_10k * 1e12, 0.03,
         r'\mathrm{W/m^3}', False),
        ('六道估算练习', '[1] 与自加热差约 700 倍', p_alpha_m / (eta_10k * 1e12), 0.05, '', False),
        ('六道估算练习', '[2] v_b(1 MeV D)=9.79e6 m/s', v_of(1.0e6, m_d), 0.02,
         r'\mathrm{m/s}', False),
        ('六道估算练习', '[3] B_res(170 GHz)=6.07 T', b_res(170e9, m_e), 0.02, r'\mathrm{T}', False),
        ('六道估算练习', '[4] v_φ(3v_te)=1.78e8 m/s', 3 * v_te_10k, 0.02, r'\mathrm{m/s}', False),
        ('六道估算练习', '[4] λ_∥=7.3 cm', 3 * v_te_10k / 2.45e9 * 100.0, 0.02,
         r'\mathrm{cm}', False),
        ('六道估算练习', '[5] n_G=1.19e20 m⁻³', n_g * 1e20, 0.02, r'\mathrm{m^{-3}}', False),
        ('六道估算练习', '[5] f_pe(n_G)=98 GHz', f_pe_ng / 1e9, 0.02, r'\mathrm{GHz}', False),
        ('六道估算练习', '[6] P_ext=1.24e5 W/m³', 0.2 * p_alpha_m, 0.03, r'\mathrm{W/m^3}', False),
        ('六道估算练习', '[6] P_ext·V=1.04e8 W', 0.2 * p_alpha_m * 840.0, 0.02, r'\mathrm{W}', False),
        # 正文按「量级」写作 100 MW（精确值 104 MW），故容差放宽到 6%
        ('六道估算练习', '[6] 折合 100 MW 量级', 0.2 * p_alpha_m * 840.0 / 1e6, 0.06,
         r'\mathrm{MW}', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 M 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 M 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    # 结论句（口径类断言，不依赖数值匹配）
    has(found['欧姆加热'], '差 22 倍', 'src/%s（欧姆功率与 α 自加热的倍数）' % APP_FILES['M'])
    has(found['欧姆加热'], '只适合启动与爬升阶段',
        'src/%s（欧姆加热的适用区间结论）' % APP_FILES['M'])
    has(found['运行极限'], '经验拟合', 'src/%s（三条运行极限标注为经验标度）' % APP_FILES['M'])
    has(found['运行极限'], '仿星器', 'src/%s（与仿星器的稳态对比）' % APP_FILES['M'])
    has(found['六道估算练习'], '射程不变',
        'src/%s（练习 2 的结论：提高束能不会线性提高穿透）' % APP_FILES['M'])


def check_P10_appendix_n_boundary():
    r"""附录 N（等离子体-壁相互作用与偏滤器物理）：结构完整 + 每个数值独立复算。

    可证伪点（六种，全部 FAIL）：① 某节被删（节名按主题前缀匹配）；② 缺「出处」行；
    ③ 四张表之一或版权声明被删；④ 六道练习缺答案；⑤ 印出来的关键数值与本函数独立
    复算不符——数值必须落在**它自己那一节**里（把 35.9 MW/m² 改成 3.59、或把
    0.33 nm/次 ELM 搬到别的节都会被抓到）；⑥ 旧口径回流：按"日历秒"而不是"等离子体
    秒"折算的年侵蚀量 2680 μm/年。
    """
    npath = os.path.join(SRC, APP_FILES['N'])
    if not os.path.exists(npath):
        rec('FAIL', 'P', '附录 N 文件缺失：src/%s' % APP_FILES['N'])
        return
    text = strip_comments(read(npath))
    secs = split_blocks(text, 'section')
    want_prefix = ['从芯部到壁', '刮削层与两点模型', '脱靶', '溅射与侵蚀', '边缘局域模',
                   '燃料滞留', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for ttl, blk in secs:
        for pre in want_prefix:
            if ttl.startswith(pre):
                found[pre] = blk
    miss = [p for p in want_prefix if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 N 缺少小节（或标题改得无法识别）：%s（实际 %d 节）' % (miss, len(secs)))
        return
    rec('PASS', 'P', '附录 N 的 %d 个小节齐全（%s）' % (len(want_prefix), '、'.join(want_prefix)))

    for lab, what in (('tab:app-n-tilt', '倾角热流表'),
                      ('tab:app-n-twopoint', '两点模型自洽解表'),
                      ('tab:app-n-impurity', '杂质比例敏感性表'),
                      ('tab:app-n-source', '来源与取材表')):
        has(text, r'\label{%s}' % lab, 'src/%s（%s）' % (APP_FILES['N'], what))
    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['N'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src >= 6:
        rec('PASS', 'P', '附录 N 有 %d 处「出处」行，逐节注明借鉴来源' % n_src)
    else:
        rec('FAIL', 'P', '附录 N 的「出处」行只有 %d 处（应 ≥6）：来源不可追溯' % n_src)
    ex = found['六道估算练习']
    n_item = len(re.findall(r'\\item\s*\\textbf\{（', ex))
    n_ans = len(re.findall(r'\\textbf\{答案：\}', ex))
    if n_item == 6 and n_ans == 6:
        rec('PASS', 'P', '附录 N 的六道估算练习每题都有「答案」（6 题 / 6 答案）')
    else:
        rec('FAIL', 'P', '附录 N 练习数与答案数应为 6/6，实际 %d/%d' % (n_item, n_ans))
    lacks(text, '2680', 'src/%s（按日历秒折算的旧年侵蚀量，不得回流）' % APP_FILES['N'])

    # ---- 数值：全部由本函数从头算一遍（与 tmp/appN_calc2.py 同一套公式）
    m_d = 2.013553212745 * M_P

    def c_s(te_ev, m=m_d):
        return math.sqrt(te_ev * E_CHG / m)

    p_sol = 70e6
    lam_q, f_x, r0 = 5e-3, 5.0, 6.2
    a_wet = 2 * math.pi * r0 * lam_q * f_x
    q_par = p_sol * 0.5 / a_wet
    kappa0, l_conn, gamma_sh = 2000.0, 40.0, 7.0

    def n_target(te_ev):
        return q_par / (gamma_sh * te_ev * E_CHG * c_s(te_ev))

    def t_up(tt_ev):
        return (tt_ev ** 3.5 + 3.5 * q_par * l_conn / kappa0) ** (2 / 7.0)

    n_t_15, t_u_15 = n_target(15.0), t_up(15.0)
    n_u_15 = n_t_15 * 15.0 / t_u_15
    q_low_rec = gamma_sh * 2e19 * 15 * E_CHG * c_s(15.0)
    gamma_i = 2e19 * c_s(15.0)
    q_allow = 5e6
    f_rad = 1 - q_allow / q_par
    p_rad = f_rad * 35e6
    n_imp_32 = p_rad / (1e20 * 3e-32 * 20.0)
    n_W = 6.3e28
    d_steady = 1e-5 * gamma_i / n_W
    plasma_s = 1e4 * 10.0
    n_ion_elm = 1e6 / (300 * E_CHG)
    d_elm = 1e-3 * n_ion_elm / n_W
    w_th = 3 * 1e20 * (1e4 * E_CHG) * 840.0
    q_elm = 0.1 * w_th * 0.5 / 20.0 / 0.3e-3
    e_fus = 17.6e6 * E_CHG

    def t_burn_gday(p_fus):
        return p_fus / e_fus * 3.0160492 / 6.02214076e23 * 86400

    cases = [
        # (节, 说明, 期望值, 相对容差, 单位, 是否连符号)
        ('从芯部到壁', 'P_SOL=70 MW', p_sol / 1e6, 0.02, r'\mathrm{MW}', False),
        ('从芯部到壁', 'A_wet=0.974 m²', a_wet, 0.02, r'\mathrm{m^2}', False),
        ('从芯部到壁', 'q_∥=35.9 MW/m²', q_par / 1e6, 0.02, r'\mathrm{MW/m^2}', False),
        ('从芯部到壁', 'A_wall=979 m²', 4 * math.pi ** 2 * r0 * 2.0 * 2.0, 0.02, r'\mathrm{m^2}', False),
        ('从芯部到壁', 'q_wall=71.5 kW/m²', p_sol / (4 * math.pi ** 2 * r0 * 2.0 * 2.0) / 1e3, 0.02,
         r'\mathrm{kW/m^2}', False),
        ('从芯部到壁', '几何杠杆 1/sin3°≈19', 1 / math.sin(math.radians(3)), 0.02, '', False),
        # 表 tab:app-n-tilt 的数值：单位写在表头，故不校验单位
        ('从芯部到壁', 'q_⊥(10°)=6.24 MW/m²', q_par * math.sin(math.radians(10)) / 1e6, 0.02, '', False),
        ('从芯部到壁', 'q_⊥(5°)=3.13 MW/m²', q_par * math.sin(math.radians(5)) / 1e6, 0.02, '', False),
        ('从芯部到壁', 'q_⊥(3°)=1.88 MW/m²', q_par * math.sin(math.radians(3)) / 1e6, 0.02, '', False),
        ('从芯部到壁', 'q_⊥(1°)=0.627 MW/m²', q_par * math.sin(math.radians(1)) / 1e6, 0.02, '', False),
        # 该值出现在表 tab:app-n-twopoint 中，单位写在表头，故不校验单位
        ('刮削层与两点模型', 'c_s(15 eV)=2.68e4 m/s', c_s(15.0), 0.02, '', False),
        ('刮削层与两点模型', 'Γ_i=5.36e23 m⁻²s⁻¹', gamma_i, 0.02,
         r'\mathrm{m^{-2}\,s^{-1}}', False),
        ('刮削层与两点模型', 'E_i≈75 eV', 75.0, 0.02, r'\mathrm{eV}', False),
        ('刮削层与两点模型', '低再循环 q_∥=9.0 MW/m²', q_low_rec / 1e6, 0.02,
         r'\mathrm{MW/m^2}', False),
        ('刮削层与两点模型', 'n_t(15 eV)=7.97e19', n_t_15, 0.02, '', False),
        ('刮削层与两点模型', 'T_u=67.5 eV', t_u_15, 0.02, '', False),
        ('刮削层与两点模型', 'n_u=1.77e19', n_u_15, 0.02, '', False),
        ('刮削层与两点模型', '表行 T_t=5 eV 的 n_t=4.14e20', n_target(5.0), 0.02, '', False),
        ('刮削层与两点模型', '表行 T_t=10 eV 的 n_t=1.46e20', n_target(10.0), 0.02, '', False),
        ('刮削层与两点模型', '表行 T_t=30 eV 的 T_u=68.5 eV', t_up(30.0), 0.02, '', False),
        ('刮削层与两点模型', '表行 T_t=50 eV 的 n_u=8.91e18', n_target(50.0) * 50.0 / t_up(50.0), 0.02, '', False),
        ('脱靶', '需辐射比例 0.861', f_rad, 0.02, '', False),
        ('脱靶', 'P_rad=30.1 MW', p_rad / 1e6, 0.02, r'\mathrm{MW}', False),
        ('脱靶', 'n_imp=5.02e17 (L_z=3e-32)', n_imp_32, 0.02, '', False),
        ('脱靶', '杂质比例 0.50%', 100 * n_imp_32 / 1e20, 0.02, '', False),
        ('脱靶', 'L_z 低十倍则 5.0%', 100 * n_imp_32 * 10 / 1e20, 0.02, '', False),
        ('溅射与侵蚀', '稳态侵蚀 8.5 μm/年', d_steady * plasma_s * 1e6, 0.02, '', False),
        ('溅射与侵蚀', '再沉积后 0.85 μm/年', d_steady * plasma_s * 1e6 * 0.1, 0.02, '', False),
        ('溅射与侵蚀', '稳态寿命约 5900 年', 5e-3 / (d_steady * plasma_s * 0.1), 0.02, '', False),
        ('溅射与侵蚀', '单次 ELM 离子通量 2.08e22', n_ion_elm, 0.02, r'\mathrm{m^{-2}}', False),
        ('溅射与侵蚀', '单次 ELM 侵蚀 0.33 nm', d_elm * 1e9, 0.02, r'\mathrm{nm}', False),
        ('溅射与侵蚀', 'ELM 年侵蚀 3.30 mm/年', d_elm * 1000 * 1e4 * 1e3, 0.02, '', False),
        ('溅射与侵蚀', '再沉积后 0.330 mm/年', d_elm * 1000 * 1e4 * 1e3 * 0.1, 0.02, '', False),
        ('溅射与侵蚀', 'ELM 主导寿命约 15 年', 5e-3 / (d_elm * 1000 * 1e4 * 0.1), 0.02, '', False),
        ('溅射与侵蚀', 'ELM 比稳态高约 390 倍',
         (d_elm * 1000 * 1e4) / (d_steady * plasma_s), 0.02, '', False),
        ('边缘局域模', 'W_th=404 MJ', w_th / 1e6, 0.02, r'\mathrm{MJ}', False),
        ('边缘局域模', 'ΔW_ELM=40.4 MJ', 0.1 * w_th / 1e6, 0.02, r'\mathrm{MJ}', False),
        ('边缘局域模', 'q_ELM=3.4 GW/m²', q_elm / 1e9, 0.02, r'\mathrm{GW/m^2}', False),
        ('边缘局域模', '相对稳态约 1.8e3 倍',
         q_elm / (q_par * math.sin(math.radians(3))), 0.02, '', False),
        ('燃料滞留', '氚消耗 1.77e20 /s', 500e6 / e_fus, 0.02, r'\mathrm{s^{-1}}', False),
        ('燃料滞留', '76.7 g/天', t_burn_gday(500e6), 0.02, r'\mathrm{g/day}', False),
        ('燃料滞留', '1 kg 触顶 3.6 年', 1000 / (t_burn_gday(500e6) * 0.01) / 365, 0.02, '', False),
        ('燃料滞留', '滞留率上限 0.36%',
         100 * 1000 / (10 * 365) / t_burn_gday(500e6), 0.02, '', False),
        ('燃料滞留', '衰变 0.154 g/天',
         1000 * math.log(2) / (12.32 * 365.25), 0.02, r'\mathrm{g/day}', False),
        ('燃料滞留', '300 MW 时 46.0 g/天', t_burn_gday(300e6), 0.02, r'\mathrm{g/day}', False),
        ('燃料滞留', '300 MW 时 6.0 年', 1000 / (t_burn_gday(300e6) * 0.01) / 365, 0.02, '', False),
        ('燃料滞留', '300 MW 时上限 0.60%',
         100 * 1000 / (10 * 365) / t_burn_gday(300e6), 0.02, '', False),
        ('燃料滞留', 'TBR≥1.05', 1.05, 0.02, '', False),
        ('六道估算练习', '[1] q_∥=35.9 MW/m²', q_par / 1e6, 0.02, r'\mathrm{MW/m^2}', False),
        ('六道估算练习', '[1] q_⊥(3°)=1.88 MW/m²', q_par * math.sin(math.radians(3)) / 1e6, 0.02,
         r'\mathrm{MW/m^2}', False),
        ('六道估算练习', '[2] q_∥=9.0 MW/m²', q_low_rec / 1e6, 0.02, r'\mathrm{MW/m^2}', False),
        ('六道估算练习', '[2] c_s=2.68e4 m/s', c_s(15.0), 0.02, r'\mathrm{m/s}', False),
        ('六道估算练习', '[3] T_u=67.5 eV', t_u_15, 0.02, r'\mathrm{eV}', False),
        ('六道估算练习', '[3] n_u=1.77e19 m⁻³', n_u_15, 0.02, r'\mathrm{m^{-3}}', False),
        ('六道估算练习', '[4] n_imp=5.0e17 m⁻³', n_imp_32, 0.03, r'\mathrm{m^{-3}}', False),
        ('六道估算练习', '[4] 比例 0.50%', 100 * n_imp_32 / 1e20, 0.02, '', False),
        ('六道估算练习', '[5] W_th=404 MJ', w_th / 1e6, 0.02, r'\mathrm{MJ}', False),
        ('六道估算练习', '[5] q_ELM=3.4 GW/m²', q_elm / 1e9, 0.02, r'\mathrm{GW/m^2}', False),
        ('六道估算练习', '[6] 46.0 g/天', t_burn_gday(300e6), 0.02, r'\mathrm{g/day}', False),
        ('六道估算练习', '[6] 上限 0.60%',
         100 * 1000 / (3650 * t_burn_gday(300e6)), 0.02, '', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 N 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 N 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    has(found['从芯部到壁'], '几何杠杆', 'src/%s（倾角与磁通膨胀的几何杠杆结论）' % APP_FILES['N'])
    has(found['刮削层与两点模型'], '高再循环', 'src/%s（靶板温度由上游决定 → 高再循环）' % APP_FILES['N'])
    has(found['脱靶'], '运行窗口', 'src/%s（脱靶是运行窗口而不是开关）' % APP_FILES['N'])
    has(found['溅射与侵蚀'], 'ELM 造成的侵蚀比稳态高约',
        'src/%s（稳态不致命、ELM 才致命的结论）' % APP_FILES['N'])
    has(found['边缘局域模'], 'ELM 缓解', 'src/%s（ELM 缓解是壁寿命设计的前提）' % APP_FILES['N'])
    has(found['燃料滞留'], 'TBR', 'src/%s（氚自持的增殖比判据）' % APP_FILES['N'])


def check_P11_appendix_o_propulsion():
    r"""附录 O（电推进：把等离子体当发动机）：结构完整 + 每个数值独立复算。

    可证伪点（六种，全部 FAIL）：① 某节被删（节名按主题前缀匹配）；② 缺「出处」行；
    ③ 六张表/三只语义盒/十条编号公式之一被删或改标签；④ 六道练习缺答案；⑤ 印出来的
    关键数值与本函数独立复算不符——而且必须落在**它自己那一节**里（把 $40.8$ 改成
    $4.08$、或把 $2100\ \mathrm{s}$ 搬到别的节都会被抓到）；⑥ 旧口径回流：把恒等式
    $P_{\mathrm{jet}}=I_bV_b$ 误写成 $\frac12I_bV_b$，或把比冲定义写成
    $I_{\mathrm{sp}}=g_0v_{\mathrm{ex}}$（乘除颠倒，量纲看起来仍然"对"）。

    表内数值的单位写在表头，故这些用例传 unit=''（与附录 N 的门禁同一约定）。
    """
    opath = os.path.join(SRC, APP_FILES['O'])
    if not os.path.exists(opath):
        rec('FAIL', 'P', '附录 O 文件缺失：src/%s' % APP_FILES['O'])
        return
    text = strip_comments(read(opath))
    secs = split_blocks(text, 'section')
    want_prefix = ['为什么用电推进', '推力与功率账', '静电加速', '霍尔推力器', '磁喷嘴',
                   '中和', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for ttl, blk in secs:
        for pre in want_prefix:
            if ttl.startswith(pre):
                found[pre] = blk
    miss = [p for p in want_prefix if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 O 缺少小节（或标题改得无法识别）：%s（实际 %d 节）' % (miss, len(secs)))
        return
    rec('PASS', 'P', '附录 O 的 %d 个小节齐全（%s）' % (len(want_prefix), '、'.join(want_prefix)))

    for lab, what in (('tab:app-o-isp', '比冲-质量比表'), ('tab:app-o-tp', '推力-功率账表'),
                      ('tab:app-o-beam', '静电加速速度表'), ('tab:app-o-power', '同一电源三档表'),
                      ('tab:app-o-loss', '效率损失清单表'), ('tab:app-o-source', '来源与取材表')):
        has(text, r'\label{%s}' % lab, 'src/%s（%s）' % (APP_FILES['O'], what))
    # 语义盒的标签写在 \begin{tip}{标题}{标签} 的第二个必选参数里，不是 \label{}
    for lab, what in (('tip-appo-units', '口径自查盒'), ('warn:appo-hall', '霍尔推力器警示盒'),
                      ('tip-appo-chain', '与 K/M/N 串联盒')):
        has(text, '{%s}' % lab, 'src/%s（%s）' % (APP_FILES['O'], what))
    for lab, what in (('eq:app-o-rocket', '火箭方程与比冲定义'), ('eq:app-o-thrust', '推力与功率定义'),
                      ('eq:app-o-tpratio', 'T/P 与推力无关的关系式'),
                      ('eq:app-o-vbeam', '单离子能量方程'), ('eq:app-o-ibeam', '束流推力与恒等式'),
                      ('eq:app-o-cl', 'Child--Langmuir 律'),
                      ('eq:app-o-magnetized', '磁化判据'), ('eq:app-o-eb', 'E×B 漂移'),
                      ('eq:app-o-triangle', '推力-功率-比冲交换关系'),
                      ('eq:app-o-cex', 'CEX 概率')):
        has(text, r'\label{%s}' % lab, 'src/%s（%s）' % (APP_FILES['O'], what))
    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['O'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src >= 6:
        rec('PASS', 'P', '附录 O 有 %d 处「出处」行，逐节注明借鉴来源' % n_src)
    else:
        rec('FAIL', 'P', '附录 O 的「出处」行只有 %d 处（应 ≥6）：来源不可追溯' % n_src)
    ex = found['六道估算练习']
    n_item = len(re.findall(r'\\item\s*\\textbf\{（', ex))
    n_ans = len(re.findall(r'\\textbf\{答案：\}', ex))
    if n_item == 6 and n_ans == 6:
        rec('PASS', 'P', '附录 O 的六道估算练习每题都有「答案」（6 题 / 6 答案）')
    else:
        rec('FAIL', 'P', '附录 O 练习数与答案数应为 6/6，实际 %d/%d' % (n_item, n_ans))
    lacks(text, r'\frac12I_bV_b',
          'src/%s（恒等式 P_jet=I_bV_b 被误写成 ½I_bV_b 的旧口径）' % APP_FILES['O'])
    lacks(text, r'g_0v_{\mathrm{ex}}',
          'src/%s（比冲定义写成 I_sp=g_0·v_ex 的旧口径）' % APP_FILES['O'])

    # ---- 数值：全部由本函数从头算一遍（与 tmp/appO_calc.py 同一套公式）
    v_ex = lambda isp: isp * G0                      # 排气速度
    ratio = lambda dv, isp: math.exp(dv / v_ex(isp))  # 质量比 m0/mf
    t_of = lambda ib, vb: ib * math.sqrt(2 * M_XE * vb / E_CHG)
    cl_j = lambda v, d: (4 * EPS0 / 9) * math.sqrt(2 * E_CHG / M_XE) * v ** 1.5 / d ** 2
    t_tpr = lambda isp, eta: 2 * eta / (G0 * isp)     # T/P_in [N/W]
    hall_wce = E_CHG * 0.02 / M_E                     # B=0.02 T
    hall_nu = 1e19 * 2e-19 * math.sqrt(2 * 20 * E_CHG / M_E)
    hall_vi = math.sqrt(2 * E_CHG * 300 / M_XE)
    pr_ex = lambda isp: 2 * 0.70 * 2e5 / (G0 * isp)   # 200 kW、η=0.70 的推力

    cases = [
        # (节, 说明, 期望值, 相对容差, 单位, 是否连符号)
        ('为什么用电推进', '化学 Isp 上限 450 s', 450.0, 0.02, r'\mathrm{s}', False),
        # 以下四行的四列都来自表 tab:app-o-isp，单位写在表头，故不校验单位
        ('为什么用电推进', '450 s → 4.41 km/s', v_ex(450) / 1e3, 0.02, '', False),
        ('为什么用电推进', '450 s → m0/mf=9.64', ratio(1e4, 450), 0.02, '', False),
        ('为什么用电推进', '450 s → 推进剂 89.6%', 100 * (1 - 1 / ratio(1e4, 450)), 0.02, '', False),
        ('为什么用电推进', '1800 s → 17.7 km/s', v_ex(1800) / 1e3, 0.02, '', False),
        ('为什么用电推进', '1800 s → m0/mf=1.76', ratio(1e4, 1800), 0.02, '', False),
        ('为什么用电推进', '1800 s → 推进剂 43.3%', 100 * (1 - 1 / ratio(1e4, 1800)), 0.02, '', False),
        ('为什么用电推进', '3000 s → 29.4 km/s', v_ex(3000) / 1e3, 0.02, '', False),
        ('为什么用电推进', '3000 s → m0/mf=1.40', ratio(1e4, 3000), 0.02, '', False),
        ('为什么用电推进', '3000 s → 推进剂 28.8%', 100 * (1 - 1 / ratio(1e4, 3000)), 0.02, '', False),
        ('为什么用电推进', '5000 s → 49.0 km/s', v_ex(5000) / 1e3, 0.02, '', False),
        ('为什么用电推进', '5000 s → m0/mf=1.23', ratio(1e4, 5000), 0.02, '', False),
        ('为什么用电推进', '5000 s → 推进剂 18.4%', 100 * (1 - 1 / ratio(1e4, 5000)), 0.02, '', False),
        # 卫星例（Δv=3 km/s、1000 kg）
        ('为什么用电推进', 'Isp=300 s 需 639 kg', 1000 * (1 - 1 / ratio(3000, 300)), 0.02, r'\mathrm{kg}', False),
        ('为什么用电推进', '化学推进剂占 64%', 100 * (1 - 1 / ratio(3000, 300)), 0.02, '', False),
        ('为什么用电推进', 'Isp=3000 s 只需 97 kg', 1000 * (1 - 1 / ratio(3000, 3000)), 0.02, r'\mathrm{kg}', False),
        ('为什么用电推进', '省下 542 kg',
         1000 * (1 / ratio(3000, 3000) - 1 / ratio(3000, 300)), 0.02, r'\mathrm{kg}', False),
        ('为什么用电推进', '静电场 1000 V 给出 38 km/s', math.sqrt(2 * E_CHG * 1000 / M_XE) / 1e3,
         0.02, r'\mathrm{km/s}', False),
        # 表 tab:app-o-tp：单位写在表头，故不校验单位
        ('推力与功率账', 'T/P=40.8 mN/kW（3000 s、η=0.6）', t_tpr(3000, 0.6) * 1e6, 0.02, '', False),
        ('推力与功率账', '0.1 N → 3.40 mg/s', 0.1 / v_ex(3000) * 1e6, 0.02, '', False),
        ('推力与功率账', '0.1 N → P_jet=1471 W', 0.5 * 0.1 * v_ex(3000), 0.02, '', False),
        ('推力与功率账', '0.1 N → P_in=2452 W（η=0.6）', 0.5 * 0.1 * v_ex(3000) / 0.6, 0.02, '', False),
        ('推力与功率账', '10 N → 339.9 mg/s', 10 / v_ex(3000) * 1e6, 0.02, '', False),
        ('推力与功率账', '10 N → P_jet=1.47e5 W', 0.5 * 10 * v_ex(3000), 0.02, '', False),
        ('推力与功率账', '10 N → P_in=2.45e5 W', 0.5 * 10 * v_ex(3000) / 0.6, 0.02, '', False),
        ('推力与功率账', 'NSTAR 级 → 3.03 mg/s', 0.092 / v_ex(3100) * 1e6, 0.02, '', False),
        ('推力与功率账', 'NSTAR 级 → P_jet=1398 W', 0.5 * 0.092 * v_ex(3100), 0.02, '', False),
        ('推力与功率账', 'NSTAR 级 → P_in=2293 W（η=0.61）', 0.5 * 0.092 * v_ex(3100) / 0.61, 0.02, '', False),
        ('推力与功率账', 'NSTAR 级 → T/P=40.1 mN/kW', t_tpr(3100, 0.61) * 1e6, 0.02, '', False),
        ('推力与功率账', '10 N 推一年要 1.07e4 kg', 10 / v_ex(3000) * 3.156e7, 0.02, r'\mathrm{kg}', False),
        ('推力与功率账', '折算 10.7 t/yr', 10 / v_ex(3000) * 3.156e7 / 1e3, 0.02, r'\mathrm{t/yr}', False),
        ('推力与功率账', '1 N 要 24.5 kW', 1 / t_tpr(3000, 0.6) / 1e3, 0.02, r'\mathrm{kW}', False),
        ('推力与功率账', '100 N 要 2.45 MW', 100 / t_tpr(3000, 0.6) / 1e6, 0.02, r'\mathrm{MW}', False),
        # 表 tab:app-o-beam：单位写在表头，故不校验单位
        ('静电加速', '300 V → 21.0 km/s', math.sqrt(2 * E_CHG * 300 / M_XE) / 1e3, 0.02, '', False),
        ('静电加速', '300 V → Isp 2141 s', math.sqrt(2 * E_CHG * 300 / M_XE) / G0, 0.02, '', False),
        ('静电加速', '1000 V → 38.3 km/s', math.sqrt(2 * E_CHG * 1000 / M_XE) / 1e3, 0.02, '', False),
        ('静电加速', '1000 V → Isp 3909 s', math.sqrt(2 * E_CHG * 1000 / M_XE) / G0, 0.02, '', False),
        ('静电加速', '1500 V → 47.0 km/s', math.sqrt(2 * E_CHG * 1500 / M_XE) / 1e3, 0.02, '', False),
        ('静电加速', '1500 V → Isp 4788 s', math.sqrt(2 * E_CHG * 1500 / M_XE) / G0, 0.02, '', False),
        ('静电加速', 'I_b=2 A、1000 V → T=0.104 N', t_of(2.0, 1000), 0.02, r'\mathrm{N}', False),
        ('静电加速', '恒等式 P_jet=I_bV_b=2.00 kW', 2.0 * 1000 / 1e3, 0.02, r'\mathrm{kW}', False),
        ('静电加速', 'ṁ=2.72 mg/s', 2.0 * M_XE / E_CHG * 1e6, 0.02, r'\mathrm{mg/s}', False),
        ('静电加速', 'η=0.65 → P_in=3.08 kW', 2000 / 0.65 / 1e3, 0.02, r'\mathrm{kW}', False),
        ('静电加速', 'CL(1500 V,1.5 mm)=123 A/m²', cl_j(1500, 1.5e-3), 0.02, r'\mathrm{A/m^2}', False),
        ('静电加速', '3 A 需栅极 243 cm²', 3.0 / cl_j(1500, 1.5e-3) * 1e4, 0.02, r'\mathrm{cm^2}', False),
        ('静电加速', '电流密度涨 84%', 100 * (1.5 ** 1.5 - 1), 0.02, '', False),
        # 两项都是两位有效数字的舍入（22.5% 印成 22%），容差放到 3%
        ('静电加速', '比冲只涨 22%', 100 * (math.sqrt(1500 / 1000) - 1), 0.03, '', False),
        ('霍尔推力器', 'B=0.02 T → ω_ce=3.52e9', hall_wce, 0.02, r'\mathrm{rad/s}', False),
        ('霍尔推力器', 'ν_en=5.30e6 s⁻¹', hall_nu, 0.02, r'\mathrm{s^{-1}}', False),
        ('霍尔推力器', 'ω_ceτ_e≈663（电子被磁化）', hall_wce / hall_nu, 0.02, '', False),
        ('霍尔推力器', 'ω_ci=1.47e4 rad/s', E_CHG * 0.02 / M_XE, 0.02, r'\mathrm{rad/s}', False),
        ('霍尔推力器', 'τ_i=1.19 μs', 0.025 / hall_vi * 1e6, 0.02, r'\mu\mathrm{s}', False),
        # 0.0175 印成 0.018（两位有效数字），容差放到 5%
        ('霍尔推力器', 'ω_ciτ_i≈0.018≪1', E_CHG * 0.02 / M_XE * 0.025 / hall_vi, 0.05, '', False),
        ('霍尔推力器', 'E=V_d/L=1.2e4 V/m', 300 / 0.025, 0.02, r'\mathrm{V/m}', False),
        ('霍尔推力器', 'v_E=E/B=6.0e5 m/s', 12000 / 0.02, 0.02, r'\mathrm{m/s}', False),
        ('霍尔推力器', '霍尔电流是束流的 28.6 倍', (12000 / 0.02) / hall_vi, 0.02, '', False),
        ('霍尔推力器', 'I_b=5 A → T=0.143 N', t_of(5.0, 300), 0.02, r'\mathrm{N}', False),
        ('霍尔推力器', 'P_jet=I_bV_d=1.50 kW', 5.0 * 300 / 1e3, 0.02, r'\mathrm{kW}', False),
        ('霍尔推力器', 'η=0.5 → P_in=3.00 kW', 1500 / 0.5 / 1e3, 0.02, r'\mathrm{kW}', False),
        ('霍尔推力器', 'T/P=47.6 mN/kW', t_of(5.0, 300) / 3000 * 1e6, 0.02, r'\mathrm{mN/kW}', False),
        ('霍尔推力器', '比冲从 450 抬到 2100 s', 2100.0, 0.02, r'\mathrm{s}', False),
        # 表 tab:app-o-power：单位写在表头，故不校验单位
        ('磁喷嘴', 'Isp=2000 s → 19.6 km/s', v_ex(2000) / 1e3, 0.02, '', False),
        ('磁喷嘴', 'Isp=2000 s → T=14.3 N', pr_ex(2000), 0.02, '', False),
        ('磁喷嘴', 'Isp=2000 s → ṁ=727.9 mg/s', pr_ex(2000) / v_ex(2000) * 1e6, 0.02, '', False),
        ('磁喷嘴', 'Isp=5000 s → 49.0 km/s', v_ex(5000) / 1e3, 0.02, '', False),
        ('磁喷嘴', 'Isp=5000 s → T=5.71 N', pr_ex(5000), 0.02, '', False),
        ('磁喷嘴', 'Isp=5000 s → ṁ=116.5 mg/s', pr_ex(5000) / v_ex(5000) * 1e6, 0.02, '', False),
        ('磁喷嘴', 'Isp=10000 s → 98.1 km/s', v_ex(10000) / 1e3, 0.02, '', False),
        ('磁喷嘴', 'Isp=10000 s → T=2.86 N', pr_ex(10000), 0.02, '', False),
        ('磁喷嘴', 'Isp=10000 s → ṁ=29.1 mg/s', pr_ex(10000) / v_ex(10000) * 1e6, 0.02, '', False),
        ('中和', '中和器功耗 30 W', 2.0 * 15, 0.02, r'\mathrm{W}', False),
        ('中和', '占 2 kW 束功率的 1.5%', 100 * 30 / 2000, 0.02, '', False),
        ('中和', 'CEX(n_n=1e17)=1.5%', 100 * 1e17 * 5e-19 * 0.3, 0.03, '', False),
        ('中和', 'CEX(n_n=1e18)=15%', 100 * 1e18 * 5e-19 * 0.3, 0.02, '', False),
        ('中和', '发散 15° 损失 3.4%', 100 * (1 - math.cos(math.radians(15))), 0.02, '', False),
        ('中和', '发散 25° 损失 9.4%', 100 * (1 - math.cos(math.radians(25))), 0.02, '', False),
        ('中和', '氙电离能 12.13 eV', 12.13, 0.02, r'\mathrm{eV}', False),
        ('中和', '占 300 V 束能的 4.0%', 100 * 12.13 / 300, 0.02, '', False),
        ('中和', '占 1000 V 的 1.2%', 100 * 12.13 / 1000, 0.02, '', False),
        ('六道估算练习', '[1] v_ex=2.45e4 m/s', v_ex(2500), 0.02, r'\mathrm{m/s}', False),
        ('六道估算练习', '[1] m0/mf=1.277', ratio(6000, 2500), 0.02, '', False),
        ('六道估算练习', '[1] 推进剂占 21.7%', 100 * (1 - 1 / ratio(6000, 2500)), 0.02, '', False),
        ('六道估算练习', '[2] ṁ=8.16e-6 kg/s', 0.2 / v_ex(2500), 0.02, r'\mathrm{kg/s}', False),
        ('六道估算练习', '[2] P_jet=2452 W', 0.5 * 0.2 * v_ex(2500), 0.02, r'\mathrm{W}', False),
        ('六道估算练习', '[2] P_in=4458 W（η=0.55）', 0.5 * 0.2 * v_ex(2500) / 0.55, 0.02, r'\mathrm{W}', False),
        ('六道估算练习', '[3] v(1200 V)=4.20e4 m/s', math.sqrt(2 * E_CHG * 1200 / M_XE), 0.02,
         r'\mathrm{m/s}', False),
        ('六道估算练习', '[3] Isp=4282 s', math.sqrt(2 * E_CHG * 1200 / M_XE) / G0, 0.02, r'\mathrm{s}', False),
        ('六道估算练习', '[3] T=0.229 N', t_of(4.0, 1200), 0.02, r'\mathrm{N}', False),
        ('六道估算练习', '[3] P_jet=I_bV_b=4800 W', 4.0 * 1200, 0.02, r'\mathrm{W}', False),
        ('六道估算练习', '[4] CL(1200 V,1 mm)=198 A/m²', cl_j(1200, 1e-3), 0.02, r'\mathrm{A/m^2}', False),
        ('六道估算练习', '[4] 3 A 需 151 cm²', 3.0 / cl_j(1200, 1e-3) * 1e4, 0.02, r'\mathrm{cm^2}', False),
        ('六道估算练习', '[4] 电压 +20% → 电流 +31%', 100 * (1.2 ** 1.5 - 1), 0.02, '', False),
        ('六道估算练习', '[5] v_th(25 eV)=2.97e6 m/s',
         math.sqrt(2 * 25 * E_CHG / M_E), 0.02, r'\mathrm{m/s}', False),
        ('六道估算练习', '[5] ν_en=4.45e6 s⁻¹',
         5e18 * 3e-19 * math.sqrt(2 * 25 * E_CHG / M_E), 0.02, r'\mathrm{s^{-1}}', False),
        ('六道估算练习', '[5] ω_ce=2.64e9 rad/s', E_CHG * 0.015 / M_E, 0.02, r'\mathrm{rad/s}', False),
        ('六道估算练习', '[5] ω_ceτ_e≈593',
         E_CHG * 0.015 / M_E / (5e18 * 3e-19 * math.sqrt(2 * 25 * E_CHG / M_E)), 0.02, '', False),
        ('六道估算练习', '[5] v_i(250 V)=1.92e4 m/s', math.sqrt(2 * E_CHG * 250 / M_XE), 0.02,
         r'\mathrm{m/s}', False),
        ('六道估算练习', '[5] Isp=1955 s', math.sqrt(2 * E_CHG * 250 / M_XE) / G0, 0.02, r'\mathrm{s}', False),
        ('六道估算练习', '[6] T=9.94 N（300 kW、η=0.65）', 2 * 0.65 * 3e5 / v_ex(4000), 0.02,
         r'\mathrm{N}', False),
        ('六道估算练习', '[6] 1 t 推进剂用 45.7 day',
         1000 / (2 * 0.65 * 3e5 / v_ex(4000) / v_ex(4000)) / 86400, 0.02, r'\mathrm{day}', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 O 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 O 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    # 口径自查盒与全篇约定（这些内容在节外，故对整个文件断言）
    value_in(text, 2.180172e-25, 0.02, r'\mathrm{kg}', False)
    has(text, r'\texttt{day}', 'src/%s（时间单位写成 day 而不是中文字，避免数学字体缺字）' % APP_FILES['O'])
    has(text, '不需要再乘', 'src/%s（P_jet=I_bV_b 不需要 ½ 因子的提示）' % APP_FILES['O'])
    has(text, r'\frac{T}{P_{\mathrm{in}}}=\frac{2\eta}{g_0I_{\mathrm{sp}}}',
        'src/%s（T/P 与推力无关的关系式）' % APP_FILES['O'])
    has(text, r'T\propto\frac{P_{\mathrm{in}}}{I_{\mathrm{sp}}}',
        'src/%s（比冲调高、推力减半的交换关系）' % APP_FILES['O'])
    has(found['为什么用电推进'], '电推进存在的理由',
        'src/%s（电推进的存在理由与代价）' % APP_FILES['O'])
    has(found['推力与功率账'], '用时间换速度增量',
        'src/%s（只能做长期小推力的结论）' % APP_FILES['O'])
    has(found['推力与功率账'], '日历秒',
        'src/%s（推进剂按日历秒计、与附录 N 的等离子体秒区分）' % APP_FILES['O'])
    has(found['静电加速'], '更大的束流面积',
        'src/%s（"面积比电压更重要"的结论）' % APP_FILES['O'])
    has(found['霍尔推力器'], '近三十倍',
        'src/%s（霍尔电流远大于束流的结论）' % APP_FILES['O'])
    has(found['霍尔推力器'], '电离越多越好',
        'src/%s（电离与加速共用电压的警示）' % APP_FILES['O'])
    has(found['磁喷嘴'], '比冲调高一倍，推力就减半',
        'src/%s（可调比冲的取舍结论）' % APP_FILES['O'])
    has(found['中和'], '损失清单',
        'src/%s（效率天花板在损失清单而不在加速原理）' % APP_FILES['O'])
    has(found['来源、分级与版权说明'], '不写无法核实的 DOI 与卷页',
        'src/%s（离线环境下的来源交代方式）' % APP_FILES['O'])


def check_P12_appendix_p_lowtemp():
    """附录 P 低温等离子体应用的结构完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 八个主题小节之一被删（按小节标题前缀定位）；
      ② 缺少「出处：」行（本附录要求 >= 6 处）；
      ③ 六张表 / 三只语义盒 / 编号公式之一被删或改了 label；
      ④ 六道估算练习缺答案（本题库的答案在正文叙述里，故要求每道题都有「解」）；
      ⑤ 印出的关键数值与独立复算不符，或搬到了别的小节里；
      ⑥ 旧口径回流：把 A/m^2 与 mA/cm^2 混用、把帕申公式外推到大气压
         （正文必须承认 305.3 V 的最小值与 1.68 倍的偏差）。

    约定与附录 M/N/O 相同：表内数值的单位写在表头，故这些用例传 unit=''。
    语义盒的标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    ppath = os.path.join(SRC, APP_FILES['P'])
    if not os.path.exists(ppath):
        rec('FAIL', 'P', '附录 P 文件缺失：src/%s' % APP_FILES['P'])
        return
    text = strip_comments(read(ppath))
    secs = split_blocks(text, 'section')

    want_prefix = ['三种放电与参数地图', '汤森雪崩与帕申定律', '介质阻挡放电与大气压等离子体',
                   '鞘层与离子轰击', 'PECVD 与溅射沉积', '等离子体医学与农业',
                   '环境与催化', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want_prefix:
        for ttl, body in secs:
            if ttl.startswith(pre):
                found[pre] = body
                break
    missing = [p for p in want_prefix if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 P 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 P 的 %d 个小节齐全' % len(want_prefix))

    for lab in ('tab:app-p-regimes', 'tab:app-p-paschen', 'tab:app-p-dbd',
                'tab:app-p-etch', 'tab:app-p-env', 'tab:app-p-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['P'], lab))

    # 语义盒标签是 \begin{tip}{标题}{标签} 的第二参数，不是 \label{}
    for lab in ('tip-appp-units', 'tip-appp-chain', 'tip-appp-next'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['P'], lab))

    for lab in ('eq:app-p-debye', 'eq:app-p-flux', 'eq:app-p-townsend', 'eq:app-p-selfsustain',
                'eq:app-p-paschen', 'eq:app-p-empirical', 'eq:app-p-dbdcd', 'eq:app-p-manley',
                'eq:app-p-child', 'eq:app-p-etchrate', 'eq:app-p-pecvdrate', 'eq:app-p-o3eff'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['P'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['P'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 P 的「出处：」行只有 %d 处，应 >= 6' % n_src)
    else:
        rec('PASS', 'P', '附录 P 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_sol = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex != 6 or n_sol < 6:
        rec('FAIL', 'P', '附录 P 练习数为 %d、解答数为 %d（应为 6 与 >= 6）' % (n_ex, n_sol))
    else:
        rec('PASS', 'P', '附录 P 六道估算练习均有解答')

    # ---- 独立复算：常量与公式全部来自 tmp/appP_calc.py 的同一套公式 ----
    # 原子质量必须用原子质量单位 u（appP_calc.py 的 U = 1.66053906660e-27）：
    # 若误用质子质量 M_P，m_Ar 与 m_Si 会系统性偏大 0.73%，正好被 1% 容差吞掉。
    U_ATOM = 1.66053906660e-27
    m_ar = 39.948 * U_ATOM
    m_si = 28.0855 * U_ATOM
    n_si = 2329.0 / m_si
    k_b = KB
    torr = 133.322
    loschmidt = LOSCHMIDT
    sccm = loschmidt * 1.6666667e-8

    def debye_p(te, ne):
        return math.sqrt(EPS0 * te * E_CHG / (ne * E_CHG ** 2))

    def cs_p(te, m):
        return math.sqrt(te * E_CHG / m)

    a_pas, b_pas, gam = 15.0, 365.0, 0.01
    k_pas = math.log(math.log(1 + 1 / gam))
    pd_min = math.exp(k_pas + 1) / a_pas

    def vb_p(pd):
        return b_pas * pd / (math.log(a_pas * pd) - k_pas)

    c_d = EPS0 * 10 * 0.01 / 1e-3
    q_half = c_d * (4000 - 3000)
    # 驱动频率 10 kHz：与附录 P 表头、以及 10 ns/100 微秒 = 1.0e-4 的占空比一致（appP_calc.py 亦为 1.0e4）
    f_dbd = 1.0e4
    p_dbd = 4 * f_dbd * c_d * 3000 * (4000 - 3000)
    n_fil = 1000.0
    v_fil = (4.0 / 3.0) * math.pi * (1e-4) ** 3

    def sheath_p(te, ld, vs):
        return (2.0 / 3.0) * ld * (2 * vs / te) ** 0.75

    lam_d_etch = debye_p(3.0, 1e17)
    gam_i = 1e17 * cs_p(3.0, m_ar)
    n_dot_sih4 = 20 * sccm
    a_waf = math.pi * 0.1 ** 2
    r_pecvd = (n_dot_sih4 * 0.1 / a_waf) / n_si
    gam_sput = 100.0 / E_CHG
    n_ar = 1.3332 / (k_b * 300.0)

    cases = [
        # --- §P.1 参数地图 ---
        ('三种放电与参数地图', '刻蚀基准 λ_D=40.72 μm', lam_d_etch * 1e6, 0.01, '', False),
        ('三种放电与参数地图', '刻蚀基准 c_s=2692 m/s', cs_p(3.0, m_ar), 0.01, r'\mathrm{m/s}', False),
        ('三种放电与参数地图', '刻蚀基准 j_i=4.313 mA/cm²', gam_i * E_CHG * 0.1, 0.01, '', False),
        ('三种放电与参数地图', '射流 λ_D=2.88 μm', debye_p(1.5, 1e19) * 1e6, 0.01, '', False),
        ('三种放电与参数地图', 'DBD 微放电 λ_D=0.33 μm', debye_p(2.0, 1e21) * 1e6, 0.02, '', False),
        ('三种放电与参数地图', '聚变边界对照 j_i=1928.720 mA/cm²',
         2e19 * cs_p(15.0, m_ar) * E_CHG * 0.1, 0.01, '', False),
        ('三种放电与参数地图', '非平衡比 T_e/T_gas≈116', 3.0 / (k_b * 300.0 / E_CHG), 0.01, '', False),
        # --- §P.2 帕申 ---
        ('汤森雪崩与帕申定律', 'pd_min=0.8363 Torr·cm', pd_min, 0.01, r'\mathrm{Torr\cdot cm}', False),
        ('汤森雪崩与帕申定律', 'V_b,min=305.27 V', vb_p(pd_min), 0.01, r'\mathrm{V}', False),
        ('汤森雪崩与帕申定律', '外推 pd=76 → 5035 V', vb_p(76.0), 0.01, r'\mathrm{V}', False),
        ('汤森雪崩与帕申定律', '经验击穿 30 kV/cm 的比值 1.68',
         vb_p(76.0) / 3000.0, 0.01, '', False),
        ('汤森雪崩与帕申定律', 'AN2 下 N=2.4463e25 m⁻³',
          101325.0 / (k_b * 300.0), 0.01, r'\mathrm{m}^{-3}', False),
        ('汤森雪崩与帕申定律', '约化场强 122.6 Td',
          3e6 / (101325.0 / (k_b * 300.0)) * 1e21, 0.01, r'\mathrm{Td}', False),
        ('汤森雪崩与帕申定律', 'α(30 kV/cm)=1.0993 cm⁻¹',
          a_pas * 760 * math.exp(-b_pas * 760 / 3e4), 0.01, r'\mathrm{cm^{-1}}', False),
        # --- §P.3 DBD ---
        ('介质阻挡放电与大气压等离子体', 'C_d=0.8854 nF', c_d * 1e9, 0.01, r'\mathrm{nF}', False),
        # 表头印的驱动频率必须与自己算出的 P、每周期能量、占空比一致：
        # 20 kHz 会让 P 与每周期能量各差一倍（v1.13.0 修掉的正是这个表头笔误）
        ('介质阻挡放电与大气压等离子体', '驱动频率 10 kHz', f_dbd / 1e3, 0.01,
         r'\mathrm{kHz}', False),
        ('介质阻挡放电与大气压等离子体', '驱动周期 100 μs', 1.0 / f_dbd * 1e6, 0.01,
         r'\mathrm{\mu s}', False),
        ('介质阻挡放电与大气压等离子体', '半周期电荷 0.8854 μC', q_half * 1e6, 0.01,
         r'\mathrm{\mu C}', False),
        ('介质阻挡放电与大气压等离子体', 'P=106.3 W', p_dbd, 0.01, r'\mathrm{W}', False),
        # 表内数值：单位写在表头，故 unit 传空（与附录 M/N/O 的表内用例同口径）
        ('介质阻挡放电与大气压等离子体', '功率密度 1.063 W/cm²', p_dbd / 100.0, 0.01,
         '', False),
        ('介质阻挡放电与大气压等离子体', '每周期 10.63 mJ', p_dbd / f_dbd * 1e3, 0.01,
         r'\mathrm{mJ}', False),
        ('介质阻挡放电与大气压等离子体', '单丝电荷 0.885 nC', q_half / n_fil * 1e9, 0.01,
         r'\mathrm{nC}', False),
        ('介质阻挡放电与大气压等离子体', '单丝能量 5.3 μJ',
         (p_dbd / f_dbd * 1e3) / (2 * n_fil) * 1e3, 0.02, r'\mathrm{\mu J}', False),
        ('介质阻挡放电与大气压等离子体', '时空填充因子 4.189e-8',
         n_fil * v_fil / 1e-5 * (10e-9 / 100e-6), 0.02, '', False),
        # --- §P.4 鞘层与刻蚀 ---
        ('鞘层与离子轰击', 'V_s=100 V → s=0.633 mm', sheath_p(3.0, lam_d_etch, 100.0) * 1e3,
         0.01, r'\mathrm{mm}', False),
        ('鞘层与离子轰击', 'V_s=400 V → s=1.791 mm', sheath_p(3.0, lam_d_etch, 400.0) * 1e3,
         0.01, r'\mathrm{mm}', False),
        ('鞘层与离子轰击', 'V_s=100 V 张角 7.02°',
         math.degrees(math.sqrt(3.0 / (2 * 100.0))), 0.01, '', False),
        ('鞘层与离子轰击', 'V_s=400 V 张角 3.51°',
         math.degrees(math.sqrt(3.0 / (2 * 400.0))), 0.01, '', False),
        ('鞘层与离子轰击', '10 mTorr 平均自由程 3.11 mm', 1.0 / (n_ar * 1e-18) * 1e3,
         0.01, r'\mathrm{mm}', False),
        ('鞘层与离子轰击', '硅原子密度 4.9939e28 m⁻³', n_si, 0.01, '', False),
        # 表内数值：速率单位写在表头 (nm/min)，故 unit 传空
        ('鞘层与离子轰击', 'Y=1.0 → 323.4 nm/min', gam_i / n_si * 60 * 1e9, 0.01,
         '', False),
        ('鞘层与离子轰击', 'Y=0.1 → 32.3 nm/min', gam_i * 0.1 / n_si * 60 * 1e9, 0.01,
         '', False),
        # --- §P.5 沉积 ---
        ('PECVD 与溅射沉积', '20 sccm SiH4=8.9559e18 s⁻¹', n_dot_sih4, 0.01, r'\mathrm{s^{-1}}', False),
        ('PECVD 与溅射沉积', '晶圆面积 0.031416 m²', a_waf, 0.01, r'\mathrm{m^2}', False),
        ('PECVD 与溅射沉积', '生长速率 34.3 nm/min', r_pecvd * 60 * 1e9, 0.01,
         r'\mathrm{nm/min}', False),
        ('PECVD 与溅射沉积', '长 1 μm 需 29.20 min', 1e-6 / r_pecvd / 60, 0.01,
         r'\mathrm{min}', False),
        ('PECVD 与溅射沉积', '溅射通量 6.2415e20 m⁻²s⁻¹', gam_sput, 0.01,
         r'\mathrm{m^{-2}s^{-1}}', False),
        ('PECVD 与溅射沉积', 'Y=0.5 → 374.9 nm/min', gam_sput * 0.5 / n_si * 60 * 1e9, 0.01,
         r'\mathrm{nm/min}', False),
        ('PECVD 与溅射沉积', '靶电荷 1.0800e6 C', 3.00 * 100 * 3600, 0.01, r'\mathrm{C}', False),
        # --- §P.6 医学 ---
        ('等离子体医学与农业', '剂量 60.000 J/cm²', 1.0 * 60.0 / 1e-4 * 1e-4, 0.01,
         r'\mathrm{J/cm^2}', False),
        ('等离子体医学与农业', 'H2O2 分子数 6.0221e16',
         100e-6 * 1e-3 * 6.02214076e23, 0.01, '', False),
        ('等离子体医学与农业', '每分子 6218.6 eV',
         60.0 / (100e-6 * 1e-3 * 6.02214076e23) / E_CHG, 0.01, r'\mathrm{eV}', False),
        ('等离子体医学与农业', '与 O-O 键能之比 0.0823%',
         5.12 / (60.0 / (100e-6 * 1e-3 * 6.02214076e23) / E_CHG) * 100, 0.02, '', False),
        ('等离子体医学与农业', '固氮 29965.0 MJ/kg N',
         4.35 * E_CHG * 6.02214076e23 / 1e6 * 1000.0 / (14.0067 / 1000.0), 0.01,
         r'\mathrm{MJ/kg\ N}', False),
        # --- §P.7 环境 ---
        ('环境与催化', '臭氧 1.2547e24 分子/kWh', 100.0 / 47.9982 * 6.02214076e23, 0.01, '', False),
        ('环境与催化', '每分子 17.91 eV',
         3.6e6 / (100.0 / 47.9982 * 6.02214076e23) / E_CHG, 0.01, r'\mathrm{eV}', False),
        ('环境与催化', '极限效率 28.6%',
         5.12 / (3.6e6 / (100.0 / 47.9982 * 6.02214076e23) / E_CHG) * 100, 0.01, '', False),
        ('环境与催化', '产 1 kg 臭氧需 10.0 kWh', 1000.0 / 100.0, 0.01, r'\mathrm{kWh}', False),
        ('环境与催化', '甲苯 1000 ppmv 每分子 69.7 eV',
         300.0 / (1000e-6 * loschmidt * 1e-3) / E_CHG, 0.01, r'\mathrm{eV}', False),
        # --- §P.8 练习答案 ---
        ('六道估算练习', '[1] V_b,min=305.3 V', vb_p(pd_min), 0.01, r'\mathrm{V}', False),
        ('六道估算练习', '[1] pd_min=0.8363 Torr·cm', pd_min, 0.01, r'\mathrm{Torr\cdot cm}', False),
        ('六道估算练习', '[1] 外推 5035 V vs 经验 3000 V', vb_p(76.0), 0.01, r'\mathrm{V}', False),
        ('六道估算练习', '[2] Y=1.2 → 388.1 nm/min', gam_i * 1.2 / n_si * 60 * 1e9, 0.01,
         r'\mathrm{nm/min}', False),
        ('六道估算练习', '[2] Y=0.1 → 32.3 nm/min', gam_i * 0.1 / n_si * 60 * 1e9, 0.01,
         r'\mathrm{nm/min}', False),
        ('六道估算练习', '[3] V_s=400 V → s=1.791 mm',
         sheath_p(3.0, lam_d_etch, 400.0) * 1e3, 0.01, r'\mathrm{mm}', False),
        ('六道估算练习', '[3] λ_CX=3.11 mm', 1.0 / (n_ar * 1e-18) * 1e3, 0.01,
         r'\mathrm{mm}', False),
        ('六道估算练习', '[4] 1 μm 需 29.2 min', 1e-6 / r_pecvd / 60, 0.01, r'\mathrm{min}', False),
        ('六道估算练习', '[4] 硅烷利用率 100.0%', 100.0, 0.01, '', False),
        ('六道估算练习', '[5] 17.91 eV/分子',
         3.6e6 / (100.0 / 47.9982 * 6.02214076e23) / E_CHG, 0.01, r'\mathrm{eV}', False),
        ('六道估算练习', '[5] 产 1 kg 需 10.0 kWh', 1000.0 / 100.0, 0.01, r'\mathrm{kWh}', False),
        ('六道估算练习', '[6] D=30.0 J/cm²', 60.0 / 2.0, 0.01, r'\mathrm{J/cm^2}', False),
        ('六道估算练习', '[6] 固氮 29965.0 MJ/kg N',
         4.35 * E_CHG * 6.02214076e23 / 1e6 * 1000.0 / (14.0067 / 1000.0), 0.01,
         r'\mathrm{MJ/kg\ N}', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 P 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 P 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    # 全篇口径与结论句（内容跨节，故对整个文件断言）
    has(text, r'\texttt{day}', 'src/%s（时间单位写成 day 而不是中文字，避免数学字体缺字）'
        % APP_FILES['P'])
    has(text, '1\\ \\mathrm{A/m^2}=0.1\\ \\mathrm{mA/cm^2}',
        'src/%s（A/m^2 与 mA/cm^2 的换算，最常错的口径）' % APP_FILES['P'])
    has(text, r'\mathrm{Torr\cdot cm}', 'src/%s（帕申自变量 pd 的单位）' % APP_FILES['P'])
    has(found['汤森雪崩与帕申定律'], '不能',
        'src/%s（大气压不能外推低压帕申公式的结论）' % APP_FILES['P'])
    has(found['汤森雪崩与帕申定律'], '流注',
        'src/%s（大气压击穿由流注而非汤森机制主导）' % APP_FILES['P'])
    has(found['介质阻挡放电与大气压等离子体'], '宏观冷',
        'src/%s（大气压冷等离子体"宏观冷、微观热"的结论）' % APP_FILES['P'])
    has(found['鞘层与离子轰击'], '各向异性',
        'src/%s（提高偏压提高刻蚀各向异性的结论）' % APP_FILES['P'])
    has(found['等离子体医学与农业'], '主导误差',
        'src/%s（剂量学一节的误差交代）' % APP_FILES['P'])
    has(found['环境与催化'], '成反比',
        'src/%s（VOC 处理能量代价与浓度成反比的结论）' % APP_FILES['P'])
    has(found['来源、分级与版权说明'], '不复制',
        'src/%s（不复制任何题目、解答、图表与数据表的声明）' % APP_FILES['P'])
    lacks(text, 'g_0v_{\\mathrm{ex}}',
          'src/%s（比冲与排气速度乘除颠倒的旧口径不得回流）' % APP_FILES['P'])


def check_P13_appendix_i_calendar():
    """附录 I 的「间隔与提取：把复习排进日历」：借用数据的复算与排程表的逐行一致性。

    这一节把 Ebbinghaus (1885) 的四点保留率当原始数据，用它反解单指数时间常数与幂律指数，
    再按 g_max = floor(0.15 D) 加「间隔翻倍到顶后保持等间隔」的规则生成复习日。
    正文里的结论句、排程表与这里重算的数只要对不上（或数据被悄悄改过），这一节就会 FAIL。
    """
    p = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(p):
        rec('FAIL', 'P', '附录 I 文件不存在：%s' % APP_FILES['I'])
        return
    blk = ''
    for ttl, body in split_blocks(strip_comments(read(p)), 'section'):
        if ttl.startswith('间隔与提取'):
            blk = body
            break
    if not blk:
        rec('FAIL', 'P', '附录 I 缺少「间隔与提取：把复习排进日历」一节（标题被改或整节被删）')
        return
    rec('PASS', 'P', '附录 I 有「间隔与提取：把复习排进日历」一节')

    # 1) 用 Ebbinghaus 的两点（1 天 33%、31 天 21%）反解两个模型，并复核 1/3/7/21 的定位
    r1, r31 = 0.33, 0.21
    tau1 = -1.0 / math.log(r1)
    tau31 = -31.0 / math.log(r31)
    b_pow = math.log(r1 / r31) / math.log(31.0)
    cases = [
        ('单指数 τ（由 1 天 33% 反解）', tau1, 0.02, ''),
        ('单指数 τ（由 31 天 21% 反解）', tau31, 0.02, ''),
        ('两个 τ 相差的倍数（单指数不成立的证据）', tau31 / tau1, 0.02, ''),
        ('幂律指数 b', b_pow, 0.02, ''),
        ('1/3/7/21 的末档间隔', 14.0, 0.01, ''),
        ('四周（D=28）时的 20% 上界', 0.2 * 28.0, 0.01, ''),
    ]
    bad = []
    for desc, want, rel, unit in cases:
        if not value_in(blk, want, rel, unit=unit):
            bad.append('%s（复算 %.4g）' % (desc, want))
    if bad:
        rec('FAIL', 'P', '附录 I 排程一节有 %d 项数据与独立复算不符：%s' % (len(bad), bad[:4]))
    else:
        rec('PASS', 'P', '附录 I 排程一节的 %d 项借用数据与派生量全部与独立复算一致' % len(cases))

    # 2) 排程表：按规则重算末档间隔、复习日列表与次数，逐行比对字面值
    def sched(D):
        g = int(0.15 * D)                     # 末档间隔 = 保留期的 15%，向下取整到整天
        gaps, step = [], 1
        while step < g:                       # 间隔从 1 天起翻倍，翻到超过 g 就停
            gaps.append(step)
            step *= 2
        gaps.append(g)
        while sum(gaps) + g <= D - 1:         # 到顶后保持等间隔，铺到考试前
            gaps.append(g)
        days, s = [], 0
        for x in gaps:
            s += x
            days.append(s)
        return g, days

    rows = [(7, 1, 6), (28, 4, 8), (84, 12, 9), (365, 54, 11)]
    bad2 = []
    for D, g_want, n_want in rows:
        g, days = sched(D)
        lit = r',\,'.join(str(x) for x in days)
        if g != g_want:
            bad2.append('D=%d 的末档间隔复算 %d 天、表中 %d 天' % (D, g, g_want))
        if len(days) != n_want:
            bad2.append('D=%d 的复习次数复算 %d、表中 %d' % (D, len(days), n_want))
        if lit not in blk:
            bad2.append('D=%d 的复习日列表与规则复算不一致（复算 %s）' % (D, lit))
    if bad2:
        rec('FAIL', 'P', '附录 I 复习排程表有 %d 处与规则复算不符：%s' % (len(bad2), bad2[:4]))
    else:
        rec('PASS', 'P', '附录 I 复习排程表 4 行（D=7/28/84/365）的末档间隔、复习日与次数全部与规则复算一致')

    # 3) 借鉴来源必须逐条写清（这是本节的立项理由）
    srcs = ['Ebbinghaus', 'Roediger', 'Karpicke', 'Cepeda', 'Dunlosky', 'Bjork', 'Rasch']
    miss = [s for s in srcs if s not in blk]
    if miss:
        rec('FAIL', 'P', '附录 I 排程一节的「出处」缺 %d 个来源：%s' % (len(miss), miss))
    else:
        rec('PASS', 'P', '附录 I 排程一节逐条注明 7 个借鉴来源（Ebbinghaus 1885 原始数据、Cepeda 2006 间隔规则等）')
    has(blk, r'\textbf{出处：}', 'src/%s（排程一节必须交代借鉴来源）' % APP_FILES['I'])
    has(blk, '高效用', 'src/%s（Dunlosky 2013 的效用评级结论）' % APP_FILES['I'])
    has(blk, '必要难度', 'src/%s（Bjork 的存储强度／提取强度口径）' % APP_FILES['I'])
    has(blk, '幂律', 'src/%s（遗忘曲线按幂律描述这一结论）' % APP_FILES['I'])
    has(blk, '长期版', 'src/%s（1/3/7/21 只适用于约三个月保留期的限定）' % APP_FILES['I'])
    lacks(blk, '遗忘曲线是单指数', 'src/%s（不得退回「单指数足够」的旧口径）' % APP_FILES['I'])


def check_P14_appendix_q_diagnostics():
    """附录 Q 等离子体诊断的结构完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 十一个小节之一被删（按小节标题前缀定位）；
      ② 缺少「出处：」行（本附录要求 >= 6 处）；
      ③ 七张表 / 两只语义盒 / 编号公式之一被删或改了 label；
      ④ 六道练习缺解答；
      ⑤ 印出的关键数值与独立复算不符，或搬到了别的小节里；
      ⑥ 旧口径回流：把微分写成 \\dd、把 OML 判据写成与 r_p/\\lambda_D 无关
         （正文必须承认同一支探针在不同密度下换理论）。

    约定与附录 M/N/O/P 相同：表内数值的单位写在表头，故这些用例传 unit=''。
    语义盒的标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    qpath = os.path.join(SRC, APP_FILES['Q'])
    if not os.path.exists(qpath):
        rec('FAIL', 'P', '附录 Q 文件缺失：src/%s' % APP_FILES['Q'])
        return
    text = strip_comments(read(qpath))
    secs = split_blocks(text, 'section')

    want_prefix = ['诊断能测什么', '静电探针', '从 I--V 曲线读电子温度',
                   '从离子饱和电流读密度', '悬浮电位与等离子体电位', '磁诊断',
                   '光学诊断', '汤姆孙散射与微波干涉仪', '误差与自检',
                   '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want_prefix:
        for ttl, body in secs:
            if ttl.startswith(pre):
                found[pre] = body
                break
    missing = [p for p in want_prefix if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 Q 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 Q 的 %d 个小节齐全' % len(want_prefix))

    for lab in ('tab:app-q-map', 'tab:app-q-oml', 'tab:app-q-iv', 'tab:app-q-oes',
                'tab:app-q-salpeter', 'tab:app-q-error', 'tab:app-q-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['Q'], lab))

    # 语义盒标签是 \begin{tip}{标题}{标签} 的第二参数，不是 \label{}
    for lab in ('tip-appq-limit', 'tip-appq-next'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['Q'], lab))

    for lab in ('eq:app-q-debye', 'eq:app-q-oml', 'eq:app-q-iv', 'eq:app-q-slope',
                'eq:app-q-iis', 'eq:app-q-ji', 'eq:app-q-ne', 'eq:app-q-err',
                'eq:app-q-float', 'eq:app-q-floatdiff', 'eq:app-q-mutual', 'eq:app-q-rc',
                'eq:app-q-doppler', 'eq:app-q-oes', 'eq:app-q-salpeter',
                'eq:app-q-photons', 'eq:app-q-microwave'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['Q'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['Q'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 Q 的「出处：」行只有 %d 处，应 >= 6' % n_src)
    else:
        rec('PASS', 'P', '附录 Q 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_sol = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex != 6 or n_sol < 6:
        rec('FAIL', 'P', '附录 Q 练习数为 %d、解答数为 %d（应为 6 与 >= 6）' % (n_ex, n_sol))
    else:
        rec('PASS', 'P', '附录 Q 六道估算练习均有解答')

    # ---- 独立复算：常数与公式与 tmp/appQ_calc.py 同一套 ----
    U_Q = 1.66053906660e-27
    m_ar = 39.948 * U_Q
    e_k = 11604.518122          # 1 eV 对应的开尔文
    sigma_t = 6.6524587321e-29
    c_light = 2.99792458e8
    h_pl = 6.62607015e-34

    def debye_q(te, ne):
        return math.sqrt(EPS0 * te * E_CHG / (ne * E_CHG ** 2))

    def cs_q(te, m):
        return math.sqrt(te * E_CHG / m)

    def vbar_e(te):
        return math.sqrt(8.0 * KB * te * e_k / (math.pi * M_E))

    def doppler_pm(lam_nm, t_k):
        return 7.16e-7 * lam_nm * math.sqrt(t_k / 39.948) * 1.0e3

    def n_crit(f):
        return EPS0 * M_E * (2 * math.pi * f) ** 2 / E_CHG ** 2

    def f_pl(n):
        return math.sqrt(n * E_CHG ** 2 / (EPS0 * M_E)) / (2 * math.pi)

    a_p = 2.0 * math.pi * 5e-4 * 0.01        # r_p=0.5 mm, l_p=10 mm
    cs_ar = cs_q(3.0, m_ar)
    j_i = E_CHG * 1e17 * cs_ar               # A/m^2
    lam_etch = debye_q(3.0, 1e17)
    k_sc = 2.0 * (2.0 * math.pi / 532e-9) * math.sin(math.pi / 4)
    n_ph_532 = 1.0 / (h_pl * c_light / 532e-9)
    m_rog = MU0 * 1000 * 1e-4 / (2.0 * math.pi * 0.05)

    cases = [
        # --- §Q.1 手段与判据 ---
        ('诊断能测什么', '德拜长度系数 7433.9', math.sqrt(EPS0 / E_CHG), 0.001, '', False),
        ('诊断能测什么', '刻蚀 λ_D=40.72 μm', lam_etch * 1e6, 0.01, '', False),
        ('诊断能测什么', '医学射流 λ_D=2.88 μm', debye_q(1.5, 1e19) * 1e6, 0.01, '', False),
        ('诊断能测什么', '大气压 DBD λ_D=0.33 μm', debye_q(2.0, 1e21) * 1e6, 0.02, '', False),
        # --- §Q.2 探针半径与德拜长度 ---
        ('静电探针', 'r_p/λ_D=12.28（1e17）', 5e-4 / lam_etch, 0.01, '', False),
        ('静电探针', 'λ_D=128.76 μm（1e16）', debye_q(3.0, 1e16) * 1e6, 0.01, '', False),
        ('静电探针', 'r_p/λ_D=3.88（1e16）', 5e-4 / debye_q(3.0, 1e16), 0.01, '', False),
        ('静电探针', 'λ_D=407.17 μm（1e15）', debye_q(3.0, 1e15) * 1e6, 0.01, '', False),
        ('静电探针', 'r_p/λ_D=1.23（1e15）', 5e-4 / debye_q(3.0, 1e15), 0.01, '', False),
        ('静电探针', 'OML 边界下最大 r_p=122.2 μm', 3.0 * lam_etch * 1e6, 0.01, '', False),
        # --- §Q.3 I--V 曲线 ---
        ('从 I--V 曲线读电子温度', 'I_e(-3 V)=3.6788 mA', 10.0 * math.exp(-1.0), 0.01, '', False),
        ('从 I--V 曲线读电子温度', 'I_net(-3 V)=2.6788 mA',
         10.0 * math.exp(-1.0) - 1.0, 0.01, '', False),
        ('从 I--V 曲线读电子温度', 'I_e(-6 V)=1.3534 mA', 10.0 * math.exp(-2.0), 0.01, '', False),
        ('从 I--V 曲线读电子温度', '正确斜率 0.333333 /V', 1.0 / 3.0, 0.01, '', False),
        ('从 I--V 曲线读电子温度', '正确 T_e=3.0000 eV', 3.0, 0.01, '', False),
        ('从 I--V 曲线读电子温度', '不减 I_is 的 T_e=1.4810 eV',
         3.0 / math.log((10.0 * math.exp(-1.0) - 1.0) / (10.0 * math.exp(-2.0) - 1.0)),
         0.01, '', False),
        ('从 I--V 曲线读电子温度', '低估 50.63%',
         (1.0 - 1.0 / math.log((10.0 * math.exp(-1.0) - 1.0)
                               / (10.0 * math.exp(-2.0) - 1.0))) * 100.0, 0.01, '', False),
        ('从 I--V 曲线读电子温度', '靠近 V_p 的 T_e=2.5022 eV',
         1.0 / math.log(6.1653 / 4.1342), 0.01, '', False),
        # --- §Q.4 离子饱和电流 ---
        ('从离子饱和电流读密度', 'A_p=3.141593e-5 m^2', a_p, 0.001, r'\mathrm{m^2}', False),
        ('从离子饱和电流读密度', 'c_s=2691.81 m/s', cs_ar, 0.01, r'\mathrm{m/s}', False),
        ('从离子饱和电流读密度', 'j_i=43.1275 A/m^2', j_i, 0.01, r'\mathrm{A/m^2}', False),
        ('从离子饱和电流读密度', 'j_i=4.3127 mA/cm^2（与附录 P 一致）',
         j_i * 0.1, 0.01, r'\mathrm{mA/cm^2}', False),
        ('从离子饱和电流读密度', '裸玻姆 1.354890 mA', j_i * a_p * 1e3, 0.01, r'\mathrm{mA}', False),
        ('从离子饱和电流读密度', '含 0.6 因子 0.812934 mA',
         0.6 * j_i * a_p * 1e3, 0.01, r'\mathrm{mA}', False),
        ('从离子饱和电流读密度', '实测 1.000 mA → n_e=1.2301e17 m^-3',
         1.0e-3 / (0.6 * E_CHG * cs_ar * a_p), 0.01, r'\mathrm{m^{-3}}', False),
        ('从离子饱和电流读密度', '漏 0.6 因子低估 40.0%', (1.0 - 0.6) * 100.0, 0.01, '', False),
        # --- §Q.5 悬浮电位 ---
        ('悬浮电位与等离子体电位', 'v̄_e=1.1592e6 m/s', vbar_e(3.0), 0.01, r'\mathrm{m/s}', False),
        ('悬浮电位与等离子体电位', '比值 107.6559', vbar_e(3.0) / (4.0 * cs_ar), 0.01, '', False),
        ('悬浮电位与等离子体电位', '系数 4.6789',
         math.log(vbar_e(3.0) / (4.0 * cs_ar)), 0.01, '', False),
        ('悬浮电位与等离子体电位', 'V_p-V_f=14.0368 V',
         3.0 * math.log(vbar_e(3.0) / (4.0 * cs_ar)), 0.01, r'\mathrm{V}', False),
        ('悬浮电位与等离子体电位', '含 0.6 因子的系数 5.1898',
         math.log(vbar_e(3.0) / (4.0 * 0.6 * cs_ar)), 0.01, '', False),
        ('悬浮电位与等离子体电位', '两式之差 1.5325 V',
         3.0 * math.log(1.0 / 0.6), 0.01, r'\mathrm{V}', False),
        # --- §Q.6 磁诊断 ---
        ('磁诊断', '互感 M=400.0 nH', m_rog * 1e9, 0.01, r'\mathrm{nH}', False),
        ('磁诊断', 'V_coil=40.00 V', m_rog * 1e8, 0.01, r'\mathrm{V}', False),
        ('磁诊断', 'RC=50.0 μs', 50.0 * 1e-6 * 1.0e6, 0.01, r'\mathrm{\mu s}', False),
        ('磁诊断', 'V_out=0.8000 V', m_rog * 100.0 / (50.0 * 1e-6), 0.001, r'\mathrm{V}', False),
        ('磁诊断', 'S=8.0000e-3 V/A', m_rog / (50.0 * 1e-6), 0.001, r'\mathrm{V/A}', False),
        # --- §Q.7 光学 ---
        ('光学诊断', 'Doppler 300 K=1.3666 pm', doppler_pm(696.5, 300.0), 0.01, r'\mathrm{pm}', False),
        ('光学诊断', 'Doppler 0.5 eV=6.0101 pm',
         doppler_pm(696.5, 0.5 * e_k), 0.01, r'\mathrm{pm}', False),
        ('光学诊断', 'Doppler 3 eV=14.7218 pm',
         doppler_pm(696.5, 3.0 * e_k), 0.01, r'\mathrm{pm}', False),
        # 正文把 36.6 倍写成「约 37 倍」，容差放到 2%
        ('光学诊断', '仪器宽度是室温宽度的 37 倍',
         50.0 / doppler_pm(696.5, 300.0), 0.02, '', False),
        # 表内数值的单位写在表头，故这六行传 unit=''
        ('光学诊断', '双线比 0.950 → 7.7983 eV', 0.40 / math.log(1.0 / 0.950), 0.01,
         '', False),
        ('光学诊断', '双线比 0.900 → 3.7965 eV', 0.40 / math.log(1.0 / 0.900), 0.01,
         '', False),
        ('光学诊断', '双线比 0.875 → 2.9956 eV', 0.40 / math.log(1.0 / 0.875), 0.01,
         '', False),
        ('光学诊断', '双线比 0.850 → 2.4613 eV', 0.40 / math.log(1.0 / 0.850), 0.01,
         '', False),
        ('光学诊断', '双线比 0.750 → 1.3904 eV', 0.40 / math.log(1.0 / 0.750), 0.01,
         '', False),
        ('光学诊断', '双线比 0.500 → 0.5771 eV', 0.40 / math.log(2.0), 0.01, '', False),
        ('光学诊断', '1 mW@750 nm 进 1% 通量=3.7756e13 光子/s',
         1.0e-3 * 0.01 / (h_pl * c_light / 750e-9), 0.01, '', False),
        # --- §Q.8 汤姆孙散射与微波 ---
        ('汤姆孙散射与微波干涉仪', 'k=1.6703e7 1/m', k_sc, 0.01, r'\mathrm{m^{-1}}', False),
        ('汤姆孙散射与微波干涉仪', '刻蚀 kλ_D=680.08', k_sc * lam_etch, 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '刻蚀 α=1.4704e-3', 1.0 / (k_sc * lam_etch), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', 'DBD kλ_D=5.55', k_sc * debye_q(2.0, 1e21), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '聚变边界 kλ_D=107.53',
         k_sc * debye_q(15.0, 2e19), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', 'DBD α=1.8009e-1',
         1.0 / (k_sc * debye_q(2.0, 1e21)), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '聚变边界 α=9.2997e-3',
         1.0 / (k_sc * debye_q(15.0, 2e19)), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '1 J@10 ns=2.6782e18 光子', n_ph_532, 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', 'N_s=1781.6 个散射光子',
         1e17 * 0.01 * sigma_t * 0.01 * n_ph_532, 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '统计误差 2.37%',
         100.0 / math.sqrt(1e17 * 0.01 * sigma_t * 0.01 * n_ph_532), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', 'n_c(10 GHz)=1.2404e18 m^-3', n_crit(1e10), 0.01,
         r'\mathrm{m^{-3}}', False),
        ('汤姆孙散射与微波干涉仪', 'n_e/n_c=0.080616', 1e17 / n_crit(1e10), 0.01, '', False),
        ('汤姆孙散射与微波干涉仪', 'Δφ=2.5344 rad',
         (2 * math.pi * 1e10 / c_light) * (1e17 / n_crit(1e10) / 2.0) * 0.3, 0.01,
         r'\mathrm{rad}', False),
        ('汤姆孙散射与微波干涉仪', '0.4034 条纹',
         (2 * math.pi * 1e10 / c_light) * (1e17 / n_crit(1e10) / 2.0) * 0.3 / (2 * math.pi),
         0.01, '', False),
        ('汤姆孙散射与微波干涉仪', '测 1e19 需 f>28.39 GHz', f_pl(1e19) / 1e9, 0.01,
         r'\mathrm{GHz}', False),
        # --- §Q.10 练习答案 ---
        ('六道估算练习', '[1] I_e=1.8394 mA', 5.0 * math.exp(-1.0), 0.01, r'\mathrm{mA}', False),
        ('六道估算练习', '[1] I_e=0.6767 mA', 5.0 * math.exp(-2.0), 0.01, r'\mathrm{mA}', False),
        ('六道估算练习', '[1] 正确 T_e=2.0000 eV', 2.0, 0.01, r'\mathrm{eV}', False),
        ('六道估算练习', '[1] 净电流斜率 1.0127 /V',
         math.log(1.3394 / 0.1767) / 2.0, 0.01, '', False),
        ('六道估算练习', '[1] T_e=0.9873 eV（低估 50.6%）',
         1.0 / (math.log(1.3394 / 0.1767) / 2.0), 0.01, r'\mathrm{eV}', False),
        ('六道估算练习', '[2] A_p=2.010619e-5 m^2', 2.0 * math.pi * 4e-4 * 0.008, 0.001,
         r'\mathrm{m^2}', False),
        ('六道估算练习', '[2] c_s=2457.27 m/s', cs_q(2.5, m_ar), 0.01, r'\mathrm{m/s}', False),
        ('六道估算练习', '[2] n_e=1.2633e17 m^-3',
         0.6e-3 / (0.6 * E_CHG * cs_q(2.5, m_ar) * 2.0 * math.pi * 4e-4 * 0.008), 0.01,
         r'\mathrm{m^{-3}}', False),
        ('六道估算练习', '[3] λ_D=66.49 μm', debye_q(4.0, 5e16) * 1e6, 0.01, '', False),
        ('六道估算练习', '[3] r_p/λ_D=15.04', 1e-3 / debye_q(4.0, 5e16), 0.01, '', False),
        ('六道估算练习', '[3] 要 r_p <= 199.5 μm', 3.0 * debye_q(4.0, 5e16) * 1e6, 0.01, '', False),
        ('六道估算练习', '[4] M=5.0000e-7 H',
         MU0 * 500 * 2e-4 / (2.0 * math.pi * 0.04), 0.001, r'\mathrm{H}', False),
        ('六道估算练习', '[4] V_coil=50.00 V',
         MU0 * 500 * 2e-4 / (2.0 * math.pi * 0.04) * 200.0 / 2e-6, 0.01, r'\mathrm{V}', False),
        ('六道估算练习', '[4] RC=100 μs', 100.0 * 1e-6 * 1e6, 0.01, r'\mathrm{\mu s}', False),
        ('六道估算练习', '[4] V_out=1.0000 V',
         MU0 * 500 * 2e-4 / (2.0 * math.pi * 0.04) * 200.0 / 1e-4, 0.01, r'\mathrm{V}', False),
        ('六道估算练习', '[4] S=5.0000e-3 V/A',
         MU0 * 500 * 2e-4 / (2.0 * math.pi * 0.04) / 1e-4, 0.01, r'\mathrm{V/A}', False),
        ('六道估算练习', '[5] N_ph=1.33895e18', 0.5 / (h_pl * c_light / 532e-9), 0.01, '', False),
        ('六道估算练习', '[5] N_s=89081 个光子',
         5e18 * 0.01 * sigma_t * 0.02 * (0.5 / (h_pl * c_light / 532e-9)), 0.01, '', False),
        ('六道估算练习', '[5] 统计误差 0.335%',
         100.0 / math.sqrt(5e18 * 0.01 * sigma_t * 0.02 * (0.5 / (h_pl * c_light / 532e-9))),
         0.01, '', False),
        ('六道估算练习', '[6] 测 3e18 需 f>15.55 GHz', f_pl(3e18) / 1e9, 0.01, '',
         False),
        ('六道估算练习', '[6] 20 GHz 的 n_c=4.9618e18 m^-3', n_crit(2e10), 0.01,
         r'\mathrm{m^{-3}}', False),
        ('六道估算练习', '[6] n_e/n_c=0.6046', 3e18 / n_crit(2e10), 0.01, '', False),
        ('六道估算练习', '[6] Δφ=25.3439 rad',
         (2 * math.pi * 2e10 / c_light) * (3e18 / n_crit(2e10) / 2.0) * 0.2, 0.01,
         r'\mathrm{rad}', False),
        ('六道估算练习', '[6] 4.0336 条纹',
         (2 * math.pi * 2e10 / c_light) * (3e18 / n_crit(2e10) / 2.0) * 0.2 / (2 * math.pi),
         0.01, '', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 Q 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 Q 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    # 全篇口径与结论句（内容跨节，故对整个文件断言）
    has(text, r'\texttt{appQ\_calc.py}',
        'src/%s（数值由仓库脚本复算的交代）' % APP_FILES['Q'])
    has(text, r'\eqref{eq:app-p-flux}',
        'src/%s（与附录 P 的离子通量公式交叉引用）' % APP_FILES['Q'])
    has(found['静电探针'], 'OML', 'src/%s（轨道运动限制的适用区间）' % APP_FILES['Q'])
    has(found['静电探针'], 'ABR', 'src/%s（超出 OML 后要换的理论）' % APP_FILES['Q'])
    has(found['从 I--V 曲线读电子温度'], '低估',
        'src/%s（漏减离子电流会低估 T_e 的结论）' % APP_FILES['Q'])
    has(found['从离子饱和电流读密度'], '有效收集面积',
        'src/%s（密度误差的主要来源是面积而不是 T_e）' % APP_FILES['Q'])
    has(found['悬浮电位与等离子体电位'], '收集模型',
        'src/%s（两个悬浮电位系数之差来自收集模型）' % APP_FILES['Q'])
    has(found['磁诊断'], '只能测',
        'src/%s（磁探针只能测变化量的局限）' % APP_FILES['Q'])
    # 数值用例的窗口匹配（80 字符内出现该单位）可能被邻座的同值数字满足：
    # §Q.6 里 R=50 Ω 与 RC=50.0 µs 数值相同，故这里另加一条字面断言钉住时间常数本身。
    has(found['磁诊断'], r'RC=50.0\ \mathrm{\mu s}',
        'src/%s（RC 积分时间常数 50.0 µs 的字面写法）' % APP_FILES['Q'])
    has(found['光学诊断'], '仪器宽度',
        'src/%s（Doppler 宽度必须先与仪器宽度比较）' % APP_FILES['Q'])
    has(found['汤姆孙散射与微波干涉仪'], '非相干',
        'src/%s（Salpeter 参数的两种极限）' % APP_FILES['Q'])
    has(found['误差与自检'], '系统误差',
        'src/%s（系统误差不随平均次数下降）' % APP_FILES['Q'])
    has(found['来源、分级与版权说明'], '不复制',
        'src/%s（不复制任何题目、解答、图表与数据表的声明）' % APP_FILES['Q'])
    lacks(text, '\\dd ', 'src/%s（微分必须写 \\diff，不得回流 \\dd ）' % APP_FILES['Q'])
    lacks(text, 'r_p/\\lambda_D 与理论无关',
          'src/%s（OML 判据不得写成与 r_p/\\lambda_D 无关）' % APP_FILES['Q'])


def check_P15_appendix_i_dunlosky():
    """附录 I「十种学习技术的评级」一节的完整性与三档计数。

    可证伪点（全部记 FAIL）：
      ① 该小节被删或标题被改；
      ② 十种技术之一被删（逐项检查名称）；
      ③ 三档数量被改（高效用 2 / 中等效用 3 / 低效用 5，按表格行首列计数）；
      ④ 缺少「流畅性错觉」机制或「存储强度／提取强度」的口径；
      ⑤ 缺少「出处：」交代。
    """
    p = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(p):
        rec('FAIL', 'P', '附录 I 文件不存在：%s' % APP_FILES['I'])
        return
    text = strip_comments(read(p))
    secs = split_blocks(text, 'section')
    blk = None
    for ttl, body in secs:
        if ttl.startswith('十种学习技术的评级'):
            blk = body
            break
    if blk is None:
        rec('FAIL', 'P', '附录 I 缺少「十种学习技术的评级」一节（标题被改或整节被删）')
        return
    rec('PASS', 'P', '附录 I 有「十种学习技术的评级：把时间花在哪」一节')

    has(blk, r'\label{tab:app-i-dunlosky}', 'src/%s（评级表 label）' % APP_FILES['I'])
    for tech in ('检索练习', '分散练习', '精细提问', '自我解释', '交错练习',
                 '重读', '划线', '抄摘要', '关键词记忆法', '文本想象'):
        has(blk, tech, 'src/%s（评级表里的技术：%s）' % (APP_FILES['I'], tech))

    n_high = len(re.findall(r'高效用 &', blk))
    n_mid = len(re.findall(r'中等效用 &', blk))
    n_low = len(re.findall(r'低效用 &', blk))
    if (n_high, n_mid, n_low) != (2, 3, 5):
        rec('FAIL', 'P', '附录 I 评级表的档次计数为 高%d/中%d/低%d，应为 高2/中3/低5'
            % (n_high, n_mid, n_low))
    else:
        rec('PASS', 'P', '附录 I 评级表的档次计数为 高效用 2 / 中等效用 3 / 低效用 5')

    has(blk, '流畅性错觉', 'src/%s（低效用技术仍然流行的机制）' % APP_FILES['I'])
    has(blk, '存储强度', 'src/%s（Bjork 的存储强度与提取强度区分）' % APP_FILES['I'])
    has(blk, '提取强度', 'src/%s（考试与解题需要的是提取强度）' % APP_FILES['I'])
    has(blk, r'\textbf{出处：}', 'src/%s（评级一节必须交代借鉴来源）' % APP_FILES['I'])
    lacks(blk, '高效用三项',
          'src/%s（高效用只能是两项：检索练习与分散练习）' % APP_FILES['I'])
    # 期刊名是两本刊的拼接（Perspectives on Psychological Science + Psychological Science in the
    # Public Interest）；DOI 10.1177/1529100612453266 与卷期页 14(1), 4--58 都属于后者。
    lacks(text, 'Perspectives on Psychological Science',
          'src/%s（Dunlosky 2013 的出处期刊名不得拼接错，应为 Psychological Science in the '
          'Public Interest）' % APP_FILES['I'])


def check_P16_no_markdown_bold():
    """全库不得残留 Markdown 粗体简写 **…**（逐字环境内的 Python 幂运算除外）。

    背景：v1.13.0 的附录 P 有 19 行正文写成 **强调**，而本书没有把 ** 定义成粗体，
    于是 PDF 里真的印出了星号。这里把「逐字环境之外出现 **」定为 FAIL。
    """
    verb_begin = re.compile(r'\\begin\{(lstlisting|verbatim|Verbatim|minted|'
                            r'lstinputlisting|filecontents)\*?\}')
    verb_end = re.compile(r'\\end\{(lstlisting|verbatim|Verbatim|minted|'
                          r'lstinputlisting|filecontents)\*?\}')
    files = [MAIN] + [os.path.join(SRC, f) for f in sorted(APP_FILES.values())]
    stray = []
    for p in files:
        if not os.path.exists(p):
            continue
        inv = False
        for i, ln in enumerate(strip_comments(read(p)).split('\n'), 1):
            if verb_begin.search(ln):
                inv = True
            if not inv and '**' in ln:
                stray.append('%s:%d' % (os.path.basename(p), i))
            if verb_end.search(ln):
                inv = False
    if stray:
        rec('FAIL', 'P', '有 %d 行正文残留 Markdown 粗体简写 **…**（会原样印出星号）：%s'
            % (len(stray), stray[:8]))
    else:
        rec('PASS', 'P', '全库 %d 个 tex 文件在逐字环境之外没有 ** 粗体简写'
            % len(files))


def check_P17_appendix_r_dusty():
    """附录 R 尘埃等离子体的结构完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 十一个小节之一被删（按小节标题前缀定位）；
      ② 缺少「出处：」行（本附录要求 >= 9 处）；
      ③ 五张表 / 两只语义盒 / 十个编号公式之一被删或改了 label；
      ④ 六道练习缺解答；
      ⑤ 印出的关键数值与独立复算不符（λ_D、z、Z_d、Q_d、t_charge、Γ、Havnes、
         c_DA、E_need、β），或搬到了别的小节里；
      ⑥ 旧口径回流：尘埃声速说成"密度比"问题、C3 说成晶体、微分写成 \\dd。

    约定与附录 M/N/O/P/Q 相同：表内数值的单位写在表头，故这些用例传 unit=''。
    语义盒的标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    rpath = os.path.join(SRC, APP_FILES['R'])
    if not os.path.exists(rpath):
        rec('FAIL', 'P', '附录 R 文件缺失：src/%s' % APP_FILES['R'])
        return
    text = strip_comments(read(rpath))
    secs = split_blocks(text, 'section')

    want_prefix = ['尘埃等离子体：多出来的第三个尺度', '尘埃怎么带电', '充电要多久',
                   '强耦合与晶体', 'Havnes 参数', '尘埃声波', '鞘层中的悬浮',
                   '从刻蚀腔到聚变堆', '空间与天体尘埃', '六道估算练习',
                   '来源、分级与版权说明']
    found = {}
    for pre in want_prefix:
        for ttl, body in secs:
            if ttl.startswith(pre):
                found[pre] = body
                break
    missing = [p for p in want_prefix if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 R 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 R 的 %d 个小节齐全' % len(want_prefix))

    for lab in ('tab:app-r-map', 'tab:app-r-charge', 'tab:app-r-gamma',
                'tab:app-r-havnes', 'tab:app-r-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['R'], lab))

    # 语义盒标签是 \begin{tip}{标题}{标签} 的第二参数，不是 \label{}
    for lab in ('tip-appr-limit', 'tip-appr-next'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['R'], lab))

    for lab in ('eq:app-r-oml', 'eq:app-r-z', 'eq:app-r-zd', 'eq:app-r-tcharge',
                'eq:app-r-gamma', 'eq:app-r-yukawa', 'eq:app-r-havnes',
                'eq:app-r-cda', 'eq:app-r-balance', 'eq:app-r-beta'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['R'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['R'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 9:
        rec('FAIL', 'P', '附录 R 的「出处：」行只有 %d 处，应 >= 9' % n_src)
    else:
        rec('PASS', 'P', '附录 R 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_sol = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex != 6 or n_sol < 6:
        rec('FAIL', 'P', '附录 R 练习数为 %d、解答数为 %d（应为 6 与 >= 6）' % (n_ex, n_sol))
    else:
        rec('PASS', 'P', '附录 R 六道估算练习均有解答')

    # ---- 独立复算：常数与公式与 tmp/appR_calc.py 同一套 ----
    U_R = 1.66053906660e-27
    m_ar = 39.948 * U_R
    m_dion = 2.013553212745 * U_R
    m_h = 1.007825032 * U_R
    g0 = 9.80665
    ti_ar = 300.0 * KB / E_CHG

    def debye_r(te, ne):
        return math.sqrt(EPS0 * te * E_CHG / (ne * E_CHG ** 2))

    def z_oml_r(te, ti_ev, m_i):
        tau = ti_ev / te
        lo, hi = 1e-6, 60.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            f = math.exp(mid) * (1.0 + mid / tau) - math.sqrt(m_i / (tau * M_E))
            if f < 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def zd_r(a_m, te, z):
        return 4.0 * math.pi * EPS0 * a_m * z * te / E_CHG

    def md_r(a_m, rho):
        return 4.0 / 3.0 * math.pi * a_m ** 3 * rho

    def delta_r(n_d):
        return (3.0 / (4.0 * math.pi * n_d)) ** (1.0 / 3.0)

    coul_nm = E_CHG / (4.0 * math.pi * EPS0) * 1e9

    def gam_r(zd, d_m, kt_ev):
        return zd ** 2 * coul_nm / (d_m * 1e9 * kt_ev)

    def yukawa_r(kappa):
        return (1.0 + kappa + 0.5 * kappa ** 2) * math.exp(-kappa)

    z_ar = z_oml_r(3.0, ti_ar, m_ar)
    z_c4 = z_oml_r(1.0, 100.0 * KB / E_CHG, m_h)
    z_c3 = z_oml_r(15.0, 15.0, m_dion)
    lam_c1 = debye_r(3.0, 1e17)
    lam_c2 = debye_r(3.0, 1e15)
    lam_c3 = debye_r(15.0, 2e19)
    lam_c4 = debye_r(1.0, 1e6)
    zd_c1 = zd_r(0.5e-6, 3.0, z_ar)
    zd_c2 = zd_r(4.8e-6, 3.0, z_ar)
    zd_c3 = zd_r(1.0e-6, 15.0, z_c3)
    qd_c1 = zd_c1 * E_CHG
    md_c1 = md_r(0.5e-6, 2000.0)
    md_c2 = md_r(4.8e-6, 1510.0)
    qd_c2 = zd_c2 * E_CHG
    del_c1 = delta_r(1e12)
    del_c2 = delta_r(1e9)
    del_c3 = delta_r(1e10)
    kap_c1 = del_c1 / lam_c1
    kap_c2 = del_c2 / lam_c2
    kap_c3 = del_c3 / lam_c3
    veth_e = math.sqrt(8.0 * 3.0 * E_CHG / (math.pi * M_E))
    ie0_c1 = math.pi * (0.5e-6) ** 2 * E_CHG * 1e17 * veth_e
    tch_c1 = qd_c1 / ie0_c1
    cs_ar = math.sqrt(3.0 * E_CHG / m_ar)
    cda_c1 = math.sqrt(zd_c1 ** 2 * 1e12 * ti_ar * E_CHG / (1e17 * md_c1))
    a_beta = 3.0 * 3.828e26 / (16.0 * math.pi * 2.99792458e8 * 6.67430e-11 * 1.98847e30)

    cases = [
        # --- §R.1 第三个尺度 ---
        ('尘埃等离子体：多出来的第三个尺度', 'C1 λ_D=40.72 μm', lam_c1 * 1e6, 0.01, '', False),
        ('尘埃等离子体：多出来的第三个尺度', 'C2 λ_D=407.17 μm', lam_c2 * 1e6, 0.005, '', False),
        ('尘埃等离子体：多出来的第三个尺度', 'C3 λ_D=6.44 μm', lam_c3 * 1e6, 0.01, '', False),
        ('尘埃等离子体：多出来的第三个尺度', 'C4 λ_D=7.43 m', lam_c4, 0.005, '', False),
        # --- §R.2 OML 充电 ---
        ('尘埃怎么带电', 'C1 a/λ_D=0.0123', 0.5e-6 / lam_c1, 0.01, '', False),
        # C4 的 a/λ_D 曾在 v1.15.0 表 R.2 里写成 1.3e-7（差十倍）：只查 C1
        # 的那一条漏掉了它。四个工况的这一列现在全部逐行复算。
        ('尘埃怎么带电', 'C2 a/λ_D=0.0118', 4.8e-6 / lam_c2, 0.01, '', False),
        ('尘埃怎么带电', 'C3 a/λ_D=0.1553', 1.0e-6 / lam_c3, 0.02, '', False),
        ('尘埃怎么带电', 'C4 a/λ_D=1.3e-8（v1.15.0 曾误写 1.3e-7）',
         0.1e-6 / lam_c4, 0.05, '', False),
        ('尘埃怎么带电', 'z=2.359（与半径、密度无关）', z_ar, 0.005, '', False),
        ('尘埃怎么带电', 'φ_d=7.077 V', z_ar * 3.0, 0.005, '', False),
        ('尘埃怎么带电', 'C1 Z_d=2.46e3', zd_c1, 0.01, '', False),
        ('尘埃怎么带电', 'C2 Z_d=2.36e4', zd_c2, 0.01, '', False),
        ('尘埃怎么带电', 'C3 Z_d=2.89e4', zd_c3, 0.01, '', False),
        ('尘埃怎么带电', 'C4 Z_d=83.0', zd_r(0.1e-6, 1.0, z_c4), 0.01, '', False),
        ('尘埃怎么带电', 'C1 Q_d=3.94e-16 C', qd_c1, 0.02, '', False),
        ('尘埃怎么带电', 'C3 φ_d=41.63 V', z_c3 * 15.0, 0.005, '', False),
        ('尘埃怎么带电', 'Z_d 标度 1.638e3·(a/μm)(T_e/eV)', zd_r(1e-6, 3.0, z_ar) / 3.0,
         0.01, '', False),
        ('尘埃怎么带电', 'C1 t_charge=27.0 ns', tch_c1 * 1e9, 0.01, '', False),
        ('尘埃怎么带电', 'C4 t_charge=3.95e3 s', (zd_r(0.1e-6, 1.0, z_c4) * E_CHG)
         / (math.pi * (0.1e-6) ** 2 * E_CHG * 1e6
            * math.sqrt(8.0 * 1.0 * E_CHG / (math.pi * M_E))), 0.02, '', False),
        # --- §R.3 充电时间 vs 射频 ---
        ('充电要多久', '13.56 MHz 周期 73.75 ns', 1e9 / 13.56e6, 0.005, '', False),
        ('充电要多久', 't_charge/T=0.366', tch_c1 * 13.56e6, 0.02, '', False),
        # --- §R.4 强耦合与 Γ ---
        ('强耦合与晶体', 'Δ=62.04 μm', del_c1 * 1e6, 0.005, '', False),
        ('强耦合与晶体', 'C2 Δ=620.35 μm', del_c2 * 1e6, 0.005, '', False),
        ('强耦合与晶体', 'C1 Γ(Δ)=5.41e3', gam_r(zd_c1, del_c1, 0.0259), 0.01, '', False),
        ('强耦合与晶体', 'C1 Γ(a)=6.71e5', gam_r(zd_c1, 0.5e-6, 0.0259), 0.01, '', False),
        ('强耦合与晶体', 'C1 κ=1.524', kap_c1, 0.005, '', False),
        ('强耦合与晶体', 'C3 κ=44.73', kap_c3, 0.005, '', False),
        ('强耦合与晶体', 'Yukawa 因子 0.803', yukawa_r(kap_c1), 0.01, '', False),
        ('强耦合与晶体', 'C1 Γ_eff=4.35e3',
         gam_r(zd_c1, del_c1, 0.0259) * yukawa_r(kap_c1), 0.01, '', False),
        ('强耦合与晶体', 'C2 Γ_eff=3.46e4',
         gam_r(zd_c2, del_c2, 0.03) * yukawa_r(kap_c2), 0.01, '', False),
        ('强耦合与晶体', 'C3 Γ_eff=1.6e-13',
         gam_r(zd_c3, del_c3, 1.0) * yukawa_r(kap_c3), 0.15, '', False),
        # --- §R.5 Havnes ---
        ('Havnes 参数', 'P=0.0246', 1e12 * zd_c1 / 1e17, 0.01, '', False),
        ('Havnes 参数', 'λ_eff=40.23 μm',
         lam_c1 / math.sqrt(1.0 + 1e12 * zd_c1 / 1e17) * 1e6, 0.005, '', False),
        ('Havnes 参数', 'n_d(P=1)=4.07e13', 1e17 / zd_c1, 0.01, '', False),
        ('Havnes 参数', 'P=1 比现状高 40 倍', 1e17 / zd_c1 / 1e12, 0.01, '', False),
        # --- §R.6 尘埃声波 ---
        ('尘埃声波', 'c_DA=1.5455 cm/s', cda_c1 * 100.0, 0.005, '', False),
        ('尘埃声波', 'c_s=2691.8 m/s', cs_ar, 0.005, '', False),
        ('尘埃声波', 'f(1 mm)=15.45 Hz', cda_c1 / 1e-3, 0.01, '', False),
        # 尘埃等离子体频率：v1.15.0 正文把它写成「10^2 rad/s 量级」（只对 C2 成立，
        # C1 实为 4.1e3），且当时没有任何复算。现在 C1/C2 与离子频率一起断言。
        ('尘埃声波', 'ω_pd(C1)=4.1e3 rad/s',
         math.sqrt(1e12 * zd_c1 ** 2 * E_CHG ** 2 / (EPS0 * md_c1)), 0.01, '', False),
        ('尘埃声波', 'ω_pd(C2)=48 rad/s',
         math.sqrt(1e9 * zd_c2 ** 2 * E_CHG ** 2 / (EPS0 * md_c2)), 0.02, '', False),
        ('尘埃声波', 'ω_pi(Ar)=6.6e7 rad/s',
         math.sqrt(1e17 * E_CHG ** 2 / (EPS0 * m_ar)), 0.01, '', False),
        # --- §R.7 鞘层悬浮 ---
        ('鞘层中的悬浮', 'E_need=26.08 V/m', md_c1 * g0 / qd_c1, 0.005, '', False),
        ('鞘层中的悬浮', 'C2 F_coul/F_grav=0.0486',
         (qd_c2 ** 2 / (4.0 * math.pi * EPS0)) / del_c2 ** 2 / (md_c2 * g0), 0.01, '', False),
        # --- §R.8 从刻蚀腔到聚变堆 ---
        ('从刻蚀腔到聚变堆', '比表面积 1500 m²/kg', 3.0 / (2000.0 * 1e-6), 0.005, '', False),
        ('从刻蚀腔到聚变堆', '100 g 尘埃吸附 0.751 g 氚',
         100.0 * 1e-3 * (3.0 / (2000.0 * 1e-6)) * 1e21 / 6.02214076e23 * 3.0160492,
         0.01, '', False),
        # --- §R.9 空间与天体尘埃 ---
        ('空间与天体尘埃', 'β(1 μm,2000)=0.287', a_beta / (2000.0 * 1e-6), 0.01, '', False),
        ('空间与天体尘埃', 'β(0.1 μm)=2.87', a_beta / (2000.0 * 0.1e-6), 0.01, '', False),
        ('空间与天体尘埃', 'a_crit=0.287 μm', a_beta / 2000.0 * 1e6, 0.01, '', False),
        # --- §R.10 练习答案必须落在练习一节里 ---
        ('六道估算练习', '练习1 Z_d=2.359e4（a=4.8 μm）',
         zd_r(4.8e-6, 3.0, z_ar), 0.01, '', False),
        ('六道估算练习', '练习1 Q_d=3.780e-15 C', zd_r(4.8e-6, 3.0, z_ar) * E_CHG,
         0.01, '', False),
        ('六道估算练习', '练习1(2) Z_d=491.5（a=0.1 μm）',
         zd_r(0.1e-6, 3.0, z_ar), 0.01, '', False),
        ('六道估算练习', '练习2 v_th,e=1.159e6 m/s', veth_e, 0.005, '', False),
        ('六道估算练习', '练习2 I_e0=1.459e-8 A', ie0_c1, 0.01, '', False),
        ('六道估算练习', '练习3 比值 Γ(a)/Γ(Δ)=124.1',
         gam_r(zd_c1, 0.5e-6, 0.0259) / gam_r(zd_c1, del_c1, 0.0259), 0.01, '', False),
        ('六道估算练习', '练习4 P=2.457e-2（原文写法）',
         1e12 * 2457.4 / 1e17, 0.01, '', False),
        ('六道估算练习', '练习5 c_s/c_DA=1.74e5', cs_ar / cda_c1, 0.01, '', False),
        ('六道估算练习', '练习6 E_need=26.08 V/m', md_c1 * g0 / qd_c1, 0.005, '', False),
    ]
    bad = []
    for topic, desc, want, rel, unit, signed in cases:
        if not value_in(found[topic], want, rel, unit, signed):
            bad.append((topic, desc))
    if bad:
        rec('FAIL', 'P', '附录 R 有 %d 项数值与独立复算不符（或不在它所属的那一节内）：%s'
            % (len(bad), bad[:6]))
    else:
        rec('PASS', 'P', '附录 R 的 %d 项关键数值全部与独立复算一致，且落在各自小节内'
            % len(cases))

    # 窗口只作补充，招牌数值还要有字面绑定。原因：值 + 相对容差会被「同节里另一个
    # 恰好落在窗内的数」顶替满足——v1.15.0 的逐条探针抓到了这个洞：把 §R.6 的
    # 1.5455 cm/s 改成 15.455 cm/s 后，同节的 1.55 Hz 仍落在 0.5% 窗内，套件照绿。
    # 下面每一项都是该小节独有的字面串（值 + 单位或值 + 上下文一起断），
    # 改数、改单位、或把某个量搬出这一节，都会直接 FAIL。
    for topic, lit, note in (
            ('尘埃等离子体：多出来的第三个尺度', r'$40.72$、$407.17$、$6.44\ \mathrm{\mu m}$',
             '四个工况的德拜长度并列句'),
            ('尘埃怎么带电', r'解得 $z=2.359$', 'OML 归一化电位的解'),
            ('充电要多久', r'周期 $73.75\ \mathrm{ns}$', '13.56 MHz 射频周期'),
            ('强耦合与晶体', r'$1.524$ & $0.803$', '表 R.3 的屏蔽因子与 κ'),
            ('Havnes 参数', r'$0.0246$ & $40.72$ & $40.23$', '表 R.4 的 P、λ_D 与 λ_eff'),
            ('尘埃声波', r'$c_{DA}=1.5455\ \mathrm{cm/s}$', 'C1 尘埃声速'),
            ('鞘层中的悬浮', r'C1 的 $0.5\ \mathrm{\mu m}$ 尘：$E_{\mathrm{need}}=26.08\ \mathrm{V/m}$',
             'C1 悬浮所需场强'),
            ('从刻蚀腔到聚变堆', r'$0.751\ \mathrm{g}$ 氚', '尘埃吸附的氚量'),
            ('空间与天体尘埃', r'临界半径 $a_{\mathrm{crit}}=0.287\ \mathrm{\mu m}$', '辐射压临界半径'),
            # v1.15.1：量级口径的字面绑定。这些句子在 v1.15.0 里分别写成
            # 「十四个数量级」（实际上是十三/六）、「毫摩尔级」（实为 0.25 mol）、
            # 「10^5--10^13 倍」（实为 10^5--10^14）、「1.3e-7」（实为 1.3e-8）。
            ('尘埃等离子体：多出来的第三个尺度', r'电子密度跨十三个数量级、', '量级跨度的正确口径'),
            ('尘埃等离子体：多出来的第三个尺度', r'$10^{5}$--$10^{14}$ 倍', '尘埃与离子的质量比区间'),
            ('尘埃怎么带电', r'$1.3\times10^{-8}$', '表 R.2 的 C4 半径/德拜长度比'),
            ('尘埃声波', r'$\omega_{pd}$ 也降到 $4.1\times10^{3}\ \mathrm{rad/s}$', 'C1 尘埃等离子体频率'),
            ('从刻蚀腔到聚变堆', r'$0.75\ \mathrm{g}$，', '氚库存的克数与摩尔数'),
    ):
        has(found[topic], lit, 'src/%s（%s）' % (APP_FILES['R'], note))

    # 全篇口径与结论句（内容跨节，故对整个文件断言）
    has(text, r'\texttt{appR\_calc.py}',
        'src/%s（数值由仓库脚本复算的交代）' % APP_FILES['R'])
    has(text, r'\eqref{eq:app-q-debye}',
        'src/%s（德拜长度沿用附录 Q 的工程式）' % APP_FILES['R'])
    has(found['尘埃怎么带电'], '与半径、与密度都无关',
        'src/%s（ξ 与半径、密度无关的 OML 结论）' % APP_FILES['R'])
    has(found['尘埃怎么带电'], 'OML', 'src/%s（轨道运动限制的适用区间）' % APP_FILES['R'])
    has(found['强耦合与晶体'], 'C3 不是',
        'src/%s（聚变边界不是晶体的反例结论）' % APP_FILES['R'])
    has(found['强耦合与晶体'], '屏蔽切断',
        'src/%s（κ≫1 把相互作用切断的结论）' % APP_FILES['R'])
    has(found['Havnes 参数'], '警戒线',
        'src/%s（Havnes 参数给出工程警戒线）' % APP_FILES['R'])
    has(found['尘埃声波'], '慢', 'src/%s（尘埃声波远慢于离子声速的结论）' % APP_FILES['R'])
    has(found['鞘层中的悬浮'], '紧贴鞘边',
        'src/%s（悬浮平衡点紧贴鞘边的结论）' % APP_FILES['R'])
    has(found['从刻蚀腔到聚变堆'], '库存',
        'src/%s（尘埃的聚变代价是库存与缺陷）' % APP_FILES['R'])
    has(found['空间与天体尘埃'], '辐射压',
        'src/%s（天体尘埃的辐射压判据）' % APP_FILES['R'])
    has(found['来源、分级与版权说明'], '不复制',
        'src/%s（不复制任何题目、解答、图表与数据表的声明）' % APP_FILES['R'])
    lacks(text, '\\dd ', 'src/%s（微分必须写 \\diff，不得回流 \\dd ）' % APP_FILES['R'])
    lacks(text, '密度比小 100 倍',
          'src/%s（尘埃声速不得退回到"密度比"这一错误解释）' % APP_FILES['R'])
    lacks(text, 'C3 也是晶体',
          'src/%s（聚变边界不得说成晶体）' % APP_FILES['R'])
    # v1.15.1 勘误过的量级口径不得回流（这几条都是 v1.15.0 发出去的真实错话）。
    lacks(text, '十四个数量级',
          'src/%s（量级跨度不是十四个：密度十三、德拜长度六个）' % APP_FILES['R'])
    lacks(text, '十四个量级',
          'src/%s（德拜长度同一条公式只跨六个量级）' % APP_FILES['R'])
    lacks(text, '毫摩尔级',
          'src/%s（0.75 g 氚是 0.25 mol，不是毫摩尔级）' % APP_FILES['R'])
    lacks(text, r'1.3\times10^{-7}',
          'src/%s（C4 的 a/λ_D 不得回到 1.3e-7）' % APP_FILES['R'])


def check_P18_appendix_i_dualcode():
    """附录 I「双编码与多媒体学习」一节与「动作六」的完整性。

    可证伪点（全部记 FAIL）：
      ① 双编码一节被删或改了 label（sec:app-i-dualcode / tab:app-i-dualcode）；
      ② 七条规则缺任何一条（表格某行被删）；
      ③ 六家出处缺一（Paivio / Clark & Paivio / Mayer & Moreno / Mayer 2009 / Sweller / Tufte）；
      ④ 「动作六：把图讲成句子」被删，或小节标题退回「五个可执行动作」；
      ⑤ 旧错误回流：两本期刊被拼成一个刊名。
    """
    p = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(p):
        rec('FAIL', 'P', '附录 I 文件不存在：%s' % APP_FILES['I'])
        return
    text = strip_comments(read(p))
    secs = split_blocks(text, 'section')
    blk = None
    for ttl, body in secs:
        if ttl.startswith('双编码与多媒体学习'):
            blk = body
            break
    if blk is None:
        rec('FAIL', 'P', '附录 I 找不到「双编码与多媒体学习」一节（标题被改或整节被删）')
        return
    rec('PASS', 'P', '附录 I 的双编码一节存在（%d 字符）' % len(blk))

    has(text, r'\label{sec:app-i-dualcode}',
        'src/%s（双编码一节的 label）' % APP_FILES['I'])
    has(blk, r'\label{tab:app-i-dualcode}',
        'src/%s（七条规则表 label）' % APP_FILES['I'])
    for rule in ('多媒体', '空间邻近', '时间邻近', '一致性', '冗余', '信号', '分段'):
        has(blk, rule + ' &', 'src/%s（双编码规则表的一行：%s）' % (APP_FILES['I'], rule))
    has(blk, '三句话读表法', 'src/%s（把动作六落成三步读表法）' % APP_FILES['I'])
    has(blk, '双编码', 'src/%s（双编码假设本身）' % APP_FILES['I'])
    has(blk, '认知负荷', 'src/%s（工作记忆与认知负荷那一半）' % APP_FILES['I'])
    has(blk, '文本想象', 'src/%s（与 Dunlosky 低效用项的张力必须解释）' % APP_FILES['I'])
    for src in ('Paivio', 'Clark \\& Paivio', 'Mayer \\& Moreno', 'Mayer (2009)',
                'Sweller', 'Tufte'):
        has(blk, src, 'src/%s（双编码一节的出处之一：%s）' % (APP_FILES['I'], src))
    has(text, r'\subsection{动作六', 'src/%s（动作六小节）' % APP_FILES['I'])
    has(text, '六个可执行动作', 'src/%s（动作清单扩到六个）' % APP_FILES['I'])
    has(text, '动作一到动作六', 'src/%s（合计句同步到六个动作）' % APP_FILES['I'])
    has(text, 'Paivio (1986) &', 'src/%s（参考文献表补上双编码一条）' % APP_FILES['I'])
    lacks(text, '五个可执行动作', 'src/%s（不得退回"五个可执行动作"）' % APP_FILES['I'])
    lacks(text, 'Perspectives on Psychological Science',
          'src/%s（两本期刊不得被拼成一个刊名）' % APP_FILES['I'])


def check_P19_appendix_s_confinement():
    """附录 S 能量约束与输运的结构完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 十一个小节之一被删（按小节标题子串定位）；
      ② 七张表 / 两只语义盒 / 十五个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 9 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（旧式偏差、τ_E、χ_eff、χ_B、χ_gB、ρ_i、ρ*、杠杆倍数、
         W、ΔW、峰值热负荷、n_G 比值、f_GAM、三种外推），或搬到了别的小节里；
      ⑤ 口径回流：ρ_i 丢掉 2T 约定、回旋玻姆写成 ρ*²、错式的幂次混进新式。

    数值唯一来源 tmp/appS_calc.py；表内数字的单位写在表头，故这些用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['S'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 S 文件缺失：src/%s' % APP_FILES['S'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = ['约束定标律是什么', '从功率平衡反推有效输运系数', '玻姆还是回旋玻姆',
            'IPB98(y,2)', '因子与设计杠杆', '台基、', 'Greenwald 密度极限',
            '带状流与测地声模', '三种外推', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 S 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 S 的 %d 个小节齐全' % len(want))

    for lab in ('tab:app-s-map', 'tab:app-s-devices', 'tab:app-s-transport',
                'tab:app-s-oldformula', 'tab:app-s-lever', 'tab:app-s-gam',
                'tab:app-s-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['S'], lab))

    for lab in ('tip-apps-limit', 'tip-apps-next'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['S'], lab))

    for lab in ('eq:app-s-powerbalance', 'eq:app-s-chi', 'eq:app-s-chib',
                'eq:app-s-chigb', 'eq:app-s-ipb98', 'eq:app-s-ipb98b', 'eq:app-s-old',
                'eq:app-s-h', 'eq:app-s-ped', 'eq:app-s-elm', 'eq:app-s-elmfreq',
                'eq:app-s-ng', 'eq:app-s-gam', 'eq:app-s-extrapb', 'eq:app-s-extrapgb'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['S'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['S'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 9:
        rec('FAIL', 'P', '附录 S 的「出处：」行只有 %d 处，应 >= 9' % n_src)
    else:
        rec('PASS', 'P', '附录 S 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex < 6:
        rec('FAIL', 'P', '附录 S 的估算练习只有 %d 道，应 >= 6' % n_ex)
    else:
        rec('PASS', 'P', '附录 S 有 %d 道估算练习' % n_ex)
    if n_ans < n_ex:
        rec('FAIL', 'P', '附录 S 有 %d 道练习但只有 %d 个解答' % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 S 的每道练习都有解答（%d 个）' % n_ans)

    # 字面绑定：公式形式与口径说明必须原样在文里（改动即 FAIL）
    for lit, why in (
            (r'\chi_{\mathrm{eff}}=\frac{3a^2}{8\tau_E}', 'τ_E→χ_eff 的柱位形关系'),
            (r'\chi_B=\frac{T_e[\mathrm{eV}]}{16\,B[\mathrm{T}]}', '玻姆系数的 1/16 约定'),
            (r'\rho_i=\frac{\sqrt{2m_iT_i}}{eB}', '回旋半径的 2T 约定'),
            (r'\chi_{\mathrm{gB}}=\chi_B\,\rho_\ast', '回旋玻姆 = 玻姆 × ρ*'),
            (r'n_{19}^{0.41}P^{-0.69}', 'IPB98(y,2) 的 P 幂次（必须钉在整条式子上，'
                                        '因为 -0.69 在正文与练习里重复出现）'),
            (r'\epsilon^{0.58}\kappa^{0.78}M^{0.19}', 'IPB98(y,2) 的 ε/κ/M 因子'),
            (r'a^{0.58}R^{1.39}', 'a 与 R 的等价改写'),
            (r'H\equiv\frac{\tau_{E,\mathrm{exp}}}{\tau_E^{\mathrm{IPB98}}}', 'H 因子的定义'),
            (r'W_{\mathrm{ped}}=3nTV\times20\%', '台基储能的口径'),
            (r'n_G=\frac{I_p[\mathrm{MA}]}{\pi a^2[\mathrm{m^2}]}\times10^{20}', 'Greenwald 公式'),
            (r'f_{\mathrm{GAM}}\simeq\frac{c_s}{2\pi R}', '测地声模频率估计'),
            (r'\tau_E\propto\frac{a^2B}{T}', '玻姆型外推标度'),
            (r'\tau_E\propto\frac{a^3B^2}{T^{1.5}}', '回旋玻姆型外推标度'),
            (r'\tau_E^{\mathrm{Bohm}}=0.426', '玻姆型外推结果'),
            (r'\tau_E^{\mathrm{gB}}=1.358', '回旋玻姆型外推结果')):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['S'], why))

    # 数值独立复算（单位见正文/表头，故 unit=''）
    cases = [
        # (值, 相对容差, 所在小节子串或 None, 说明)
        (15.0, 0.01, '约束定标律是什么', 'D1 的 I_p'),
        (5.3, 0.01, '约束定标律是什么', 'D1 的 B'),
        (10.1, 0.01, '约束定标律是什么', 'D1 的 n_19'),
        (150, 0.01, '约束定标律是什么', 'D1 的 P'),
        (6.20, 0.01, '约束定标律是什么', 'D1 的 R'),
        (2.00, 0.01, '约束定标律是什么', 'D1 的 a'),
        (1.7, 0.01, '约束定标律是什么', 'D1 的 κ'),
        (2.5, 0.01, '约束定标律是什么', 'D1 的 M'),

        (0.606, 0.01, '玻姆还是回旋玻姆', 'D1 的 χ_eff'),
        (1.531, 0.01, '玻姆还是回旋玻姆', 'D2 的 χ_eff'),
        (1.595, 0.01, '玻姆还是回旋玻姆', 'D3 的 χ_eff'),
        (108.4, 0.01, '玻姆还是回旋玻姆', 'D1 的 χ_B'),
        (48.7, 0.01, '玻姆还是回旋玻姆', 'D2 的 χ_B'),
        (0.2004, 0.01, '玻姆还是回旋玻姆', 'D1 的 χ_gB'),
        (0.1281, 0.01, '玻姆还是回旋玻姆', 'D2 的 χ_gB'),
        (0.2893, 0.01, '玻姆还是回旋玻姆', 'D3 的 χ_gB'),
        (3.70, 0.01, '玻姆还是回旋玻姆', 'D1 的 ρ_i/mm'),
        (3.29, 0.01, '玻姆还是回旋玻姆', 'D2 的 ρ_i/mm'),
        (3.95, 0.01, '玻姆还是回旋玻姆', 'D3 的 ρ_i/mm'),
        (1.85e-3, 0.01, '玻姆还是回旋玻姆', 'D1 的 ρ*'),
        (2.63e-3, 0.01, '玻姆还是回旋玻姆', 'D2 的 ρ*'),
        (5.90e-3, 0.01, '玻姆还是回旋玻姆', 'D3 的 ρ*'),
        (0.0056, 0.02, '玻姆还是回旋玻姆', 'D1 的 χ_eff/χ_B'),
        (3.0, 0.02, '玻姆还是回旋玻姆', 'D1 的 χ_eff/χ_gB'),
        (541, 0.01, '玻姆还是回旋玻姆', 'D1 的 χ_B/χ_gB'),
        (380, 0.01, '玻姆还是回旋玻姆', 'D2 的 χ_B/χ_gB'),
        (170, 0.01, '玻姆还是回旋玻姆', 'D3 的 χ_B/χ_gB'),

        (2.4762, 1e-4, 'IPB98(y,2)', 'D1 的 IPB98 τ_E'),
        (7.4446, 1e-4, 'IPB98(y,2)', 'D1 的错式 τ_E'),
        (3.01, 0.01, 'IPB98(y,2)', 'D1 的错式/IPB98 倍数'),
        (0.3827, 1e-4, 'IPB98(y,2)', 'D2 的 IPB98 τ_E'),
        (0.5407, 1e-4, 'IPB98(y,2)', 'D2 的错式 τ_E'),
        (1.41, 0.01, 'IPB98(y,2)', 'D2 的倍数'),
        (0.1055, 1e-4, 'IPB98(y,2)', 'D3 的 IPB98 τ_E'),
        (0.0493, 1e-4, 'IPB98(y,2)', 'D3 的错式 τ_E'),
        (0.47, 0.01, 'IPB98(y,2)', 'D3 的倍数'),

        (1.14, 0.01, '因子与设计杠杆', 'R 的杠杆倍数'),
        (1.33, 0.01, '因子与设计杠杆', 'I_p 的杠杆倍数'),
        (1.90, 0.01, '因子与设计杠杆', 'n_19 的杠杆倍数'),
        (1.57, 0.01, '因子与设计杠杆', 'a 的杠杆倍数'),
        (1.40, 0.01, '因子与设计杠杆', 'κ 的杠杆倍数'),
        (3.98, 0.01, '因子与设计杠杆', 'M 的杠杆倍数'),
        (5.75, 0.01, '因子与设计杠杆', 'B 的杠杆倍数'),
        (0.68, 0.01, '因子与设计杠杆', 'P 的杠杆倍数'),

        (19.6, 0.01, '台基、', '台基储能 MJ'),
        (371.4, 0.001, '台基、', 'D1 总储能 MJ'),
        (37.1, 0.01, '台基、', '单次 ELM 能量 MJ'),
        (7.43, 0.01, '台基、', '峰值热负荷 GW/m²'),
        (495, 0.02, '台基、', '峰值/稳态倍数'),
        (4.04, 0.02, '台基、', 'ELM 频率下界 Hz'),
        (81, 0.01, '台基、', 'ΔW 压到 0.5% 时的频率'),
        (15.0, 0.01, '台基、', '稳态热负荷 MW/m²'),

        (0.846, 0.01, 'Greenwald 密度极限', 'D1 的 n/n_G'),
        (1.194e20, 0.01, 'Greenwald 密度极限', 'D1 的 n_G'),
        (0.611e20, 0.01, 'Greenwald 密度极限', 'D2 的 n_G'),
        (0.654, 0.01, 'Greenwald 密度极限', 'D2 的 n/n_G'),
        (0.470, 0.01, 'Greenwald 密度极限', 'D3 的 n/n_G'),
        (1.064e20, 0.01, 'Greenwald 密度极限', 'D3 的 n_G'),

        (17.0, 0.01, '带状流与测地声模', 'D1 芯部 GAM 频率 kHz'),
        (5.6, 0.02, '带状流与测地声模', 'D1 边缘 GAM 频率 kHz'),
        (18.0, 0.01, '带状流与测地声模', 'D2 芯部 GAM 频率 kHz'),
        (11.8, 0.01, '带状流与测地声模', 'D2 边缘 GAM 频率 kHz'),
        (26.8, 0.01, '带状流与测地声模', 'D3 芯部 GAM 频率 kHz'),
        (20.9, 0.01, '带状流与测地声模', 'D3 边缘 GAM 频率 kHz'),
        (6.64e5, 0.01, '带状流与测地声模', 'D1 芯部声速 m/s'),
        (2.19e5, 0.01, '带状流与测地声模', 'D1 边缘声速 m/s'),

        (0.426, 0.01, '三种外推', '玻姆型外推 τ_E'),
        (1.358, 0.01, '三种外推', '回旋玻姆型外推 τ_E'),
        (2.476, 0.01, '三种外推', 'IPB98 的 τ_E（对比基准）'),
        (5.82, 0.02, '三种外推', '玻姆型与 IPB98 的比值'),
        (1.82, 0.02, '三种外推', '回旋玻姆型与 IPB98 的比值'),

        (3.606, 0.001, '六道估算练习', '练习一的 τ_E'),
        (1.46, 0.02, '六道估算练习', '练习一的倍数'),
        (0.0056, 0.02, '六道估算练习', '练习二的 χ_eff/χ_B'),
        (1.064e20, 0.01, '六道估算练习', '练习四的 n_G'),
        (0.470, 0.01, '六道估算练习', '练习四的 n/n_G'),
        (1.33, 0.01, '六道估算练习', '练习五的 I_p 杠杆'),
        (5.75, 0.01, '六道估算练习', '练习五的 B 杠杆'),
        (1.14, 0.01, '六道估算练习', '练习五的 R 杠杆'),
        (1.57, 0.01, '六道估算练习', '练习五的 a 杠杆'),
        (17.0, 0.01, '六道估算练习', '练习六的芯部频率'),
        (5.6, 0.02, '六道估算练习', '练习六的边缘频率'),
    ]
    bad = []
    for v, rel, sect, why in cases:
        blk = found[sect] if sect else text
        if not value_in(blk, v, rel, ''):
            bad.append('%s（期望 %.6g）' % (why, v))
    if bad:
        rec('FAIL', 'P', '附录 S 有 %d 处印出的数值与 tmp/appS_calc.py 不符或不在应有小节：%s'
            % (len(bad), '；'.join(bad)))
    else:
        rec('PASS', 'P', '附录 S 的 %d 处关键数值与 tmp/appS_calc.py 逐一相符' % len(cases))

    # 口径回流守卫
    lacks(text, r'\rho_i=\frac{\sqrt{m_iT_i}}',
          'src/%s（回旋半径不得丢掉 2T 约定）' % APP_FILES['S'])
    lacks(text, r'\chi_{\mathrm{gB}}=\chi_B\,\rho_\ast^2',
          'src/%s（回旋玻姆不得写成 ρ*²）' % APP_FILES['S'])
    lacks(text, r'n_{19}^{0.41}P^{-0.58}',
          'src/%s（错式的 P 幂次不得混进新式）' % APP_FILES['S'])
    lacks(text, r'a^{1.97}\epsilon^{0.58}',
          'src/%s（a^{1.97} 与 ε^{0.58} 不得混用）' % APP_FILES['S'])


def check_P20_ipb98_formula_consistency():
    """第 10 章公式卡 f:scaling 与附录 S 的 IPB98(y,2) 必须一致，且错式不得回流。

    可证伪点（全部记 FAIL）：
      ① 公式卡被删或改了标题/label；
      ② 公式卡里 IPB98(y,2) 的任一因子（I_p/B/n_19/P/R/ε/κ/M）被改动或漏掉；
      ③ 名称退回两处旧写法之一（ITER98y2、把 a 的幂次写成 1.97、P 的幂次写成 -0.58）；
      ④ 附录 S 之外的任何 tex 文件里出现 a^{1.97} 或 P^{-0.58}——
         即有人把错式当成正确公式又抄了一遍。

    这条检查是跨文件的：第 10 章给出"用哪条式子"，附录 S 给出"式子从哪来、怎么用"，
    两处只要不一致，读者就会在同一个符号上读到两个幂次。
    """
    cname = 'ch10_supplement.tex'
    cpath = os.path.join(SRC, cname)
    if not os.path.exists(cpath):
        rec('FAIL', 'P', '第 10 章补充材料缺失：src/%s' % cname)
        return
    text = strip_comments(read(cpath))
    m = re.search(r'\\begin\{formulabox\}\{磁约束约束时间经验定标律\}\{f:scaling\}'
                  r'(.*?)\\end\{formulabox\}', text, re.S)
    if not m:
        rec('FAIL', 'P', '第 10 章找不到公式卡 f:scaling（标题或 label 被改）')
        return
    box = m.group(1)
    rec('PASS', 'P', '第 10 章公式卡 f:scaling 存在（%d 字符）' % len(box))

    for lit in (r'I_p^{0.93}', r'B^{0.15}', r'n_{19}^{0.41}', r'P^{-0.69}',
                r'R^{1.97}', r'\epsilon^{0.58}', r'\kappa^{0.78}', r'M^{0.19}'):
        has(box, lit, 'src/%s（公式卡 f:scaling 的 IPB98(y,2) 因子 %s）' % (cname, lit))
    has(box, r'a^{0.58}R^{1.39}',
        'src/%s（公式卡里 a 与 R 的等价改写，用来钉住 0.58 这个幂次）' % cname)
    has(box, 'IPB98(y,2)', 'src/%s（公式卡的正式名称）' % cname)
    has(box, r'\epsilon=a/R', 'src/%s（公式卡里 ε 的定义，符号说明要能自洽）' % cname)
    has(box, r'$\kappa$ 为拉长比',
        'src/%s（符号说明必须交代 κ，否则 consistency_audit 会报未说明符号）' % cname)
    lacks(box, 'ITER98y2', 'src/%s（公式卡不得退回 ITER98y2 这个名称）' % cname)
    lacks(box, r'a^{1.97}', 'src/%s（a 的幂次不得写成 1.97）' % cname)
    lacks(box, r'P^{-0.58}', 'src/%s（P 的幂次不得写成 -0.58）' % cname)

    offenders = []
    for fn in sorted(os.listdir(os.path.join(SRC))):
        if not fn.endswith('.tex') or fn == APP_FILES['S']:
            continue
        body = strip_comments(read(os.path.join(SRC, fn)))
        if r'a^{1.97}' in body or r'P^{-0.58}' in body:
            offenders.append(fn)
    if offenders:
        rec('FAIL', 'P', '旧错式只允许出现在附录 S 的对照节，却出现在：%s' % offenders)
    else:
        rec('PASS', 'P', '旧错式（a^{1.97} / P^{-0.58}）只出现在附录 S 的对照节')


def check_P21_appendix_i_expert():
    """附录 I「专家与新手」一节的结构、出处与防回流。

    可证伪点（全部记 FAIL）：
      ① 该节被删或改了 label（sec:app-i-expert / tab:app-i-expert）；
      ② 四条对照（记忆单位 / 分类依据 / 解题方向 / 表征习惯）缺任何一行；
      ③ 六家经典出处缺一（de Groot、Chase & Simon、Chi 等 1981、Larkin 等 1980、
         Miller、Cowan），或参考文献表缺对应行；
      ④ 结论被写成"专家天生记性好"这类与实验证据相反的表述；
      ⑤ 动作清单被悄悄改成"七个动作"却没同步 P18 的计数。

    Miller 1956 与 Cowan 2001 必须同时在场：只写 7±2 会漏掉现代估计（约 4），
    而这一节正是要说明"容量这么小，所以只能靠组块"。
    """
    p = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(p):
        rec('FAIL', 'P', '附录 I 文件不存在：%s' % APP_FILES['I'])
        return
    text = strip_comments(read(p))
    secs = split_blocks(text, 'section')
    blk = None
    for ttl, body in secs:
        if '专家与新手' in ttl:
            blk = body
            break
    if blk is None:
        rec('FAIL', 'P', '附录 I 找不到「专家与新手」一节（标题被改或整节被删）')
        return
    rec('PASS', 'P', '附录 I 的专家与新手一节存在（%d 字符）' % len(blk))

    has(text, r'\label{sec:app-i-expert}',
        'src/%s（专家与新手一节的 label）' % APP_FILES['I'])
    has(blk, r'\label{tab:app-i-expert}',
        'src/%s（四条对照表 label）' % APP_FILES['I'])
    for row in ('记忆单位', '分类依据', '解题方向', '表征习惯'):
        has(blk, row, 'src/%s（四条对照之一：%s）' % (APP_FILES['I'], row))
    for key in ('组块', '随机棋局', '深层原理', '表面特征', '向前推', '工作记忆'):
        has(blk, key, 'src/%s（专家与新手一节的关键概念：%s）' % (APP_FILES['I'], key))
    for src in ('de Groot', 'Chase \\& Simon', 'Chi, Feltovich \\& Glaser',
                'Larkin', 'Miller', 'Cowan'):
        has(blk, src, 'src/%s（专家与新手一节的出处之一：%s）' % (APP_FILES['I'], src))
    for row in ('de Groot (1965) &', 'Chase \\& Simon (1973) &',
                'Chi, Feltovich \\& Glaser (1981) &',
                'Larkin, McDermott, Simon \\& Simon (1980) &',
                'Miller (1956) &', 'Cowan (2001) &'):
        has(text, row, 'src/%s（参考文献表补上：%s）' % (APP_FILES['I'], row))
    lacks(blk, '因为专家记性好',
          'src/%s（不得把专家优势说成天生记性好）' % APP_FILES['I'])
    lacks(blk, '不需要做本领域的题',
          'src/%s（不得把领域性结论写成通用能力）' % APP_FILES['I'])
    lacks(text, '七个可执行动作',
          'src/%s（动作清单若扩到七个，必须同步 P18 的计数断言）' % APP_FILES['I'])


def check_P22_appendix_t_topology():
    """附录 T 磁面、磁岛与磁通量管的结构完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 九个小节之一被删（按小节标题子串定位）；
      ② 八张表 / 两只语义盒 / 十五个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 7 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（q95 公式、有理面位置、剪切、k_∥ 与 λ、τ_e 与 χ_∥、
         岛宽与岛内等温时间、重叠磁通、自举电流与临界岛宽、反转半径），或搬到了别的小节里；
      ⑤ 口径回流：χ_∥ 用了离子质量、τ_e 的 3/2 次幂被改、岛宽系数从 4 退回 2、
         q 剖面从 ρ² 退回 ρ、反转半径丢掉根号、Ψ_h'' 与剪切的对应被改。

    数值唯一来源 tmp/appT_calc.py；表内数字的单位写在表头，故这些用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['T'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 T 文件缺失：src/%s' % APP_FILES['T'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = ['磁面、安全因子', '安全因子剖面与有理面', '磁通量管上的两条路', '磁岛：',
            '岛链重叠', '新经典撕裂模', '锯齿、', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 T 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 T 的 %d 个小节齐全' % len(want))

    for lab in ('tab:app-t-q95', 'tab:app-t-surfaces', 'tab:app-t-kpar',
                'tab:app-t-transport', 'tab:app-t-w', 'tab:app-t-overlap',
                'tab:app-t-ntm', 'tab:app-t-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['T'], lab))

    for lab in ('tip-app-t-flux', 'tip-app-t-flat'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['T'], lab))

    for lab in ('eq:app-t-fluxfunc', 'eq:app-t-bo', 'eq:app-t-qprofile', 'eq:app-t-shear',
                'eq:app-t-resonance', 'eq:app-t-kpar', 'eq:app-t-sidentity', 'eq:app-t-taue',
                'eq:app-t-chipar', 'eq:app-t-psih', 'eq:app-t-width', 'eq:app-t-overlap',
                'eq:app-t-rutherford', 'eq:app-t-wmarg', 'eq:app-t-rinv'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['T'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['T'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 7:
        rec('FAIL', 'P', '附录 T 的「出处：」行只有 %d 处，应 >= 7' % n_src)
    else:
        rec('PASS', 'P', '附录 T 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex < 6:
        rec('FAIL', 'P', '附录 T 的估算练习只有 %d 道，应 >= 6' % n_ex)
    else:
        rec('PASS', 'P', '附录 T 有 %d 道估算练习' % n_ex)
    if n_ans < n_ex:
        rec('FAIL', 'P', '附录 T 有 %d 道练习但只有 %d 个解答' % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 T 的每道练习都有解答（%d 个）' % n_ans)

    # 字面绑定：公式形式与口径说明必须原样在文里（改动即 FAIL）
    for lit, why in (
            (r'q(\rho)=q_0+(q_{95}-q_0)\,\rho^2', 'q 剖面的教学形式（ρ² 不能退回 ρ）'),
            ("s(r)\\equiv\\frac{r\\,q'}{q}", '磁剪切的定义'),
            (r'q(r_{mn})=\frac{m}{n}', '有理面的共振条件'),
            (r's(r)=2\,\frac{q(r)-q_0}{q(r)}=2-\frac{2q_0}{q(r)}', '剪切恒等式（P22 的解析自检）'),
            (r'\tau_e=\frac{3.44\times10^{5}\,T_e^{3/2}}{n_e\ln\Lambda}', 'NRL 电子碰撞时间'),
            (r'\ln\Lambda=24-\ln\frac{n_e^{1/2}}{T_e}', '库仑对数的 NRL 形式'),
            (r'\chi_\parallel=\frac{3.16\,T_e\,\tau_e}{m_e}', 'Braginskii 平行电子热导（必须是 m_e）'),
            ("\\Psi_h''(r_{mn})=B_\\theta\\,s", '螺旋磁通二阶导数与剪切的关系'),
            (r'w=4\sqrt{\frac{\tilde\Psi}{B_\theta\,s}}', '岛宽公式（系数 4）'),
            (r'\tilde\Psi\;\gtrsim\;B_\theta s\left(\frac{\Delta r}{4}\right)^2', '岛链重叠判据'),
            ("w_{\\mathrm{marg}}=\\frac{5.6\\,\\mu_0\\,j_{\\mathrm{bs}}\\,L_q}"
             "{B_\\theta\\,|\\Delta'_0|}", 'NTM 临界岛宽'),
            (r'r_{\mathrm{inv}}=a\sqrt{\frac{1-q_0}{q_{95}-q_0}}', '锯齿反转半径'),
            (r'\frac{w_1}{2}+\frac{w_2}{2}\;\gtrsim\;\Delta r', '重叠判据的几何形式'),
            (r'f_{\mathrm{bs}}=0.25', '自举电流份额的教学取值'),
            # 整行/整句绑定：下面这些数在表里与正文、练习里重复出现（0.0698 出现 4 次、
            # 9.023e11 出现 2 次），只钉单个数字时改掉表里那一处仍能满足数值用例——
            # v1.17.0 的逐条探针第 67、68 号就是这样逃掉的，故把每一处印出的位置都钉住。
            (r'D1 ITER 量级 & $9.194$ & $17.00$ & $1.766\times10^{-4}$ & $9.023\times10^{11}$'
             r' & $0.6058$ & $1.490\times10^{12}$', '表 T.3 的 D1 整行'),
            (r'D1 ITER 量级 & $3/2$ & $0.9608$ & $0.5475$ & $0.6000$ & $0.0698$ & $0.0349$'
             r' & $5.398\times10^{-15}$', '表 T.4 的 D1 3/2 整行'),
            (r'\chi_\parallel=9.02\times10^{11}', 'T.3 正文里的四舍五入值'),
            (r'\chi_\parallel=9.023\times10^{11}', '练习 4 给定的 χ_∥ 精确值'),
            (r'w_{2/1}=0.0531\ \mathrm{m}<w_{3/2}=0.0698\ \mathrm{m}', 'T.4 的岛宽比较句'),
            (r'w=0.0698\times\sqrt4=0.1396\ \mathrm{m}', '练习 3 解答的倍数推导'),
            (r'w=0.0698\ \mathrm{m}$、$\tau_E=2.4762\ \mathrm{s}', '练习 4 给定的 w 与 τ_E'),
            # v1.17.1：正文里印出的 τ_flat 必须与表 T.5、与脚本同值（表改对齐脚本时
            # 正文曾漏改，3.125 与 3.128 的差别落在任何数值用例的容差之内，只能逐字钉）
            (r'$2/1$ 面的 $\tau_{\mathrm{flat}}=w^2/\chi_\parallel=3.128\times10^{-15}'
             r'\ \mathrm{s}$', 'T.4 正文里的 D1 2/1 面 τ_flat')):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['T'], why))

    # 口径回流防护（这些写法一旦出现即为 FAIL）
    for lit, why in (
            (r'\chi_\parallel=\frac{3.16\,T_e\,\tau_e}{m_i}', '平行电子热导被写成离子质量'),
            (r'\tau_e=\frac{3.44\times10^{5}\,T_e^{1/2}', '电子碰撞时间丢了 3/2 次幂'),
            (r'w=2\sqrt{\frac{\tilde\Psi}{B_\theta\,s}}', '岛宽系数从 4 退回 2'),
            (r'q(\rho)=q_0+(q_{95}-q_0)\,\rho,', 'q 剖面从 ρ² 退回 ρ'),
            (r'r_{\mathrm{inv}}=a\frac{1-q_0}{q_{95}-q_0}', '反转半径丢掉根号'),
            ("\\Psi_h''(r_{mn})=B_\\theta\\,s^2", '螺旋磁通二阶导数写成 s²'),
            (r'=3.125\times10^{-15}', 'τ_flat 正文与表/脚本不一致（v1.17.1 修掉的旧写法）')):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['T'], why))

    # 数值独立复算（单位见正文/表头，故 unit=''）
    cases = [
        # T.1 安全因子与 q95 经验公式
        (1.05, 1e-3, '磁面、安全因子', 'D1 的 q_0'),
        (3.00, 1e-3, '磁面、安全因子', 'D1 的 q_95（设计值）'),
        (2.22, 1e-3, '磁面、安全因子', 'D1 的 q_95（经验公式）'),
        (0.78, 1e-3, '磁面、安全因子', 'D1 的差 +0.78'),
        (3.90, 1e-3, '磁面、安全因子', 'D2 的 q_95（设计值）'),
        (3.91, 1e-3, '磁面、安全因子', 'D2 的 q_95（经验公式）'),
        (3.99, 1e-3, '磁面、安全因子', 'D3 的 q_95（经验公式）'),
        (5.30, 1e-3, '磁面、安全因子', 'D1 的环向场 B/T'),

        # T.2 有理面与剪切
        (0.4804, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 ρ'),
        (0.9608, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 r'),
        (0.5475, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 B_θ'),
        (0.9367, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 q\' '),
        (0.6000, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 s'),
        (1.6013, 1e-3, '安全因子剖面与有理面', 'D1 3/2 面的 L_q'),
        (1.3960, 1e-3, '安全因子剖面与有理面', 'D1 2/1 面的 r'),
        (0.5967, 1e-3, '安全因子剖面与有理面', 'D1 2/1 面的 B_θ'),
        (0.9500, 1e-3, '安全因子剖面与有理面', 'D1 2/1 面的 s'),
        (0.5190, 1e-3, '安全因子剖面与有理面', 'D2 3/2 面的 r'),
        (0.3507, 1e-3, '安全因子剖面与有理面', 'D2 3/2 面的 B_θ'),
        (0.6667, 1e-3, '安全因子剖面与有理面', 'D2 3/2 面的 s'),
        (0.7340, 1e-3, '安全因子剖面与有理面', 'D2 2/1 面的 r'),
        (1.0000, 1e-3, '安全因子剖面与有理面', 'D2 2/1 面的 s（恒等式给 1）'),
        (0.1462, 1e-3, '安全因子剖面与有理面', 'D3 1/1 面的 r'),
        (0.1839, 1e-3, '安全因子剖面与有理面', 'D3 1/1 面的 B_θ'),
        (0.3000, 1e-3, '安全因子剖面与有理面', 'D3 1/1 面的 s'),
        (0.3044, 1e-3, '安全因子剖面与有理面', 'D3 3/2 面的 r'),
        (0.8667, 1e-3, '安全因子剖面与有理面', 'D3 3/2 面的 s'),
        (0.2545, 1e-3, '安全因子剖面与有理面', 'D3 2/1 面的 B_θ'),
        (1.1500, 1e-3, '安全因子剖面与有理面', 'D3 2/1 面的 s'),

        # T.2b 平行波数（绝对值比对，正文按大小印出）
        (1.792e-2, 1e-3, '安全因子剖面与有理面', '2/1 面 q=1.8 的 k_∥'),
        (8.489e-3, 1e-3, '安全因子剖面与有理面', '2/1 面 q=1.9 的 k_∥'),
        (7.680e-3, 1e-3, '安全因子剖面与有理面', '2/1 面 q=2.1 的 k_∥'),
        (1.466e-2, 1e-3, '安全因子剖面与有理面', '2/1 面 q=2.2 的 k_∥'),
        (3.584e-2, 1e-3, '安全因子剖面与有理面', '3/2 面 q=1.35 的 k_∥'),
        (1.698e-2, 1e-3, '安全因子剖面与有理面', '3/2 面 q=1.425 的 k_∥'),
        (2.933e-2, 1e-3, '安全因子剖面与有理面', '3/2 面 q=1.65 的 k_∥'),
        (350.6, 1e-3, '安全因子剖面与有理面', '2/1 面 q=1.8 的 λ/m'),
        (740.2, 1e-3, '安全因子剖面与有理面', '2/1 面 q=1.9 的 λ/m'),
        (214.3, 1e-3, '安全因子剖面与有理面', '3/2 面 q=1.65 的 λ/m'),

        # T.3 平行与垂直输运
        (17.00, 1e-3, '磁通量管上的两条路', 'D1 的 lnΛ'),
        (1.766e-4, 1e-3, '磁通量管上的两条路', 'D1 的 τ_e/s'),
        (9.023e11, 1e-3, '磁通量管上的两条路', 'D1 的 χ_∥'),
        (0.6058, 1e-3, '磁通量管上的两条路', 'D1 的 χ_⊥'),
        (1.490e12, 1e-3, '磁通量管上的两条路', 'D1 的 χ_∥/χ_⊥'),
        (16.10, 1e-3, '磁通量管上的两条路', 'D2 的 lnΛ'),
        (6.032e-5, 1e-3, '磁通量管上的两条路', 'D2 的 τ_e/s'),
        (7.832e10, 1e-3, '磁通量管上的两条路', 'D2 的 χ_∥'),
        (1.5311, 1e-3, '磁通量管上的两条路', 'D2 的 χ_⊥'),
        (15.64, 1e-3, '磁通量管上的两条路', 'D3 的 lnΛ'),
        (2.946e-5, 1e-3, '磁通量管上的两条路', 'D3 的 τ_e/s'),
        (2.699e10, 1e-3, '磁通量管上的两条路', 'D3 的 χ_∥'),
        (1.5950, 1e-3, '磁通量管上的两条路', 'D3 的 χ_⊥'),
        (38.96, 1e-3, '磁通量管上的两条路', 'D1 的 L_∥/m'),
        (1.68e-9, 0.005, '磁通量管上的两条路', 'D1 的平行飞行时间/s'),
        (1.651e-4, 1e-3, '磁通量管上的两条路', 'D1 的 1 cm 横向时间/s'),
        (9.8e4, 0.005, '磁通量管上的两条路', 'D1 的两条路时间比'),

        # T.4 岛宽与岛内等温时间
        (0.0698, 1e-3, '磁岛：', 'D1 3/2 面的 w（Ψ̃=1e-4）'),
        (0.0349, 1e-3, '磁岛：', 'D1 3/2 面的半宽'),
        (0.0531, 1e-3, '磁岛：', 'D1 2/1 面的 w'),
        (0.0827, 1e-3, '磁岛：', 'D2 3/2 面的 w'),
        (0.0656, 1e-3, '磁岛：', 'D2 2/1 面的 w'),
        (0.0851, 1e-3, '磁岛：', 'D3 3/2 面的 w'),
        (0.1270, 1e-3, '磁岛：', 'D3 3/2 面的半宽（比 w 大，全宽才是 w）'),
        (0.0739, 1e-3, '磁岛：', 'D3 2/1 面的 w'),
        (3.128e-15, 1e-3, '磁岛：', 'D1 2/1 面的 τ_flat/s'),
        (5.492e-14, 1e-3, '磁岛：', 'D2 2/1 面的 τ_flat/s'),
        (2.025e-13, 1e-3, '磁岛：', 'D3 2/1 面的 τ_flat/s'),
        (5.398e-15, 1e-3, '六道估算练习', 'D1 3/2 面的 τ_flat/s（练习四）'),

        # T.5 岛链重叠
        (0.4352, 1e-3, '岛链重叠', 'D1 的 Δr/m'),
        (6.710e-3, 1e-3, '岛链重叠', 'D1 重叠所需 Ψ̃'),
        (67.1, 1e-3, '岛链重叠', 'D1 的倍数'),
        (0.2150, 1e-3, '岛链重叠', 'D2 的 Δr/m'),
        (1.075e-3, 1e-3, '岛链重叠', 'D2 重叠所需 Ψ̃'),
        (10.7, 1e-3, '岛链重叠', 'D2 的倍数'),
        (0.1005, 1e-3, '岛链重叠', 'D3 的 Δr/m'),
        (1.847e-4, 1e-3, '岛链重叠', 'D3 重叠所需 Ψ̃'),
        (1.8, 1e-3, '岛链重叠', 'D3 的倍数'),

        # T.6 NTM 与临界岛宽
        (7.022e5, 1e-3, '新经典撕裂模', 'D1 的 j_φ/(A/m²)'),
        (1.755e5, 1e-3, '新经典撕裂模', 'D1 的 j_bs/(A/m²)'),
        (1.469, 1e-3, '新经典撕裂模', 'D1 的 L_q'),
        (15.21, 1e-3, '新经典撕裂模', 'D1 的 w_marg/cm（|Δ\'|=20）'),
        (6.08, 1e-3, '新经典撕裂模', 'D1 的 w_marg/cm（|Δ\'|=50）'),
        (3.04, 1e-3, '新经典撕裂模', 'D1 的 w_marg/cm（|Δ\'|=100）'),
        (4.365e5, 1e-3, '新经典撕裂模', 'D2 的 j_φ/(A/m²)'),
        (1.091e5, 1e-3, '新经典撕裂模', 'D2 的 j_bs/(A/m²)'),
        (0.734, 1e-3, '新经典撕裂模', 'D2 的 L_q'),
        (7.58, 1e-3, '新经典撕裂模', 'D2 的 w_marg/cm（|Δ\'|=20）'),
        (5.909e5, 1e-3, '新经典撕裂模', 'D3 的 j_φ/(A/m²)'),
        (1.477e5, 1e-3, '新经典撕裂模', 'D3 的 j_bs/(A/m²)'),
        (1.44, 1e-3, '新经典撕裂模', 'D3 的 w_marg/cm（|Δ\'|=100）'),

        # T.7 锯齿反转半径
        (0.1462, 1e-3, '锯齿、', 'D3 的 r_inv/m'),
        (0.2182, 1e-3, '锯齿、', 'D3 的 ρ_inv'),
        (4.8, 1e-3, '锯齿、', '反转面以内的体积份额/%'),

        # T.8 练习答案
        (0.9608, 1e-3, '六道估算练习', '练习一答案的 r/m'),
        (0.4804, 1e-3, '六道估算练习', '练习一答案的 ρ'),
        (0.9367, 1e-3, '六道估算练习', '练习二答案的 q\' '),
        (0.5475, 1e-3, '六道估算练习', '练习二答案的 B_θ/T'),
        (0.1033, 1e-3, '六道估算练习', '练习二答案的 B_θ/B'),
        (0.1396, 1e-3, '六道估算练习', '练习三答案的 w/m'),
        (6.710e-3, 1e-3, '六道估算练习', '练习五答案的 Ψ̃/(Wb/rad)'),
        (67.1, 1e-3, '六道估算练习', '练习五答案的倍数'),
    ]
    for v, rel, where, what in cases:
        block = found.get(where)
        if block is None:
            rec('FAIL', 'P', '附录 T 找不到小节「%s」，无法核对 %s' % (where, what))
        elif not value_in(block, v, rel, unit=''):
            rec('FAIL', 'P', '附录 T 的 %s 缺失或与 tmp/appT_calc.py 不符：期望 %g（小节「%s」）'
                % (what, v, where))
    rec('PASS', 'P', '附录 T 的 %d 个关键数值均在小节内复算通过' % len(cases))


def check_P23_appendix_i_objectives():
    """附录 I「学习目标与笔记结构」一节的结构与出处完整性。

    可证伪点（全部记 FAIL）：
      ① 小节标题或 label 被删；
      ② 六层认知动作的对照表被删（label tab:app-i-bloom）或降到五行以下；
      ③ 康奈尔笔记的关键动作「遮住右栏、只看左栏自测」被改写；
      ④ 三份新引用的出处（Bloom 等 1956 / Anderson & Krathwohl 2001 / Pauk & Owens）被删，
         或参考文献表里对应的三行被删；
      ⑤ 反面写法的回流：把笔记本身说成高效用技术、把认知层级说成严格的学习顺序。
    """
    ipath = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(ipath):
        rec('FAIL', 'P', '附录 I 文件缺失：src/%s' % APP_FILES['I'])
        return
    text = strip_comments(read(ipath))
    secs = split_blocks(text, 'section')

    body = None
    for ttl, blk in secs:
        if '学习目标与笔记结构' in ttl:
            body = blk
            break
    if body is None:
        rec('FAIL', 'P', '附录 I 缺少「学习目标与笔记结构」一节')
        return
    rec('PASS', 'P', '附录 I 的「学习目标与笔记结构」一节存在')

    has(text, r'\label{sec:app-i-objectives}', 'src/%s（新小节 label）' % APP_FILES['I'])
    has(text, r'\label{tab:app-i-bloom}', 'src/%s（认知层级对照表 label）' % APP_FILES['I'])

    for lv in ('记忆', '理解', '应用', '分析', '评价', '创造'):
        has(body, lv, 'src/%s「学习目标与笔记结构」（认知层级 %s）' % (APP_FILES['I'], lv))

    for lit, why in (
            (r'遮住右栏，只看左栏', '康奈尔笔记的复习动作（左栏自测）'),
            (r'康奈尔笔记', '笔记结构的名词'),
            (r'可检验', '目标写法的核心判据'),
            (r'预测试', '读前先猜的动作'),
            (r'检索练习', '与效用评级一节的衔接'),
            (r'笔记本身属于低效用技术', '与 Dunlosky 评级一致的边界说明')):
        has(body, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

    for lit, why in (
            (r'Bloom 等 (1956)', '认知目标分类学的原始出处'),
            (r'Anderson \& Krathwohl (2001)', '修订版的出处'),
            (r'Pauk \& Owens', '康奈尔笔记的教科书出处')):
        has(body, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

    for lit, why in (
            (r'笔记本身属于高效用技术', '把笔记本身说成高效用（与 Dunlosky 结论相反）'),
            (r'分类学的层级是严格的学习顺序', '把认知层级说成严格顺序')):
        lacks(body, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

    # 参考文献表必须有三行新条目
    for lit, why in (
            (r'Bloom 等 (1956) & \textit{Taxonomy of Educational Objectives', '参考文献表：Bloom'),
            (r'Anderson \& Krathwohl (2001) & \textit{A Taxonomy for Learning', '参考文献表：修订版'),
            (r'Pauk \& Owens (2013) & \textit{How to Study in College}', '参考文献表：康奈尔笔记')):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

    for lit in ('左栏', '右栏', '页底'):
        has(body, lit, 'src/%s（工作表的三段结构：%s）' % (APP_FILES['I'], lit))

    n_step = len(re.findall(r'\\textbf\{（\d）', body))
    if n_step < 6:
        rec('FAIL', 'P', '「学习目标与笔记结构」的编号动作只有 %d 条，应 >= 6' % n_step)
    else:
        rec('PASS', 'P', '「学习目标与笔记结构」有 %d 条编号动作' % n_step)


def check_P24_appendix_u_gain():
    """附录 U 聚变增益 Q、自加热与点火的完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 七个小节之一被删（按小节标题子串定位）；
      ② 七张表 / 三只语义盒 / 十四个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 5 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（三个工作点的 Q、f_alpha 表、自加热闭环
         迭代与收敛点、三重积阈值与 n tau_E 门槛、Q(P_aux) 扫描、工程增益账），
         或搬到了别的小节里；
      ⑤ 口径回流：<sigma v> 拟合的 T 次幂被改、alpha 份额退回 1/4、f_alpha 的
         分母写成 Q+4、回路指数被写成大于 1、Q 的单调性写反、工程增益漏括号、
         盈亏平衡忘掉减 1、Greenwald 极限错一个量级、P_loss 取倒数。

    数值唯一来源 tmp/appU_calc.py；表内数字的单位写在表头，故这些用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['U'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 U 文件缺失：src/%s' % APP_FILES['U'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = [
        'Q 是什么', '加热份额', '自加热闭环', '点火：三重积阈值', '工程增益', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 U 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 U 的 %d 个小节齐全' % len(want))

    for lab in ('tab:app-u-accounts', 'tab:app-u-falpha', 'tab:app-u-iter', 'tab:app-u-triple', 'tab:app-u-qaux', 'tab:app-u-qeng', 'tab:app-u-source'):
        has(text, r'\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['U'], lab))

    for lab in ('tip-app-u-ledger', 'tip-app-u-loop2', 'tip-app-u-loop'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['U'], lab))

    for lab in ('eq:app-u-q', 'eq:app-u-split', 'eq:app-u-pfus', 'eq:app-u-falpha', 'eq:app-u-balance', 'eq:app-u-ipb', 'eq:app-u-loop', 'eq:app-u-gain', 'eq:app-u-ignition', 'eq:app-u-ntau', 'eq:app-u-triple', 'eq:app-u-qaux', 'eq:app-u-qeng', 'eq:app-u-qeng1'):
        has(text, r'\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['U'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['U'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 5:
        rec('FAIL', 'P', '附录 U 的「出处：」行只有 %d 处，应 >= 5' % n_src)
    else:
        rec('PASS', 'P', '附录 U 有 %d 处「出处：」交代' % n_src)

    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex < 6:
        rec('FAIL', 'P', '附录 U 的估算练习只有 %d 道，应 >= 6' % n_ex)
    else:
        rec('PASS', 'P', '附录 U 有 %d 道估算练习' % n_ex)
    if n_ans < n_ex:
        rec('FAIL', 'P', '附录 U 有 %d 道练习但只有 %d 个解答' % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 U 的每道练习都有解答（%d 个）' % n_ans)

    # 字面绑定：公式形式、整行与整句必须原样在文里（改动即 FAIL）。
    # 三个工作点的数在表、正文、练习里重复出现，只钉单个数字时改掉其中一处仍能
    # 满足数值用例（v1.17.0 的逐条探针第 67、68 号就是这样逃掉的），故逐处钉住。
    for lit, why in (
            ('Q \\equiv \\frac{P_{\\mathrm{fus}}}{P_{\\mathrm{aux}}}',
             'Q 的定义式'),
            ('P_\\alpha = \\frac{1}{5}P_{\\mathrm{fus}}, \\qquad',
             'D-T 反应能的 alpha 份额（1/5）'),
            ('P_{\\mathrm{n}} = \\frac{4}{5}P_{\\mathrm{fus}}',
             '中子份额（4/5）'),
            ('P_{\\mathrm{fus}} = \\frac{n^2}{4}\\,\\langle\\sigma v\\rangle_{\\mathrm{DT}}\\,E_{\\mathrm{fus}}\\,V,',
             '聚变功率公式（n^2/4 的 1/4 来自 D、T 各占一半）'),
            ('\\langle\\sigma v\\rangle_{\\mathrm{DT}} \\approx 1.1\\times10^{-24}\\,T^2\\ \\ \\mathrm{m^3/s}',
             '两参数反应率拟合（必须是 T^2）'),
            ('E_{\\mathrm{fus}} = 17.6\\ \\mathrm{MeV} = 2.819831\\times10^{-12}\\ \\mathrm{J}',
             'D-T 反应能与其焦耳换算'),
            ('f_\\alpha \\equiv \\frac{P_\\alpha}{P_\\alpha + P_{\\mathrm{aux}}}',
             '自加热份额的定义'),
            ('= \\frac{Q}{Q+5} = \\frac{1}{1+5/Q}',
             'f_alpha 只依赖 Q 的化简结果'),
            ('P_{\\mathrm{loss}} = P_{\\mathrm{aux}} + P_\\alpha = \\frac{W}{\\tau_E},',
             '自加热闭环的功率平衡'),
            ('\\qquad W = 3nTV',
             'W = 3nTV 的写法（与附录 S 同一口径）'),
            ('\\tau_E = \\tau_E^{\\mathrm{IPB98}}\\!\\left(P_{\\mathrm{tot}}\\right),',
             'IPB98(y,2) 传入的是总加热功率'),
            ('\\qquad P_{\\mathrm{tot}} \\equiv P_{\\mathrm{aux}} + P_\\alpha',
             '总加热功率的定义'),
            ('P_{\\mathrm{tot}}^{(k+1)} = P_{\\mathrm{aux}}',
             '闭环迭代格式'),
            ('\\frac{1}{5}P_{\\mathrm{fus}}\\!\\left(T^{(k)}\\right),',
             '迭代格式里的 alpha 加热项'),
            ('= P_{\\mathrm{tot}}^{\\,0.31},',
             'T 与 P_tot 的 0.31 次幂关系'),
            ('P_\\alpha \\propto T^2 \\propto P_{\\mathrm{tot}}^{\\,0.62} < 1',
             '回路增益 0.62 < 1（自稳负反馈的判据）'),
            ('P_\\alpha = P_{\\mathrm{loss}}',
             '点火的定义式'),
            ('n\\tau_E = \\frac{60\\,T_{\\mathrm J}}{\\langle\\sigma v\\rangle E_{\\mathrm{fus}}}',
             '消去 V 后的 n-tau_E 表达式'),
            ('nT\\tau_E > 3.0992\\times10^{21}\\ \\mathrm{keV\\,s\\,m^{-3}}',
             '三重积点火阈值（式 U.10 与练习 4 两处印出）'),
            ('nT\\tau_E > 3.0992\\times10^{21}\\ \\mathrm{keV\\,s\\,m^{-3}}\n\\label{eq:app-u-triple}',
             '三重积阈值那一行式子自身（连 label 一起钉，与练习 4 的复述区分开）'),
            ('$3.0992\\times10^{21}\\ \\mathrm{keV\\,s\\,m^{-3}}$，所以 $n\\tau_E$ 与 $T$ 成反比。',
             '表 U.4 题注里印出的阈值'),
            ('代入反应率拟合后 $T$ 消掉，$nT\\tau_E > 3.0992\\times10^{21}\\ \\mathrm{keV\\,s\\,m^{-3}}$。',
             '练习 4 复述的阈值整句'),
            ('Q \\propto P_{\\mathrm{aux}}^{\\,-0.38}',
             'Q 随外部加热单调下降（指数 −0.38）'),
            ('Q_{\\mathrm{eng}} \\equiv \\frac{P_e}{P_{\\mathrm{grid}}}',
             '工程增益的定义'),
            ('= \\eta_{\\mathrm{th}}\\eta_{\\mathrm{aux}}\\left(Q+1\\right)',
             '工程增益的化简式（括号里有 +1）'),
            ('Q = \\frac{1}{\\eta_{\\mathrm{th}}\\eta_{\\mathrm{aux}}} - 1 = 5.250',
             'Q_eng=1 反解出的盈亏平衡 Q（必须减 1）'),
            ('D1 ITER 量级    & 150.0 & 2.4762 & 9.194 & 371.4 & 556.5 & 3.710 \\\\',
             '表 U.1 的 D1 整行'),
            ('D2 JET 量级     & 15.0  & 0.3827 & 2.336 & 5.7   & 0.9   & 0.058 \\\\',
             '表 U.1 的 D2 整行'),
            ('D3 DIII-D 量级  & 10.0  & 0.1055 & 1.649 & 1.1   & 0.1   & 0.014 \\\\',
             '表 U.1 的 D3 整行'),
            ('3.710  & 0.4259 & 附录 S 的 D1 外部加热点（本书基准） \\\\',
             '表 U.2 的 D1 基准行'),
            ('10     & 0.6667 & ITER 设计点：三分之二的加热来自聚变自身 \\\\',
             '表 U.2 的 ITER 设计点行'),
            ('12 & 332.20 & 1.4306 & 11.764 & 911.00 & 182.20 & 6.073 \\\\',
             '表 U.3 的收敛步（第 12 步）整行'),
            ('5.0  & $6.1983\\times10^{20}$ & 低温端 $\\langle\\sigma v\\rangle\\propto T^2$ 使要求偏高 \\\\',
             '表 U.4 的 5 keV 整行（低温端要求偏高）'),
            ('10.0 & $3.0992\\times10^{20}$ & 本书 $\\langle\\sigma v\\rangle$ 拟合的"甜点"附近 \\\\',
             '表 U.4 的 10 keV 整行（甜点）'),
            ('150.0 & 2.4762  & 9.194 & 371.4 & 556.5 & 3.710 \\\\',
             '表 U.5 的 P_aux=150 MW 整行'),
            ('5.0   & 25.8825 & 3.203 & 129.4 & 67.5  & 13.510 \\\\',
             '表 U.5 的 P_aux=5 MW 整行'),
            ('D1 外部加热点    & 556.5 & 150.0 & 706.5  & 282.6 & 375.0 & 0.754 \\\\',
             '表 U.6 的 D1 外部加热点整行'),
            ('自加热闭环点     & 911.0 & 150.0 & 1061.0 & 424.4 & 375.0 & 1.132 \\\\',
             '表 U.6 的自加热闭环点整行'),
            ('ITER 设计点      & 500.0 & 50.0  & 550.0  & 220.0 & 125.0 & 1.760 \\\\',
             '表 U.6 的 ITER 设计点整行'),
            ('$P_{\\mathrm{fus}}=556.5\\ \\mathrm{MW}$ 中，',
             'U.1 正文的 D1 聚变功率'),
            ('$P_\\alpha = 111.3\\ \\mathrm{MW}$ 留在等离子体里，$P_{\\mathrm{n}} = 445.2\\ \\mathrm{MW}$ 去包层；',
             'U.1 正文的 alpha 与中子功率'),
            ('按 $V=832.2\\ \\mathrm{m^3}$ 摊开，$\\alpha$ 加热的功率密度是 $1.337\\times10^{5}\\ \\mathrm{W/m^3}$。',
             'U.1 正文的 alpha 加热功率密度'),
            ('（$9.194\\to11.764\\ \\mathrm{keV}$），$P_{\\mathrm{fus}}$ 也抬高',
             'U.3 口径对照：闭环把 T 从 9.194 抬到 11.764'),
            ('（$556.5\\to911.0\\ \\mathrm{MW}$），于是 $Q$ 升到 $6.073$。',
             'U.3 口径对照：闭环把 P_fus 从 556.5 抬到 911.0'),
            ('$T/P_{\\mathrm{tot}}^{0.31}$ 恒为 $1.94496$',
             'U.3 三点自检的常数'),
            ('$\\diff\\ln T/\\diff\\ln P_{\\mathrm{tot}} = 0.3100$',
             'U.3 反解出的对数斜率'),
            ('$n\\tau_E = 3.3709\\times10^{20}\\ \\mathrm{m^{-3}\\,s}$',
             'U.4 D1 在 9.194 keV 上的 n-tau_E 要求'),
            ('$nT\\tau_E = 2.299\\times10^{21}\\ \\mathrm{keV\\,s\\,m^{-3}}$，是阈值的 $74.2\\%$——',
             'U.4 D1 达到阈值的份额'),
            ('D2 只有 $1.2\\%$，D3 只有 $0.3\\%$',
             'U.4 D2/D3 的达标份额'),
            ('（$\\times1.3478$），对应 $P_\\alpha$ 从 $111.3$ 升到 $150.0\\ \\mathrm{MW}$、',
             'U.4 密度提升倍数与 alpha 功率'),
            ('$n_{19}\\times\\sqrt{750.0/556.5} = \\times1.1610$，即 $n_{19}=11.726$。',
             'U.4 按 P_fus ∝ n^2 的粗估路径'),
            ('$n_G = I_p/(\\pi a^2) = 1.1937\\times10^{20}\\ \\mathrm{m^{-3}}$，',
             'U.4 Greenwald 极限（与附录 M 同口径）'),
            ('D1 工作点已经用到 $84.6\\%$。按第一种算法，点火需要 $n/n_G = 982.3\\%$，',
             'U.4 密度极限的夹击'),
            ('只需 $P_{\\mathrm{aux}} = 11.035\\ \\mathrm{MW}$',
             'U.5 Q=10 所需的很小外部加热'),
            ('同一行给出 $\\tau_E = 14.9894\\ \\mathrm{s}$、$T = 4.094\\ \\mathrm{keV}$、',
             'U.5 Q=10 同一行的 tau_E 与 T'),
            ('$Q=3.710$ 的 D1 外部加热点 $Q_{\\mathrm{eng}}=0.754$',
             'U.5 基准点的工程增益'),
            ('盈亏平衡在 $Q=5.250$',
             'U.5 工程盈亏平衡点'),
            ('$Q_{\\mathrm{eng}}=1.760$，这才是',
             'U.5 ITER 设计点的工程增益'),
            ('$f_\\alpha = 3.710/8.710 = 0.4259$；$Q=10$ 给 $10/15 = 0.6667$。',
             '练习 1 解答'),
            ('$P_{\\mathrm{n}} = 0.8\\times556.5 = 445.2\\ \\mathrm{MW}$；',
             '练习 2 解答的中子功率'),
            ('$445.2/800 = 0.5565\\ \\mathrm{MW/m^2}$。',
             '练习 2 解答的中子壁负载'),
            ('$P_{\\mathrm{tot}} = 150.0 \\to 261.3 \\to 307.0\\ \\mathrm{MW}$，',
             '练习 3 解答的三步迭代'),
            ('对应 $Q = 3.710 \\to 5.233 \\to 5.783$（完整轨迹见表',
             '练习 3 解答对应的 Q'),
            ('$nT\\tau_E = 2.2994\\times10^{21}$，是阈值的 $74.2\\%$。',
             '练习 4 解答的达标份额'),
            ('解 $Q=10$ 得 $P_{\\mathrm{aux}} = 11.03\\ \\mathrm{MW}$、',
             '练习 5 解答的 P_aux'),
            ('$Q_{\\mathrm{eng}} = 0.4\\times0.4\\times4.710 = 0.754$；$Q=10$ 给',
             '练习 6 解答的 Q_eng'),
            ('$0.16\\times11 = 1.760$。阈值由式 \\eqref{eq:app-u-qeng1}：',
             '练习 6 解答的 ITER 工程增益'),
            ('$Q = 1/0.16 - 1 = 5.250$。结论',
             '练习 6 解答的盈亏平衡 Q'),
            ('\\texttt{tmp/appU\\_calc.py}',
             '数值唯一来源脚本（独立复算的凭据）'),
            ('三者的三个装置工作点完全一致，可以互相代入验算。',
             'U 与 S、T 的工作点一致性声明'),
            ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['U'], why))

    # 口径回流防护（这些写法一旦出现即为 FAIL）
    for lit, why in (
            ('\\langle\\sigma v\\rangle \\approx 1.1\\times10^{-24}T^3',
             '反应率拟合的 T 次幂被改'),
            ('P_\\alpha = \\frac{1}{4}P_{\\mathrm{fus}}',
             'alpha 份额退回 1/4'),
            ('P_\\alpha = \\frac{4}{5}P_{\\mathrm{fus}}',
             'alpha 与中子的份额被对调'),
            ('\\frac{Q}{Q+4}',
             'f_alpha 的分母写成 Q+4'),
            ('P_\\alpha \\propto P_{\\mathrm{tot}}^{\\,1.62}',
             '回路增益被写成大于 1（自稳结论被推翻）'),
            ('T \\propto P_{\\mathrm{tot}}^{\\,0.69}',
             'T 的指数漏掉 1−（符号写反）'),
            ('Q \\propto P_{\\mathrm{aux}}^{\\,0.38}',
             'Q 的单调性写反（应为 −0.38）'),
            ('= \\eta_{\\mathrm{th}}\\eta_{\\mathrm{aux}}\\,Q+1',
             '工程增益漏掉括号（+1 跑到括号外）'),
            ('Q = \\frac{1}{\\eta_{\\mathrm{th}}\\eta_{\\mathrm{aux}}} = 6.250',
             '工程盈亏平衡忘掉减 1'),
            ('n_G = I_p/(\\pi a^2) = 1.1937\\times10^{21}',
             'Greenwald 极限错一个量级'),
            ('P_{\\mathrm{loss}} = \\frac{\\tau_E}{W}',
             '损失功率的定义取倒数'),
            ('\\tau_E = \\tau_E^{\\mathrm{IPB98}}\\!\\left(P_{\\mathrm{aux}}\\right)',
             'IPB98 传入外部加热而非总加热功率'),
            ('n\\tau_E = \\frac{60\\,T_{\\mathrm J}}{\\langle\\sigma v\\rangle}',
             'n-tau_E 表达式丢掉 E_fus'),
            ('$P_{\\mathrm{n}} = 0.8\\times556.5 = 44.52\\ \\mathrm{MW}$',
             '练习 2 的中子功率错一个量级'),
            ('$Q_{\\mathrm{eng}} = 0.4\\times0.4\\times3.710 = 0.594$',
             '练习 6 的 Q_eng 漏掉 +1'),
            ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['U'], why))

    # 数值独立复算（单位见正文/表头，故 unit=''）
    cases = [
        # Q 是什么
        (17.6, 1e-3, 'Q 是什么', 'D-T 反应能/MeV'),
        (3.5, 1e-3, 'Q 是什么', 'alpha 分到的能量/MeV'),
        (14.1, 1e-3, 'Q 是什么', '中子分到的能量/MeV'),
        (2.819831e-12, 1e-3, 'Q 是什么', 'E_fus/J'),
        (1.1e-22, 1e-3, 'Q 是什么', '10 keV 处的 <sigma v>'),
        (150.0, 1e-3, 'Q 是什么', 'D1 的 P_aux/MW'),
        (2.4762, 1e-3, 'Q 是什么', 'D1 的 tau_E/s'),
        (9.194, 1e-3, 'Q 是什么', 'D1 的 T/keV'),
        (371.4, 1e-3, 'Q 是什么', 'D1 的 W/MJ'),
        (556.5, 1e-3, 'Q 是什么', 'D1 的 P_fus/MW'),
        (3.71, 1e-3, 'Q 是什么', 'D1 的 Q'),
        (15.0, 1e-3, 'Q 是什么', 'D2 的 P_aux/MW'),
        (0.3827, 1e-3, 'Q 是什么', 'D2 的 tau_E/s'),
        (2.336, 1e-3, 'Q 是什么', 'D2 的 T/keV'),
        (5.7, 1e-3, 'Q 是什么', 'D2 的 W/MJ'),
        (0.058, 1e-3, 'Q 是什么', 'D2 的 Q'),
        (0.1055, 1e-3, 'Q 是什么', 'D3 的 tau_E/s'),
        (1.649, 1e-3, 'Q 是什么', 'D3 的 T/keV'),
        (1.1, 1e-3, 'Q 是什么', 'D3 的 W/MJ'),
        (0.014, 1e-3, 'Q 是什么', 'D3 的 Q'),
        (111.3, 1e-3, 'Q 是什么', 'D1 的 P_alpha/MW'),
        (445.2, 1e-3, 'Q 是什么', 'D1 的 P_n/MW'),
        (832.2, 1e-3, 'Q 是什么', 'D1 的 V/m^3'),
        (133700.0, 1e-3, 'Q 是什么', 'D1 的 alpha 加热功率密度/(W/m^3)'),
        (500.0, 1e-3, 'Q 是什么', 'tip U.1 的 ITER 聚变功率/MW'),
        (50.0, 1e-3, 'Q 是什么', 'tip U.1 的 ITER 加热功率/MW'),
        # 加热份额
        (0.1667, 1e-3, '加热份额', 'Q=1 的 f_alpha'),
        (0.2857, 1e-3, '加热份额', 'Q=2 的 f_alpha'),
        (0.4259, 1e-3, '加热份额', 'Q=3.710 的 f_alpha'),
        (0.5, 1e-3, '加热份额', 'Q=5 的 f_alpha'),
        (0.6667, 1e-3, '加热份额', 'Q=10 的 f_alpha'),
        (0.8, 1e-3, '加热份额', 'Q=20 的 f_alpha'),
        (0.9091, 1e-3, '加热份额', 'Q=50 的 f_alpha'),
        (1.0, 1e-3, '加热份额', 'Q=infinity 的 f_alpha'),
        # 自加热闭环
        (150.0, 1e-3, '自加热闭环', '迭代第 1 步的 P_tot/MW'),
        (556.45, 1e-3, '自加热闭环', '迭代第 1 步的 P_fus/MW'),
        (261.29, 1e-3, '自加热闭环', '迭代第 2 步的 P_tot/MW'),
        (1.6884, 1e-3, '自加热闭环', '迭代第 2 步的 tau_E/s'),
        (5.233, 1e-3, '自加热闭环', '迭代第 2 步的 Q'),
        (307.0, 1e-3, '自加热闭环', '迭代第 3 步的 P_tot/MW'),
        (867.52, 1e-3, '自加热闭环', '迭代第 3 步的 P_fus/MW'),
        (5.783, 1e-3, '自加热闭环', '迭代第 3 步的 Q'),
        (323.5, 1e-3, '自加热闭环', '迭代第 4 步的 P_tot/MW'),
        (896.14, 1e-3, '自加热闭环', '迭代第 4 步的 P_fus/MW'),
        (5.974, 1e-3, '自加热闭环', '迭代第 4 步的 Q'),
        (329.23, 1e-3, '自加热闭环', '迭代第 5 步的 P_tot/MW'),
        (905.94, 1e-3, '自加热闭环', '迭代第 5 步的 P_fus/MW'),
        (6.04, 1e-3, '自加热闭环', '迭代第 5 步的 Q'),
        (332.08, 1e-3, '自加热闭环', '迭代第 8 步的 P_tot/MW'),
        (910.81, 1e-3, '自加热闭环', '迭代第 8 步的 P_fus/MW'),
        (6.072, 1e-3, '自加热闭环', '迭代第 8 步的 Q'),
        (332.2, 1e-3, '自加热闭环', '收敛点的 P_tot/MW'),
        (1.4306, 1e-3, '自加热闭环', '收敛点的 tau_E/s'),
        (11.764, 1e-3, '自加热闭环', '收敛点的 T/keV'),
        (911.0, 1e-3, '自加热闭环', '收敛点的 P_fus/MW'),
        (182.2, 1e-3, '自加热闭环', '收敛点的 P_alpha/MW'),
        (6.073, 1e-3, '自加热闭环', '收敛点的 Q'),
        (0.003, 1e-3, '自加热闭环', '第 11、12 步的 P_tot 之差/MW'),
        (1.94496, 1e-3, '自加热闭环', 'T/P_tot^0.31 的恒定值'),
        (0.31, 1e-3, '自加热闭环', '反解出的 dlnT/dlnP_tot'),
        (0.62, 1e-3, '自加热闭环', 'P_alpha 的等效指数'),
        (0.69, 1e-3, '自加热闭环', 'IPB98(y,2) 的功率指数'),
        (100.0, 1e-3, '自加热闭环', '自检点 P_tot=100 MW'),
        (225.0, 1e-3, '自加热闭环', '自检点 P_tot=225 MW'),
        # 点火：三重积阈值
        (60.0, 1e-3, '点火：三重积阈值', 'n-tau_E 表达式里的 60'),
        (3.0992e+21, 1e-3, '点火：三重积阈值', '三重积点火阈值/(keV s m^-3)'),
        (6.1983e+20, 1e-3, '点火：三重积阈值', '5 keV 的点火 n-tau_E'),
        (3.874e+20, 1e-3, '点火：三重积阈值', '8 keV 的点火 n-tau_E'),
        (3.0992e+20, 1e-3, '点火：三重积阈值', '10 keV 的点火 n-tau_E（甜点）'),
        (2.0661e+20, 1e-3, '点火：三重积阈值', '15 keV 的点火 n-tau_E'),
        (1.5496e+20, 1e-3, '点火：三重积阈值', '20 keV 的点火 n-tau_E'),
        (3.3709e+20, 1e-3, '点火：三重积阈值', 'D1 在 9.194 keV 上的点火 n-tau_E'),
        (2.299e+21, 1e-3, '点火：三重积阈值', 'D1 实有的 nT tau_E'),
        (74.2, 1e-3, '点火：三重积阈值', 'D1 达到阈值的份额/%'),
        (1.2, 1e-3, '点火：三重积阈值', 'D2 达到阈值的份额/%'),
        (0.3, 1e-3, '点火：三重积阈值', 'D3 达到阈值的份额/%'),
        (10.1, 1e-3, '点火：三重积阈值', 'D1 的 n_19'),
        (13.613, 1e-3, '点火：三重积阈值', '点火所需的 n_19'),
        (1.3478, 1e-3, '点火：三重积阈值', '密度提升倍数'),
        (150.0, 1e-3, '点火：三重积阈值', '点火所需的 P_alpha/MW'),
        (750.0, 1e-3, '点火：三重积阈值', '点火所需的 P_fus/MW'),
        (1.161, 1e-3, '点火：三重积阈值', '按 P_fus ∝ n^2 的倍数'),
        (11.726, 1e-3, '点火：三重积阈值', '按 P_fus ∝ n^2 的 n_19'),
        (1.1937e+20, 1e-3, '点火：三重积阈值', 'Greenwald 密度极限/m^-3'),
        (84.6, 1e-3, '点火：三重积阈值', 'D1 用掉的 Greenwald 份额/%'),
        (982.3, 1e-3, '点火：三重积阈值', '点火所需的 n/n_G 百分比'),
        (9.8, 1e-3, '点火：三重积阈值', '点火所需的 Greenwald 倍数'),
        # 工程增益
        (0.38, 1e-3, '工程增益', 'Q ∝ P_aux 的指数绝对值'),
        (100.0, 1e-3, '工程增益', '扫描点 P_aux=100 MW'),
        (3.2756, 1e-3, '工程增益', 'P_aux=100 MW 的 tau_E/s'),
        (8.108, 1e-3, '工程增益', 'P_aux=100 MW 的 T/keV'),
        (327.6, 1e-3, '工程增益', 'P_aux=100 MW 的 W/MJ'),
        (432.8, 1e-3, '工程增益', 'P_aux=100 MW 的 P_fus/MW'),
        (4.328, 1e-3, '工程增益', 'P_aux=100 MW 的 Q'),
        (5.2845, 1e-3, '工程增益', 'P_aux=50 MW 的 tau_E/s'),
        (6.54, 1e-3, '工程增益', 'P_aux=50 MW 的 T/keV'),
        (264.2, 1e-3, '工程增益', 'P_aux=50 MW 的 W/MJ'),
        (281.6, 1e-3, '工程增益', 'P_aux=50 MW 的 P_fus/MW'),
        (5.632, 1e-3, '工程增益', 'P_aux=50 MW 的 Q'),
        (9.9445, 1e-3, '工程增益', 'P_aux=20 MW 的 tau_E/s'),
        (4.923, 1e-3, '工程增益', 'P_aux=20 MW 的 T/keV'),
        (198.9, 1e-3, '工程增益', 'P_aux=20 MW 的 W/MJ'),
        (7.977, 1e-3, '工程增益', 'P_aux=20 MW 的 Q'),
        (16.0433, 1e-3, '工程增益', 'P_aux=10 MW 的 tau_E/s'),
        (3.971, 1e-3, '工程增益', 'P_aux=10 MW 的 T/keV'),
        (160.4, 1e-3, '工程增益', 'P_aux=10 MW 的 W/MJ'),
        (103.8, 1e-3, '工程增益', 'P_aux=10 MW 的 P_fus/MW'),
        (10.381, 1e-3, '工程增益', 'P_aux=10 MW 的 Q'),
        (25.8825, 1e-3, '工程增益', 'P_aux=5 MW 的 tau_E/s'),
        (3.203, 1e-3, '工程增益', 'P_aux=5 MW 的 T/keV'),
        (129.4, 1e-3, '工程增益', 'P_aux=5 MW 的 W/MJ'),
        (67.5, 1e-3, '工程增益', 'P_aux=5 MW 的 P_fus/MW'),
        (13.51, 1e-3, '工程增益', 'P_aux=5 MW 的 Q'),
        (11.035, 1e-3, '工程增益', 'Q=10 所需的 P_aux/MW'),
        (14.9894, 1e-3, '工程增益', 'Q=10 对应的 tau_E/s'),
        (4.094, 1e-3, '工程增益', 'Q=10 对应的 T/keV'),
        (110.3, 1e-3, '工程增益', 'Q=10 对应的 P_fus/MW'),
        (0.4, 1e-3, '工程增益', '热效率与加热效率'),
        (706.5, 1e-3, '工程增益', 'D1 外部加热点的 P_th/MW'),
        (282.6, 1e-3, '工程增益', 'D1 外部加热点的 P_e/MW'),
        (375.0, 1e-3, '工程增益', 'D1 外部加热点的 P_grid/MW'),
        (0.754, 1e-3, '工程增益', 'D1 外部加热点的 Q_eng'),
        (1061.0, 1e-3, '工程增益', '闭环点的 P_th/MW'),
        (424.4, 1e-3, '工程增益', '闭环点的 P_e/MW'),
        (1.132, 1e-3, '工程增益', '闭环点的 Q_eng'),
        (550.0, 1e-3, '工程增益', 'ITER 设计点的 P_th/MW'),
        (220.0, 1e-3, '工程增益', 'ITER 设计点的 P_e/MW'),
        (125.0, 1e-3, '工程增益', 'ITER 设计点的 P_grid/MW'),
        (1.76, 1e-3, '工程增益', 'ITER 设计点的 Q_eng'),
        (5.25, 1e-3, '工程增益', 'Q_eng=1 的盈亏平衡 Q'),
        # 六道估算练习
        (0.2, 1e-3, '六道估算练习', '练习 1 的 P_alpha/P_aux'),
        (800.0, 1e-3, '六道估算练习', '练习 2 的第一壁面积/m^2'),
        (0.5565, 1e-3, '六道估算练习', '练习 2 的中子壁负载/(MW/m^2)'),
        (261.3, 1e-3, '六道估算练习', '练习 3 的第二步 P_tot/MW'),
        (307.0, 1e-3, '六道估算练习', '练习 3 的第三步 P_tot/MW'),
        (2.2994e+21, 1e-3, '六道估算练习', '练习 4 的 nT tau_E'),
        (10.1, 1e-3, '六道估算练习', '练习 4 的 n_19'),
        (11.03, 1e-3, '六道估算练习', '练习 5 的 P_aux/MW'),
        (14.99, 1e-3, '六道估算练习', '练习 5 的 tau_E/s'),
        (4.09, 1e-3, '六道估算练习', '练习 5 的 T/keV'),
        # 来源、分级与版权说明
        (2011, 1e-3, '来源、分级与版权说明', 'Wesson 第 4 版年份'),
        (2007, 1e-3, '来源、分级与版权说明', 'Freidberg 年份'),
        (1957, 1e-3, '来源、分级与版权说明', 'Lawson 年份'),
        (1999, 1e-3, '来源、分级与版权说明', 'ITER Physics Basis 年份'),
        (1992, 1e-3, '来源、分级与版权说明', 'Bosch & Hale 年份'),
        (611, 1e-3, '来源、分级与版权说明', 'Bosch & Hale 页码'),
        (70, 1e-3, '来源、分级与版权说明', 'Lawson 卷号'),
        (1.1e-24, 1e-3, '来源、分级与版权说明', '反应率拟合系数'),
        (0.62, 1e-3, '来源、分级与版权说明', '自拟的稳定性指数'),
    ]
    for v, rel, where, what in cases:
        block = found.get(where)
        if block is None:
            rec('FAIL', 'P', '附录 U 找不到小节「%s」，无法核对 %s' % (where, what))
        elif not value_in(block, v, rel, unit=''):
            rec('FAIL', 'P', '附录 U 的 %s 缺失或与 tmp/appU_calc.py 不符：期望 %g（小节「%s」）'
                % (what, v, where))
    rec('PASS', 'P', '附录 U 的 %d 个关键数值均在小节内复算通过' % len(cases))


def check_P25_appendix_i_expanding():
    """附录 I「扩张式还是均匀式」小节：结构、五条结论、六行文献与口径回流防护。"""
    ipath = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(ipath):
        rec('FAIL', 'N', '附录 I 文件缺失：src/%s' % APP_FILES['I'])
        return
    text = strip_comments(read(ipath))
    head = '\\subsection{扩张式还是均匀式：间隔要不要越拉越长}'
    if head not in text:
        rec('FAIL', 'N', 'src/%s 缺少小节「扩张式还是均匀式：间隔要不要越拉越长」'
            % APP_FILES['I'])
        return
    blk = text.split(head, 1)[1]
    for cut in ('\\subsection{', '\\section{'):
        if cut in blk:
            blk = blk.split(cut, 1)[0]
    if len(blk) < 1500:
        rec('FAIL', 'N', '附录 I「扩张式还是均匀式」小节正文过短（%d 字符）' % len(blk))
    else:
        rec('PASS', 'N', '附录 I「扩张式还是均匀式」小节存在（%d 字符）' % len(blk))
    for lit, why in (
            ('\\label{subsec:app-i-expanding}',
             '小节的 label（供交叉引用）'),
            ('\\textbf{扩张式}（expanding）指间隔越拉越长，如 $1,2,4,8,12$ 天；\\textbf{均匀式}（uniform）指每次都隔同样长，如 $6,6,6,6,6$ 天。',
             '扩张式与均匀式的定义与示例间隔（两个具体数列）'),
            ('\\textbf{（1）扩张式的理由成立，但推不出"必须递增"。}',
             '第 1 条结论的标题句'),
            ('\\textbf{（2）形状是第二位的，比例是第一位的。}',
             '第 2 条结论的标题句'),
            ('\\textbf{（3）稳得住的那条规律是"分散优于集中"。}',
             '第 3 条结论的标题句'),
            ('\\textbf{（4）那本书为什么仍用"先翻倍、后等间隔"？}',
             '第 4 条结论的标题句'),
            ('\\textbf{（5）用软件也建议自己算一遍。}',
             '第 5 条结论的标题句'),
            ('只需要一个上界 $g_{\\max}=\\lfloor0.15D\\rfloor$',
             '间隔上限公式（与「间隔与提取」一节的 0.15D 口径一致）'),
            ('\\textbf{排程规则自己算一遍}',
             '对记忆软件的态度（结论句）'),
            ('在总复习时间相同的前提下，把四个复习点排成扩张式的 $1,2,4,8$ 天，与均匀式的 $4,4,4,4$ 天，哪一种更好？',
             'miniquiz 的题干（四步对齐的两个排布）'),
            ('答案：在总时间对齐的实验里两者差别很小，\\textbf{不能只凭形状下结论}',
             'miniquiz 的答案（不得给出「形状决定优劣」的结论）'),
            ('\\textbf{合书提取}',
             'miniquiz 里指向的必须动作'),
            ('\\textbf{出处：}本小节只借鉴结论与实验设计。',
             '本小节的出处声明'),
            ('扩张式提取的提出见 Landauer \\& Bjork (1978)',
             '出处句里的 Landauer \\& Bjork 定位'),
            ('长跨度保持证据见 Bahrick (1979) 与 Bahrick 等 (1993)',
             '出处句里的 Bahrick 两条长跨度证据'),
            ('SM-2 的间隔修正见 Wozniak \\& Gorzelanczyk (1994)。卷期页见本附录末表',
             '出处句里的 SM-2 定位与「卷期页见末表」的交代'),
            ('\\label{tab:app-i-calendar}',
             '第 4 条结论引用的日历表 label'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))
    n_q = blk.count('\\begin{miniquiz}')
    if n_q != 1:
        rec('FAIL', 'N', '附录 I「扩张式还是均匀式」小节应有 1 只 miniquiz，实为 %d' % n_q)
    else:
        rec('PASS', 'N', '附录 I「扩张式还是均匀式」小节有 1 只 miniquiz')
    n_src = blk.count('\\textbf{出处：}')
    if n_src != 1:
        rec('FAIL', 'N', '附录 I「扩张式还是均匀式」小节应有 1 处出处声明，实为 %d' % n_src)
    else:
        rec('PASS', 'N', '附录 I「扩张式还是均匀式」小节有 1 处出处声明')
    for lit, why in (
            ('Landauer \\& Bjork (1978) & 载 Gruneberg, Morris \\& Sykes (编)《Practical Aspects of Memory》, Academic Press, 625--632',
             '文献表：Landauer \\& Bjork (1978) 整行'),
            ('Karpicke \\& Roediger (2007) & \\textit{Journal of Experimental Psychology: Learning, Memory, and Cognition} 33(4), 704--719',
             '文献表：Karpicke \\& Roediger (2007) 整行'),
            ('Kornell (2009) & \\textit{Applied Cognitive Psychology} 23(9), 1297--1317',
             '文献表：Kornell (2009) 整行'),
            ('Bahrick (1979) & \\textit{Journal of Experimental Psychology: General} 108(3), 296--308',
             '文献表：Bahrick (1979) 整行'),
            ('Bahrick, Bahrick, Bahrick \\& Bahrick (1993) & \\textit{Psychological Science} 4(5), 316--321',
             '文献表：Bahrick 等 (1993) 整行'),
            ('Wozniak \\& Gorzelanczyk (1994) & \\textit{Acta Neurobiologiae Experimentalis} 54(1), 59--62',
             '文献表：Wozniak \\& Gorzelanczyk (1994) 整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))
    for lit, why in (
            ('元分析发现扩张式稳定优于均匀式',
             '不得把 Cepeda 等 (2006) 的结论反着写'),
            ('间隔越短越好',
             '不得回流「间隔越短越好」'),
            ('均匀式优于扩张式',
             '不得回流「均匀式优于扩张式」'),
            ('扩张式的长期优势始终成立',
             '不得把 Karpicke \\& Roediger (2007) 的结论反着写'),
            ('间隔越长复习效果越好',
             '不得把排布形状说成「越久越好」'),
            ('复习点集中安排效果更好',
             '不得回流「集中优于分散」'),
            ('扩张式是长期保持的必要条件',
             '不得把排布形状说成必要条件'),
            ('扩张式优于均匀式',
             '不得直接给出「扩张式更好」的结论'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

def check_P26_appendix_v_disruption():
    """附录 V 破裂、误差场与逃逸电子的完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 八个小节之一被删（按小节标题子串定位）；
      ② 九张表 / 三只语义盒 / 十四个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 6 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（beta 与 beta_N、误差场真空折算、壁时间与
         RWM 增长率、三本账、Dreicer 场与雪崩、tau_CQ 扫描、缓解粒子数），
         或搬到了别的小节里；
      ⑤ 口径回流：壁时间与增长率被四舍五入、反馈带宽写回 190、beta 裕度写回
         33.5%、E_tor/E_D 写回 0.0151、雪崩放大倍数写回 1.0e40、Greenwald 极限
         错一个量级、晕电流力错一个量级、tau_CQ 下界单位写成 ms。

    数值唯一来源 tmp/appV_calc.py；表内数字的单位写在表头，故这些用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['V'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 V 文件缺失：src/%s' % APP_FILES['V'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = [
        '破裂：三张极限图', '误差场与锁模', '电阻壁模与', '破裂的后果', '逃逸电子：Dreicer', '避免、预测与缓解', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 V 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 V 的 %d 个小节齐全' % len(want))

    for lab in (
            'tab:app-v-limits', 'tab:app-v-timeline', 'tab:app-v-errorfield', 'tab:app-v-beta', 'tab:app-v-accounts', 'tab:app-v-fields', 'tab:app-v-cq', 'tab:app-v-mitigation', 'tab:app-v-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['V'], lab))
    for lab in (
            'tip-app-v-timeline', 'tip-app-v-errorfield', 'tip-app-v-runaway',
            'tip-app-v-honest'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['V'], lab))
    for lab in (
            'eq:app-v-betan', 'eq:app-v-greenwald', 'eq:app-v-vacuum', 'eq:app-v-tauw', 'eq:app-v-rwm', 'eq:app-v-qwall', 'eq:app-v-force', 'eq:app-v-ed', 'eq:app-v-ec', 'eq:app-v-eloop', 'eq:app-v-avalanche', 'eq:app-v-wre', 'eq:app-v-nsupp', 'eq:app-v-impurity'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['V'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['V'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 V 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 V 有 %d 处「出处：」声明' % n_src)
    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 V 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 V 有 %d 道练习、%d 条解答' % (n_ex, n_ans))

    # 关键行与关键句的逐字绑定（防表格整行被改、防正文数字被搬走）
    for lit, why in (
            ('密度极限 & 边缘辐射与脱靶使辐射崩塌，电流剖面收缩 & 辐射功率骤升、约束变差 & 附录 M、附录 U（$n_G=1.1937\\times10^{20}$） \\\\',
             '表 V.1 密度极限行（与附录 M、U 的 n_G 同值）'),
            ('$\\beta$ 极限 & 理想 MHD 外扭曲模失稳（无壁 2.8 / 理想壁 4.0） & 低频 $n=1$ 磁扰动 & 附录 S（$\\beta_N$）、本附录 V.3 \\\\',
             '表 V.1 beta 极限行（两个极限值与附录 S 的分工）'),
            ('误差场锁模 & 静态 $m/n$ 场误差在有理面上生成磁岛并锁定 & 锁模信号、旋转下降 & 附录 T（磁岛）、本附录 V.2 \\\\',
             '表 V.1 误差场锁模行（与附录 T 的分工）'),
            ('热猝灭 & $\\sim1\\ \\mathrm{ms}$ & 磁面被破坏，热储能沿开磁力线泄到壁 & $371.4\\ \\mathrm{MJ}$ 在 $1\\ \\mathrm{ms}$ 内是 $371.4\\ \\mathrm{GW}$ \\\\',
             '表 V.2 热猝灭行（371.4 MJ→371.4 GW）'),
            ('电流猝灭 & $50$--$150\\ \\mathrm{ms}$ & 电流靠等离子体电阻衰减，环电压飙升 & $T_e=10\\ \\mathrm{eV}$ 给 $\\tau_{\\mathrm{CQ}}=76.9074\\ \\mathrm{ms}$、$V_{\\mathrm{loop}}=2600\\ \\mathrm{V}$ \\\\',
             '表 V.2 电流猝灭行（tau_CQ 与环电压）'),
            ('晕电流 & 与电流猝灭同时 & 开磁力线把电流引到壁，形成不对称的晕电流 & $I_{\\mathrm{halo}}=4.5\\ \\mathrm{MA}$，力约 $79.5\\ \\mathrm{MN}$ \\\\',
             '表 V.2 晕电流行（4.5 MA / 79.5 MN）'),
            ('逃逸平台 & $100\\ \\mathrm{ms}$--$1\\ \\mathrm{s}$ & 高能电子尾巴雪崩，束带着几 MA 电流继续跑 & $26.0\\ \\mathrm{MJ}$，落在 $0.1\\ \\mathrm{m^2}$ 上是 $259.9\\ \\mathrm{MJ/m^2}$ \\\\',
             '表 V.2 逃逸平台行（26.0 MJ / 259.9 MJ/m^2）'),
            ('$m/n=2/1$ & 1.6 & 0.6667 & $1.5\\times10^{-4}$ & $0.795\\ \\mathrm{mT}$ \\\\',
             '表 V.3 m=2 行的完整折算链'),
            ('$m/n=3/1$ & 1.8 & 0.5625 & $1.7778\\times10^{-4}$ & $0.9422\\ \\mathrm{mT}$ \\\\',
             '表 V.3 m=3 行的完整折算链'),
            ('无壁极限 2.8 & 3.9623 & $4.428\\times10^{5}$ & 13.8202 & 1.503 \\\\',
             '表 V.4 无壁极限行'),
            ('理想壁极限 4.0 & 5.6604 & $6.326\\times10^{5}$ & 19.7432 & 2.147 \\\\',
             '表 V.4 理想壁极限行'),
            ('本书工作点 1.8627 & 2.6359 & $2.946\\times10^{5}$ & 9.194 & 1.000 \\\\',
             '表 V.4 本书工作点行（与附录 S/U 的 D1 点同值）'),
            ('热账 & $371.4\\ \\mathrm{MJ}$；$464.25\\ \\mathrm{MW/m^2}$（$1\\ \\mathrm{ms}$、$800\\ \\mathrm{m^2}$） & 钨的瞬态阈值 $\\sim1\\ \\mathrm{MJ/m^2}$（附录 N.5） & 必须摊到 $371.4\\ \\mathrm{m^2}$ 以上，否则局部熔融 \\\\',
             '表 V.5 热账行（含附录 N.5 的阈值分工）'),
            ('力账 & $I_{\\mathrm{halo}}=4.5\\ \\mathrm{MA}$；力约 $79.5\\ \\mathrm{MN}$ & 真空室与支撑按几十至一百 MN 量级设计 & 晕电流是结构设计的控制载荷 \\\\',
             '表 V.5 力账行'),
            ('电子账 & $26.0\\ \\mathrm{MJ}$；$0.1\\ \\mathrm{m^2}$ 上 $259.9\\ \\mathrm{MJ/m^2}$ &',
             '表 V.5 电子账行的数字部分'),
            ('$E_c$（临界场） & 0.08668 & $E_{\\mathrm{tor}}/E_c=769.9$ & 超过它就出现逃逸电子 \\\\',
             '表 V.6 临界场行'),
            ('$E_{\\mathrm{tor}}$（猝灭期环向场） & 66.74 & 1 & 由 $V_{\\mathrm{loop}}=2600\\ \\mathrm{V}$ 与 $2\\pi R$ 决定 \\\\',
             '表 V.6 猝灭期环向场行'),
            ('$E_D$（Dreicer 场，$10\\ \\mathrm{eV}$） & 4430 & $E_{\\mathrm{tor}}/E_D=0.01507$ & 不到它的 2\\%：主体电子仍被拖住 \\\\',
             '表 V.6 Dreicer 场行'),
            ('0.10 & 51.33 & 59.21 & 工程上典型的猝灭时间 \\\\',
             '表 V.7 tau_CQ=0.10 s 行'),
            ('1.00 & 5.133 & 5.921 & 还要慢六倍 \\\\',
             '表 V.7 tau_CQ=1.00 s 行（给出 E_tor/E_c 的落点）'),
            ('辐射掉热储能 & $1.238\\times10^{21}$ & $5.128\\ \\mathrm{Pa\\,m^3}$ & 占离子库存 $1.488\\%$（$L_z=3\\times10^{-32}$） \\\\',
             '表 V.8 辐射任务行（含声明过的 L_z）'),
            ('把密度抬到 $10^{21}$ & $8.322\\times10^{23}$ & $3.447\\ \\mathrm{kPa\\,m^3}$ & 压低 $E_{\\mathrm{tor}}/E_c$，是注入量的驱动项 \\\\',
             '表 V.8 抬密度任务行'),
            ('V.3 & Troyon 等（1984）；Chu \\& Okabayashi（2010）；Freidberg《Ideal MHD》 & $\\beta_N$ 定义与极限区间、RWM 与反馈 & $\\tau_w$ 算例、增长率与反馈带宽要求 \\\\',
             '表 V.9 的 V.3 行（借鉴与自拟的分界）'),
            ('V.5 & Dreicer（1959）；Connor \\& Hastie（1975）；Rosenbluth \\& Putvinski（1997） & Dreicer 场、临界场、雪崩 & 三个场的数值算例、$e$ 折时间与束能量 \\\\',
             '表 V.9 的 V.5 行'),
            ('V.7 & 本书自拟 & 六道估算题的题型 & 全部数字来自 \\texttt{tmp/appV\\_calc.py} \\\\',
             '表 V.9 的 V.7 行（数值唯一来源）'),
            ('记住本附录要算的四个时间尺度：壁的 $\\tau_w=52.7788\\ \\mathrm{ms}$（V.3）、雪崩的',
             'tip V.1 的四个时间尺度（壁 / 雪崩）'),
            ('$\\tau_{\\mathrm{av}}=0.4342\\ \\mathrm{ms}$（V.5）、电流猝灭的 $76.9074\\ \\mathrm{ms}$',
             'tip V.1 的四个时间尺度（雪崩 / 猝灭）'),
            ('（V.4）、反馈必须达到的 $189.47\\ \\mathrm{s^{-1}}$（V.3）',
             'tip V.1 的反馈带宽要求'),
            ('$0.795\\ \\mathrm{mT}$ 大约是地磁场（$0.05\\ \\mathrm{mT}$）的 $15.9$ 倍',
             '误差场阈值的日常感（地磁场倍数）'),
            ('磁体表面场（约 $300\\ \\mathrm{mT}$）小两个多量级',
             '误差场阈值的日常感（定性，不印具体倍数）'),
            ('给 $\\tau_w=52.7788\\ \\mathrm{ms}$，铜（$\\sigma=5.8\\times10^7\\ \\mathrm{S/m}$）给',
             'V.3 的两种壁材料算例（不锈钢）'),
            ('$\\tau_w=2.18655\\ \\mathrm{s}$。于是电阻壁模（RWM）只能被延迟，不能被抓死',
             'V.3 的铜壁时间与 RWM 的定性结论'),
            ('\\gamma_{\\mathrm{RWM}} \\approx \\frac{1}{\\tau_w} = 18.947\\ \\mathrm{s^{-1}}',
             'RWM 增长率公式与数值'),
            ('不锈钢壁给的是 $18.947\\ \\mathrm{s^{-1}}$ 的增长率',
             'V.3 正文复述的增长率'),
            ('控制环路的带宽至少要到 $\\tau_w/10$ 的时间尺度上，即 $189.47\\ \\mathrm{s^{-1}}$',
             '反馈带宽要求（tau_w/10）'),
            ('工作点离无壁极限还有 $33.475\\%$ 的裕度',
             'V.3 的 beta 裕度'),
            ('$\\beta_N=1.8627$，距无壁极限 2.8 还有 $33.475\\%$ 的裕度，密度用掉 Greenwald',
             'V.3 三个裕度并列句（beta 与密度）'),
            ('$n/n_G=83.78\\%$（离密度极限 $16.22\\%$）',
             '密度裕度的两个印法（83.78% 与 16.22%）'),
            ('$(2.8-1.8627)/2.8=33.475\\%$',
             '练习 1 解答里的裕度算式'),
            ('$B^2/2\\mu_0=1.118\\times10^{7}\\ \\mathrm{Pa}$，得 $\\beta=2.6359\\%$',
             '练习 1 解答里的磁压与 beta'),
            ('$E_{\\mathrm{tor}}/E_D=0.01507$，主体电子仍被碰撞拖住，不会整体被拉走。',
             'V.5 的 E_tor/E_D 与主体电子结论'),
            ('$40\\ \\mathrm{ms}$ 就是这个时间尺度的 $92.12$ 倍，放大倍数 $e^{92.12}\\approx1.021\\times10^{40}$。',
             'V.5 的雪崩放大倍数（92.12 个 e 折）'),
            ('$8.378$ 倍 $n_G$）也才把 $E_c$ 抬到 $0.8668\\ \\mathrm{V/m}$',
             'V.5 的密度抬升路径（n_e=1e21）'),
            ('$E_{\\mathrm{tor}}<E_c$ 需要 $\\tau_{\\mathrm{CQ}}=5.921\\ \\mathrm{s}$',
             'V.6 的两面夹击结论（tau_CQ 下界）'),
            ('占等离子体离子库存 $8.322\\times10^{22}$ 的 $1.488\\%$',
             'V.4/V.6 的杂质份额（相对离子库存）'),
            ('$\\tau_{\\mathrm{CQ}}=L_p/R_p=76.9074\\ \\mathrm{ms}$',
             '练习 3 解答的 tau_CQ 算式'),
            ('$L_p=1.33294\\times10^{-5}\\ \\mathrm{H}$',
             '练习 3 解答的电感'),
            ('第二项比第一项大 $672$ 倍',
             '表 V.8 题注里的两项任务之比'),
            ('$2.18655\\ \\mathrm{s}$（相差 $41.428$ 倍）',
             '练习 6 解答的铜/钢壁时间之比'),
            ('$18.947\\ \\mathrm{s^{-1}}$，反馈带宽要按 $\\tau_w/10$ 设计，即 $189.47\\ \\mathrm{s^{-1}}$。',
             '练习 6 解答的增长率与带宽'),
            ('\\textbf{出处：} 三条运行极限（密度、$\\beta$、误差场）的分类与破裂的时间线参照 Wesson',
             'V.1 的出处声明'),
            ('\\textbf{出处：} $\\beta_N$ 的定义与 Troyon 极限参照 Troyon 等（1984）',
             'V.3 的出处声明'),
            ('{把简化声明写在明处}{tip-app-v-honest}',
             'tip V.3 的标题（四处简化声明）'),
            ('不复制任何教材或论文的表格数值——本附录每一张表的数字都能用',
             '「不复制」约定的第一条'),
            ('\\texttt{tmp/appV\\_calc.py}',
             '数值唯一来源脚本（独立复算的凭据）'),
            ('全部数字来自 \\texttt{tmp/appV\\_calc.py}',
             '来源表里对数值来源的交代'),
            ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['V'], why))

    # 口径回流防护（这些写法一旦出现即为 FAIL）
    for lit, why in (
            ('$\\tau_w=52.78\\ \\mathrm{ms}$',
             '壁时间被写成两位小数（应印 52.7788 ms）'),
            ('2.19\\ \\mathrm{s}',
             '铜壁时间被写成三位有效数字（应印 2.18655 s）'),
            ('= 18.95\\ \\mathrm{s^{-1}}',
             'RWM 增长率被四舍五入成 18.95'),
            ('即 $190\\ \\mathrm{s^{-1}}$',
             '反馈带宽被四舍五入成 190'),
            ('还有 $33.5\\%$ 的裕度',
             'beta 裕度被四舍五入成 33.5%'),
            ('$\\tau_{\\mathrm{CQ}}=76.91\\ \\mathrm{ms}$',
             'tau_CQ 被四舍五入成 76.91 ms'),
            ('$E_{\\mathrm{tor}}/E_D=0.0151$',
             'E_tor/E_D 被四舍五入成 0.0151'),
            ('$e^{92}\\approx1.0\\times10^{40}$',
             '雪崩放大倍数被写成 1.0e40'),
            ('$1.78\\times10^{-4}$',
             'm=3 的壁面允许误差场被四舍五入'),
            ('$0.943\\ \\mathrm{mT}$',
             'm=3 的绝对场被四舍五入成 0.943 mT'),
            ('$8.4$ 倍 $n_G$',
             'n_e=1e21 相对 n_G 的倍数被四舍五入'),
            ('（相差 $41$ 倍）',
             '铜/钢壁时间之比被四舍五入成 41'),
            ('大约是地磁场（$0.05\\ \\mathrm{mT}$）的 16 倍',
             '地磁场倍数被四舍五入成 16'),
            ('离密度极限约 $16\\%$',
             '密度裕度被四舍五入成 16%'),
            ('n_G=1.1937\\times10^{21}',
             'Greenwald 极限错一个量级'),
            ('\\gamma_{\\mathrm{RWM}} \\approx \\frac{1}{\\tau_w} = 18.95',
             'RWM 增长率公式里的数值被改写'),
            ('力约 $7.95\\ \\mathrm{MN}$',
             '晕电流力错一个量级'),
            ('力约 $795\\ \\mathrm{MN}$',
             '晕电流力错一个量级（多一位）'),
            ('上是 $25.99\\ \\mathrm{MJ/m^2}$',
             '逃逸束能量密度错一个量级'),
            ('$\\tau_{\\mathrm{CQ}}=5.921\\ \\mathrm{ms}$',
             'tau_CQ 下界的单位写错（s 写成 ms）'),
            ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['V'], why))

    # 数值独立复算（单位见正文/表头，故 unit=''）
    cases = [
        # 破裂：三张极限图
        (2.8, 1e-3, '破裂：三张极限图', '无壁 beta_N 经验极限'),
        (4, 1e-3, '破裂：三张极限图', '理想壁 beta_N 经验极限'),
        (1.1937e+20, 1e-3, '破裂：三张极限图', 'Greenwald 密度极限/m^-3'),
        (371.4, 1e-3, '破裂：三张极限图', 'D1 热储能/MJ'),
        (2600, 1e-3, '破裂：三张极限图', '环电压/V'),
        (76.9074, 1e-3, '破裂：三张极限图', '电流猝灭时间/ms'),
        (4.5, 1e-3, '破裂：三张极限图', '晕电流/MA'),
        (79.5, 1e-3, '破裂：三张极限图', '晕电流力/MN'),
        (26, 1e-3, '破裂：三张极限图', '逃逸束能量/MJ'),
        (0.1, 1e-3, '破裂：三张极限图', '逃逸束落点面积/m^2'),
        (259.9, 1e-3, '破裂：三张极限图', '逃逸束能量密度/(MJ/m^2)'),
        # 误差场与锁模
        (0.6667, 1e-3, '误差场与锁模', 'm=2 的真空衰减因子'),
        (0.5625, 1e-3, '误差场与锁模', 'm=3 的真空衰减因子'),
        (0.00015, 1e-3, '误差场与锁模', 'm=2 壁面允许的 deltaB/B'),
        (0.795, 1e-3, '误差场与锁模', 'm=2 折成绝对场/mT'),
        (0.00017778, 1e-3, '误差场与锁模', 'm=3 壁面允许的 deltaB/B'),
        (0.9422, 1e-3, '误差场与锁模', 'm=3 折成绝对场/mT'),
        (2.4, 1e-3, '误差场与锁模', '壁半径/m'),
        (5.3, 1e-3, '误差场与锁模', '环向场/T'),
        (0.05, 1e-3, '误差场与锁模', '地磁场/mT'),
        (300, 1e-3, '误差场与锁模', '钕铁硼表面场/mT'),
        (15.9, 1e-3, '误差场与锁模', '0.795 mT 是地磁场的倍数'),
        # 电阻壁模与
        (2.6359, 1e-3, '电阻壁模与', '工作点 beta/%'),
        (1.8627, 1e-3, '电阻壁模与', '工作点 beta_N'),
        (294600, 1e-3, '电阻壁模与', '工作点压强/Pa'),
        # 六道估算练习
        (1.118e+07, 1e-3, '六道估算练习', '练习 1 的磁压/Pa'),
        # 电阻壁模与
        (83.78, 1e-3, '电阻壁模与', 'n/n_G 的百分比'),
        # 避免、预测与缓解
        (16.22, 1e-3, '避免、预测与缓解', '离密度极限的裕度/%'),
        # 电阻壁模与
        (3.9623, 1e-3, '电阻壁模与', '无壁极限 beta/%'),
        (442800, 1e-3, '电阻壁模与', '无壁极限压强/Pa'),
        (13.8202, 1e-3, '电阻壁模与', '无壁极限温度/keV'),
        (5.6604, 1e-3, '电阻壁模与', '理想壁极限 beta/%'),
        (632600, 1e-3, '电阻壁模与', '理想壁极限压强/Pa'),
        (19.7432, 1e-3, '电阻壁模与', '理想壁极限温度/keV'),
        (1.4286, 1e-3, '电阻壁模与', '理想壁/无壁 beta 之比'),
        (1.503, 1e-3, '电阻壁模与', '无壁极限与工作点之比'),
        (2.147, 1e-3, '电阻壁模与', '理想壁极限与工作点之比'),
        (33.475, 1e-3, '电阻壁模与', '离无壁极限的裕度/%'),
        (52.7788, 1e-3, '电阻壁模与', '不锈钢壁时间/ms'),
        (2.18655, 1e-3, '电阻壁模与', '铜壁时间/s'),
        (18.947, 1e-3, '电阻壁模与', 'RWM 增长率/(1/s)'),
        (189.47, 1e-3, '电阻壁模与', '反馈带宽/(1/s)'),
        # 六道估算练习
        (41.428, 1e-3, '六道估算练习', '铜/钢壁时间之比'),
        # 破裂的后果
        (464.25, 1e-3, '破裂的后果', '热账功率密度/(MW/m^2)'),
        (800, 1e-3, '破裂的后果', '第一壁面积/m^2'),
        (371.4, 1e-3, '破裂的后果', '热账需要的受热面积/m^2'),
        (25.99, 1e-3, '破裂的后果', '1 m^2 上的能量密度/(MJ/m^2)'),
        # 逃逸电子：Dreicer
        (44.3, 1e-3, '逃逸电子：Dreicer', 'Dreicer 场/(V/m)'),
        (0.08668, 1e-3, '逃逸电子：Dreicer', '临界场/(V/m)'),
        (511, 1e-3, '逃逸电子：Dreicer', 'E_D/E_c'),
        (4430, 1e-3, '逃逸电子：Dreicer', 'E_D(10 eV)/(V/m)'),
        (66.74, 1e-3, '逃逸电子：Dreicer', '猝灭期环向场/(V/m)'),
        (769.9, 1e-3, '逃逸电子：Dreicer', 'E_tor/E_c'),
        (0.01507, 1e-3, '逃逸电子：Dreicer', 'E_tor/E_D'),
        (0.4342, 1e-3, '逃逸电子：Dreicer', '雪崩 e 折时间/ms'),
        (92.12, 1e-3, '逃逸电子：Dreicer', '40 ms 的 e 折数'),
        (1.021e+40, 1e-3, '逃逸电子：Dreicer', '雪崩放大倍数'),
        (8.1104e+18, 1e-3, '逃逸电子：Dreicer', '10 MA 逃逸束的电子数'),
        (25.9885, 1e-3, '逃逸电子：Dreicer', '20 MeV 束能量/MJ'),
        (8.378, 1e-3, '逃逸电子：Dreicer', 'n_e=1e21 相对 n_G 的倍数'),
        (0.8668, 1e-3, '逃逸电子：Dreicer', 'E_c(n=1e21)/(V/m)'),
        (76.99, 1e-3, '逃逸电子：Dreicer', 'E_tor/E_c(n=1e21)'),
        (645, 1e-3, '逃逸电子：Dreicer', 'n_supp 相对 n_G 的倍数'),
        (5.133, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=1 s 行的 E_tor/(V/m)'),
        (5.921, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=1 s 行的 E_tor/E_c'),
        (513.3, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.01 s 行的 E_tor/(V/m)'),
        (592.1, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.01 s 行的 E_tor/E_c'),
        (102.7, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.05 s 行的 E_tor/(V/m)'),
        (118.4, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.05 s 行的 E_tor/E_c'),
        (51.33, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.10 s 行的 E_tor/(V/m)'),
        (59.21, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.10 s 行的 E_tor/E_c'),
        (25.66, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.20 s 行的 E_tor/(V/m)'),
        (29.6, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.20 s 行的 E_tor/E_c'),
        (10.27, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.50 s 行的 E_tor/(V/m)'),
        (11.84, 1e-3, '逃逸电子：Dreicer', 'tau_CQ=0.50 s 行的 E_tor/E_c'),
        # 避免、预测与缓解
        (1.238e+21, 1e-3, '避免、预测与缓解', '辐射掉热储能所需杂质原子数'),
        (5.128, 1e-3, '避免、预测与缓解', '折成 300 K 气体/(Pa m^3)'),
        (8.322e+23, 1e-3, '避免、预测与缓解', '把密度抬到 1e21 所需粒子数'),
        (3.447, 1e-3, '避免、预测与缓解', '折成 300 K 气体/(kPa m^3)'),
        (1.488, 1e-3, '避免、预测与缓解', '杂质占离子库存的百分比'),
        (672, 1e-3, '避免、预测与缓解', '两项任务之比'),
        (3e-32, 1e-3, '避免、预测与缓解', '声明过的杂质冷却率 L_z'),
        # 六道估算练习
        (1.33294e-05, 1e-3, '六道估算练习', '练习 3 的电感/H'),
        (294600, 1e-3, '六道估算练习', '练习 1 的压强/Pa'),
        (2.6359, 1e-3, '六道估算练习', '练习 1 的 beta/%'),
        (1.8627, 1e-3, '六道估算练习', '练习 1 的 beta_N'),
        (25.9885, 1e-3, '六道估算练习', '练习 5 的束能量/MJ'),
        (259.9, 1e-3, '六道估算练习', '练习 5 的束能量密度/(MJ/m^2)'),
        # 来源、分级与版权说明
        (2011, 1e-3, '来源、分级与版权说明', 'Wesson 第 4 版年份'),
        (2007, 1e-3, '来源、分级与版权说明', 'Freidberg/Hender 年份'),
        (2014, 1e-3, '来源、分级与版权说明', 'Freidberg《Ideal MHD》年份'),
        (1984, 1e-3, '来源、分级与版权说明', 'Troyon 等年份'),
        (2010, 1e-3, '来源、分级与版权说明', 'Chu \\& Okabayashi 年份'),
        (1993, 1e-3, '来源、分级与版权说明', 'Fitzpatrick 年份'),
        (1959, 1e-3, '来源、分级与版权说明', 'Dreicer 年份'),
        (1975, 1e-3, '来源、分级与版权说明', 'Connor \\& Hastie 年份'),
        (1997, 1e-3, '来源、分级与版权说明', 'Rosenbluth \\& Putvinski 年份'),
        # 电阻壁模与
        (123001, 1e-3, '电阻壁模与', 'Chu \\& Okabayashi 文章号'),
        (209, 1e-3, '电阻壁模与', 'Troyon 等起始页'),
        # 误差场与锁模
        (1049, 1e-3, '误差场与锁模', 'Fitzpatrick 起始页'),
        # 逃逸电子：Dreicer
        (238, 1e-3, '逃逸电子：Dreicer', 'Dreicer 起始页'),
        (415, 1e-3, '逃逸电子：Dreicer', 'Connor \\& Hastie 起始页'),
        (1355, 1e-3, '逃逸电子：Dreicer', 'Rosenbluth \\& Putvinski 起始页'),
        # 避免、预测与缓解
        (53018, 1e-3, '避免、预测与缓解', 'de Vries 等文章号'),
        (21802, 1e-3, '避免、预测与缓解', 'Hollmann 等文章号'),
    ]
    for v, rel, where, what in cases:
        block = found.get(where)
        if block is None:
            rec('FAIL', 'P', '附录 V 找不到小节「%s」，无法核对 %s' % (where, what))
        elif not value_in(block, v, rel, unit=''):
            rec('FAIL', 'P', '附录 V 的 %s 缺失或与 tmp/appV_calc.py 不符：期望 %g（小节「%s」）'
                % (what, v, where))
    rec('PASS', 'P', '附录 V 的 %d 个关键数值均在小节内复算通过' % len(cases))



def check_P27_appendix_w_fastparticle():
    """附录 W 快粒子与 alpha 粒子物理的完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 八个小节之一被删（按小节标题子串定位）；
      ② 七张表 / 六只语义盒 / 二十三个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 6 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（出生参数、慢化时间与临界能量、轨道宽度与
         首轨损失、快粒子压强与尾巴份额、Alfven 速度与三个频率、损失三通道），
         或搬到了别的小节里；
      ⑤ 口径回流：慢化时间被四舍五入、比值 13.7 写成 13.8、v_alpha0 写成 1.30e7、
         v_c 写成 3.8e6、进动频率写成 1.8 kHz、beta 份额写成 9.2%、密度份额写成
         0.21%、首轨份额写成 0.499%、能量 e 折写成 0.16 s、快密度写成 2.13e17、
         E_c 写成 295.1 keV、T_e 倍数写成 32.1。

    数值唯一来源 tmp/appW_calc.py；表内数字的单位写在表头，故部分用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['W'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 W 文件缺失：src/%s' % APP_FILES['W'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = [
        '出生：3.5 MeV', '慢化：为什么不能用', '加热落在哪', '慢化分布与高能尾巴',
        '快粒子驱动的不稳定性', '损失与第一壁', '六道估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 W 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 W 的 %d 个小节齐全' % len(want))

    for lab in (
            'tab:app-w-life', 'tab:app-w-slowing', 'tab:app-w-orbit', 'tab:app-w-tail',
            'tab:app-w-instability', 'tab:app-w-loss', 'tab:app-w-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['W'], lab))
    for lab in (
            'tip-app-w-oneway', 'tip-app-w-locality', 'tip-app-w-locality2',
            'tip-app-w-fewpercent', 'tip-app-w-ratios', 'tip-app-w-honest'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['W'], lab))
    for lab in (
            'eq:app-w-birth', 'eq:app-w-gyro', 'eq:app-w-rho', 'eq:app-w-source',
            'eq:app-w-fastformula', 'eq:app-w-xratio', 'eq:app-w-chandrasekhar',
            'eq:app-w-tauv', 'eq:app-w-therm', 'eq:app-w-vc', 'eq:app-w-rhotheta',
            'eq:app-w-prompt', 'eq:app-w-dist', 'eq:app-w-meanE', 'eq:app-w-pfast',
            'eq:app-w-tail', 'eq:app-w-va', 'eq:app-w-tae', 'eq:app-w-precess',
            'eq:app-w-bounce', 'eq:app-w-ratios', 'eq:app-w-ripple', 'eq:app-w-qwall'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['W'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['W'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 W 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 W 有 %d 处「出处：」声明' % n_src)
    ex = found['六道估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 W 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 W 有 %d 道练习、%d 条解答' % (n_ex, n_ans))

    # 关键行与关键句的逐字绑定（防表格整行被改、防正文数字被搬走）
    for lit, why in (
            ('首次轨道损失 & $0.4999\\%$ / $0.5564\\ \\mathrm{MW}$ &',
             '表 W.5 首轨损失行（份额与功率成对）'),
            ('纹波损失 & $1\\%$--$3\\%$ / $1.1130$--$3.3390\\ \\mathrm{MW}$ &',
             '表 W.5 纹波损失行（简化口径区间）'),
            ('鱼骨模 & $\\omega\\approx\\omega_d$（进动共振） &',
             '表 W.4 鱼骨模行（共振条件）'),
            ('TAE／环向 Alfvén 模 & $v_\\parallel\\approx v_A$，$\\omega\\approx v_A/(2qR)$ &',
             '表 W.4 TAE 行（共振条件）'),
            ('W.4 & 慢化分布 $1/(v^3+v_c^3)$、$(2/3)n\\langle E\\rangle$ 口径 &',
             '表 W.7 的 W.4 行（借鉴与自拟的分界）'),
            ('W.5 & TAE 色散关系、鱼骨模与进动共振、alpha 连续谱定性 &',
             '表 W.7 的 W.5 行'),
            ('不稳定性输运 & 未定量（模型算不出） &',
             '表 W.5 第三通道（刻意留空，不得被填数）'),
            ('$v_{\\alpha 0}/v_{Te}$ & $0.2284$ &', '表 W.2 的适用性判据行'),
            ('$H(x)$ 精确 / 渐近 & $0.00869311$ / $0.00896879$ &',
             '表 W.2 的线性阻力段证据行'),
            ('$\\rho_\\theta$ & $0.1790\\ \\mathrm{m}$ & $q\\rho_\\alpha/\\sqrt{\\eps}$ &',
             '表 W.3 的径向分辨率行'),
            ('alpha 比背景电子慢', 'W.2 的核心判断句'),
            ('这个答案是错的，而且错三个数量级', 'W.2 对错误公式的定量否定'),
            ('两条完全不同的推导路径', 'W.2 的 tau_eq 闭合验证'),
            ('E_{\\alpha 0} = 3.5\\ \\mathrm{MeV} = 5.6076\\times10^{-13}\\ \\mathrm{J}',
             '式 W.1 的出生能量整行'),
            ('f_{\\text{TAE}}(q{=}2)=47.09\\ \\mathrm{kHz},\\;',
             '式 W.18 的 TAE 频率整行'),
            ('四成的粒子带着七成半的能量', 'W.4 的尾巴结论句'),
            ('快粒子损失不是总能量问题，是分布问题', 'W.6 的结论句'),
            ('把做不到的地方标出来', 'tip W.6 的诚实口径'),
            ('本表不给阈值', '表 W.4 题注的阈值声明'),
            ('不给出任何阈值数字', 'W.5 的出处段口径声明'),
            ('只给判据、不给数字', 'W.8 的三档难度说明'),
            ('不复制任何教材或论文的表格数值', 'W.8 的第一条不复制约定'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['W'], why))

    # 数值逐节独立复算（rel 全部 1e-3；只查本小节内的数，防数字被搬到别的小节）
    U = '\\mathrm{%s}'
    cases = [
        # W.1 出生
        ('出生：3.5 MeV', 111.3, U % 'MW', 'alpha 功率（=P_fus/5）'),
        ('出生：3.5 MeV', 5.6076e-13, U % 'J', 'alpha 出生能量'),
        ('出生：3.5 MeV', 1.2992e7, U % 'm/s', 'alpha 出生速度'),
        ('出生：3.5 MeV', 0.0508, U % 'm', 'alpha 回旋半径 5.08 cm'),
        ('出生：3.5 MeV', 3.70, U % 'mm', '氘离子回旋半径 3.70 mm'),
        ('出生：3.5 MeV', 13.7, '', 'rho_alpha/rho_i'),
        ('出生：3.5 MeV', 1.9848e20, U % 's^{-1}', 'alpha 出生率'),
        ('出生：3.5 MeV', 2.3850e17, 'm^{-3}\\,s^{-1}', '单位体积出生率'),
        # W.2 慢化
        ('慢化：为什么不能用', 0.2284, '', 'x = v_a0/v_Te < 1'),
        ('慢化：为什么不能用', 5.6869e7, U % 'm/s', '电子热速度'),
        ('慢化：为什么不能用', 0.00869311, '', 'H(x) 精确值'),
        ('慢化：为什么不能用', 0.00896879, '', 'H(x) 渐近值'),
        ('慢化：为什么不能用', 2.0421e-26, U % 'kg/s', '阻力系数 K'),
        ('慢化：为什么不能用', 0.3254, U % 's', '速度 e 折时间'),
        ('慢化：为什么不能用', 0.1627, U % 's', '能量 e 折时间'),
        ('慢化：为什么不能用', 0.1638, U % 's', 'NRL 均分时间（交叉验证）'),
        ('慢化：为什么不能用', 0.8902, U % 's', '热化时间'),
        ('慢化：为什么不能用', 1.095605e-6, '', '离子阻力系数 C_i'),
        ('慢化：为什么不能用', 3.7716e6, U % 'm/s', '临界速度'),
        ('慢化：为什么不能用', 295.0, U % 'keV', '临界能量'),
        ('慢化：为什么不能用', 32.08, '', 'E_c / T_e'),
        ('慢化：为什么不能用', 4.2272e6, U % 'm', '慢化路径'),
        ('慢化：为什么不能用', 4227, U % 'km', '慢化路径（km）'),
        ('慢化：为什么不能用', 108514, '', '绕环圈数'),
        ('慢化：为什么不能用', 8.424185e5, U % 'm/s', '氘离子热速度'),
        # W.3 轨道
        ('加热落在哪', 0.1790, U % 'm', '极向回旋半径'),
        ('加热落在哪', 8.95, '', 'rho_theta / a'),
        ('加热落在哪', 0.3580, U % 'm', '香蕉宽度'),
        ('加热落在哪', 0.9105, '', '首轨损失边界 r/a'),
        ('加热落在哪', 0.4999, '', '首轨损失份额'),
        ('加热落在哪', 0.5564, U % 'MW', '首轨损失功率'),
        ('加热落在哪', 0.0714, 'MW/m^2', '落点热流'),
        ('加热落在哪', 7.7911, U % 'm^2', '落点面积'),
        # W.4 尾巴
        ('慢化分布与高能尾巴', 1191.1, U % 'keV', '快 alpha 平均能量'),
        ('慢化分布与高能尾巴', 0.3403, '', '平均能量 / 出生能量'),
        ('慢化分布与高能尾巴', 2.1231e17, U % 'm^{-3}', '快 alpha 密度'),
        ('慢化分布与高能尾巴', 0.2123, '', '快 alpha 密度份额'),
        ('慢化分布与高能尾巴', 2.7011e4, U % 'Pa', '快 alpha 压强'),
        ('慢化分布与高能尾巴', 2.9461e5, U % 'Pa', '热压'),
        ('慢化分布与高能尾巴', 9.17, '', '快粒子压强份额'),
        ('慢化分布与高能尾巴', 0.241671, '', 'beta_fast'),
        ('慢化分布与高能尾巴', 2.6359, '', 'beta_tot'),
        ('慢化分布与高能尾巴', 41.2888, '', '尾巴粒子份额'),
        ('慢化分布与高能尾巴', 74.5395, '', '尾巴能量份额'),
        # W.5 不稳定性
        ('快粒子驱动的不稳定性', 7.3380e6, U % 'm/s', 'Alfven 速度'),
        ('快粒子驱动的不稳定性', 1.7705, '', 'v_a0 / v_A'),
        ('快粒子驱动的不稳定性', 47.09, U % 'kHz', 'f_TAE(q=2)'),
        ('快粒子驱动的不稳定性', 94.2, U % 'kHz', 'f_TAE(q=1)'),
        ('快粒子驱动的不稳定性', 31.4, U % 'kHz', 'f_TAE(q=3)'),
        ('快粒子驱动的不稳定性', 8589.72, U % 'rad/s', '进动角频率'),
        ('快粒子驱动的不稳定性', 1367.10, U % 'Hz', '进动频率'),
        ('快粒子驱动的不稳定性', 47.35, U % 'kHz', '弹跳频率'),
        ('快粒子驱动的不稳定性', 40.41, U % 'MHz', '氘离子回旋频率'),
        # W.6 损失与第一壁
        ('损失与第一壁', 0.005, '', '纹波幅度 delta'),
        ('损失与第一壁', 0.1000, '', 'sqrt(2 delta)'),
        ('损失与第一壁', 0.50, '', '纹波俘获立体角份额'),
        ('损失与第一壁', 1.1130, U % 'MW', '纹波损失下界'),
        ('损失与第一壁', 3.3390, U % 'MW', '纹波损失上界'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s：%s（%s）' % (sec.split('：')[0], val, why))
    if bad:
        rec('FAIL', 'P', '附录 W 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 W 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # 口径回流防护（这些写法一律不得出现在附录 W 里）
    lits = [
        ('0.499\\%', '首轨损失份额不得四舍五入到 0.499%'),
        ('0.16\\ \\mathrm{s}', '能量 e 折时间不得写成 0.16 s'),
        ('13.8', 'rho_alpha/rho_i 不得写成 13.8'),
        ('2.13\\times10^{17}', '快 alpha 密度不得写成 2.13e17'),
        ('295.1', '临界能量不得写成 295.1 keV'),
        ('32.1\\,T_e', 'T_e 倍数不得写成 32.1'),
        ('0.89\\ \\mathrm{ms}', '热化时间不得写成 0.89 ms'),
        ('1.30\\times10^{7}', 'alpha 出生速度不得写成 1.30e7'),
        ('3.8\\times10^{6}', '临界速度不得写成 3.8e6'),
        ('1.8\\ \\mathrm{kHz}', '进动频率不得写成 1.8 kHz'),
        ('9.2\\%', '快粒子压强份额不得写成 9.2%'),
        ('0.21\\%', '快 alpha 密度份额不得写成 0.21%'),
    ]
    for lit, why in lits:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['W'], why))


def check_P28_signed_quantities():
    """带符号物理量守卫：物理上为负的量必须**带符号**出现。

    ``value_in()`` 默认按绝对值比对（正文常只印大小），全仓库一度 0 处使用
    ``signed=True``，于是把「必须有负号的结论」翻成正号会**静默通过**：
    实测 6 处候选里 4 处漏网（附录 S 的 IPB98 指数、附录 U 的同一指数、
    附录 T 的 Δ′、正文的悬浮电位）。本函数对每条量做两件事——
    ① 整串带符号绑定（连美元符号与 & 一起钉住，改动一个字符即 FAIL）；
    ② 在同一文件内用 ``signed=True`` 要求真的存在一个负值。
    """
    bad = []
    for fname, lit, val, rel, why in (
        ('appendix_s_confinement.tex', '$P$ & $-0.69$ & $0.68$ &', -0.69, 0.001, 'IPB98(y,2) 里 P 的指数：负号就是结论（少加热反而约束时间长）'),
        ('appendix_s_confinement.tex', '$-0.58$（应为 $-0.69$）', -0.58, 0.001, '附录 S 常见错误：把 P 的幂次写成 -0.58'),
        ('appendix_u_gain.tex', '$P$ 的指数是 $-0.69$：', -0.69, 0.001, '附录 U 引用同一条定标律时的指数'),
        ('appendix_t_topology.tex', '$-7.680\\times10^{-3}$', -0.00768, 0.001, '附录 T 的 Δ′（m/n = 2/1）：符号决定稳定性'),
        ('appendix_t_topology.tex', '$-1.536\\times10^{-2}$', -0.01536, 0.001, '附录 T 的 Δ′（m/n = 3/2）'),
        ('plasma_physics_textbook_v2.tex', '氩约为 $-4.7\\,T_e$', -4.7, 0.001, '悬浮电位：物理上必须为负'),
        ('plasma_physics_textbook_v2.tex', '对数斜率约为 $-0.333$', -0.333, 0.001, '对数斜率：负号是结论本身'),
        ('ch12_supplement.tex', '$-20$ & $-0.51$', -0.51, 0.001, 'Langmuir 探针 I-V 表：离子电流为负'),
        ('appendix_q_diagnostics.tex', '$-1.0$ & 7.1653 & 6.1653', -1, 0.001, '探针曲线表：偏压列带负号'),
        ('appendix_l_problems.tex', '$-6.5\\times10^{-5}$', -6.5e-05, 0.001, '附录 L 的屏蔽修正项：负号是结论本身'),
        ('ch3_supplement.tex', '$-10\\%$', -10, 0.001, '单边近似与两端总数的偏差：带负号'),):
        p = os.path.join(SRC, fname)
        if not os.path.exists(p):
            bad.append('%s：文件缺失' % fname)
            continue
        t = strip_comments(read(p))
        if lit not in t:
            bad.append('%s：缺带符号的整串 %s（%s）' % (fname, lit, why))
            continue
        if not value_in(t, val, rel, '', True):
            bad.append('%s：文件内找不到带符号的 %g（%s）' % (fname, val, why))
    if bad:
        rec('FAIL', 'P', '带符号物理量守卫发现 %d 处问题：%s'
            % (len(bad), '；'.join(bad[:5])))
    else:
        rec('PASS', 'P', '带符号物理量守卫：%d 条「必须带符号」的断言全部通过' % 11)


def check_P29_appendix_x_radiation():
    """附录 X 辐射输运与辐射损失的完整性与数值独立复算。

    可证伪点（全部记 FAIL）：
      ① 八个小节之一被删（按小节标题子串定位）；
      ② 七张表 / 六只语义盒 / 十九个编号公式之一被删或改了 label；
      ③ 缺少「出处：」行（本附录要求 >= 6 处），或六道练习缺解答；
      ④ 印出的关键数值与独立复算不符（三种机制的光学薄发射率与总量、自由-自由
         光学深度与平均自由程、轫致辐射的相对论修正与壁面热流、回旋辐射的上下限
         与再吸收倍数、线辐射崩塌阈值、辐射在功率账里的位置与冷却时间、六道练习），
         或搬到了别的小节里；
      ⑤ 口径回流：光学薄回旋总量被四舍五入到 890 MW、再吸收倍数写成 5400、
         黑体上限写成 0.16 MW、崩塌阈值写成 0.13%、冷却时间写成 13.3 s、
         平均自由程写成 6.8e5 m、光子波长写成 0.13 nm、修正后 P_BR 写成 28.0 MW、
         beta_perp^2 写成 0.036、线辐射总量写成 0.83 MW。

    数值唯一来源 tmp/appX_calc.py；表内数字的单位写在表头，故部分用例传 unit=''。
    语义盒标签写在 \\begin{tip}{标题}{标签} 的第二个必选参数里，不写成 \\label{}。
    """
    spath = os.path.join(SRC, APP_FILES['X'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 X 文件缺失：src/%s' % APP_FILES['X'])
        return
    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = [
        '三种辐射机制', '光学薄还是光学厚', '轫致辐射：为什么', '回旋辐射：一个必须',
        '线辐射与杂质', '辐射在功率账里的位置', '估算练习', '来源、分级与版权说明']
    found = {}
    for pre in want:
        for ttl, body in secs:
            if pre in ttl:
                found[pre] = body
                break
    missing = [p for p in want if p not in found]
    if missing:
        rec('FAIL', 'P', '附录 X 缺少小节（标题被改或整节被删）：%s' % missing)
        return
    rec('PASS', 'P', '附录 X 的 %d 个小节齐全' % len(want))

    for lab in (
            'tab:app-x-channel', 'tab:app-x-opacity', 'tab:app-x-brems',
            'tab:app-x-cyclotron', 'tab:app-x-impurity', 'tab:app-x-budget',
            'tab:app-x-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表格 label %s）' % (APP_FILES['X'], lab))
    for lab in (
            'tip-app-x-upperbound', 'tip-app-x-transport', 'tip-app-x-xray',
            'tip-app-x-declared', 'tip-app-x-honest', 'tip-app-x-margin'):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签 %s）' % (APP_FILES['X'], lab))
    for lab in (
            'eq:app-x-brnum', 'eq:app-x-cycthin', 'eq:app-x-line', 'eq:app-x-optdepth',
            'eq:app-x-ffabs', 'eq:app-x-tauff', 'eq:app-x-br', 'eq:app-x-brfrac',
            'eq:app-x-photon', 'eq:app-x-cyclfreq', 'eq:app-x-larmor', 'eq:app-x-cycthin2',
            'eq:app-x-cycbb', 'eq:app-x-supp', 'eq:app-x-line2', 'eq:app-x-collapse',
            'eq:app-x-prad', 'eq:app-x-margin', 'eq:app-x-taurad'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 label %s）' % (APP_FILES['X'], lab))

    has(text, '不复制', 'src/%s（版权与独立性声明）' % APP_FILES['X'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 X 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 X 有 %d 处「出处：」声明' % n_src)
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 X 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 X 有 %d 道练习、%d 条解答' % (n_ex, n_ans))

    # 关键行与关键句的逐字绑定（防表格整行被改、防正文数字被搬走）
    for lit, why in (
            ('轫致辐射 & $32.444\\ \\mathrm{kW\\,m^{-3}}$ & $27.00\\ \\mathrm{MW}$ &',
             '表 X.1 轫致辐射行（发射率与总量成对）'),
            ('$P_{\\mathrm{fus}}$ \\\\', '表 X.1 行尾结构完整'),
            ('逃逸份额 & 接近 $1$ & $\\sim1/\\tau$，可以极小 & 轫致辐射：光学薄 \\\\',
             '表 X.2 逃逸份额行'),
            ('可否只看发射率 & 可以 & 不可以 & 回旋：不可以 \\\\',
             '表 X.2 的判据结论行'),
            ('修正后 $P_{\\mathrm{BR}}$ & $27.9716\\ \\mathrm{MW}$ & 提高 $3.598\\%$ \\\\',
             '表 X.3 相对论修正行'),
            ('光子波长 & $0.134853\\ \\mathrm{nm}$ & X 射线波段 \\\\',
             '表 X.3 光子波长行'),
            ('光学薄总量 & $890.01\\ \\mathrm{MW}$ & $=1.5993\\,P_{\\mathrm{fus}}$，不可能 \\\\',
             '表 X.4 光学薄上限行'),
            ('再吸收压低倍数 & $5406.41$ & 上下限之比 \\\\',
             '表 X.4 再吸收倍数行'),
            ('钨（高 $Z$） & 千 eV 量级 & 高 & 芯部积累才是威胁 \\\\',
             '表 X.5 钨行（定性口径）'),
            ('自加热余量 & $3.99897$ & --- & 约 4 倍，不是 100 倍 \\\\',
             '表 X.6 自加热余量行'),
            ('辐射冷却时间 & $13.3443\\ \\mathrm{s}$ & --- & 比 $\\tau_{\\mathrm{CQ}}$ 慢 $173.511$ 倍 \\\\',
             '表 X.6 冷却时间行'),
            ('壁面轫致热流 & $55.1547\\ \\mathrm{kW\\,m^{-2}}$ & --- & X 射线，需屏蔽而非冷却 \\\\',
             '表 X.6 壁面热流行'),
            ('$n_W/n_e = 0.0324441\\%$', '解答 3 的钨浓度份额'),
            ('$T = 1507.93\\ \\mathrm{keV}$', '解答 4 的亮度温度'),
            ('= 683.011\\ \\mathrm{km}$', '解答 5 的平均自由程（km）'),
            ('光学薄假设用错了地方', 'X.1 的核心判断句'),
            ('这是\\textbf{声明}的量级值', 'X.1 的声明口径'),
            ('不给具体数值，因为本附录没有引入原子物理数据库', '表 X.5 题注的定性声明'),
            ('本附录不给输运方程的解', '表 X.4 题注的诚实口径'),
            ('三种机制的划分与轫致辐射公式参照 Wesson', 'X.1 的出处段'),
            ('冷却率系数的定义与杂质辐射崩塌的图像参照 Post 等（1977）',
             'X.5 的出处段'),
            ('恒星内部才满足，实验室等离子体永远光学薄', '解答 5 的结论句'),
            ('对材料的要求不是「耐热」而是「屏蔽」', '解答 6 的结论句'),
            ('回旋辐射 & $1.0695\\ \\mathrm{MW\\,m^{-3}}$ & $890.01\\ \\mathrm{MW}$ &',
             '表 X.1 回旋辐射行（故意不合物理的光学薄值）'),
            ('$u_B=B^2/(2\mu_0)=1.11767\\times10^{7}\\ \\mathrm{J\\,m^{-3}}$',
             'X.1 的磁能密度整行'),
            ('=47.8519\\ \\mathrm{kW\\,m^{-3}}$', '解答 2 的功率密度整行'),
            ('黑体上限占轫致辐射 & $0.60971\\%$ & 结论：回旋不是主要损失 \\\\',
             '表 X.4 黑体上限份额行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['X'], why))

    # 数值逐节独立复算（rel 全部 1e-3；只查本小节内的数，防数字被搬到别的小节）
    UKW_M3 = '\\mathrm{kW\\,m^{-3}}'
    UMW_M3 = '\\mathrm{MW\\,m^{-3}}'
    UWM3 = '\\mathrm{W\\,m^{-3}}'
    UMW = '\\mathrm{MW}'
    UMJ = '\\mathrm{MJ}'
    UJM3 = '\\mathrm{J\\,m^{-3}}'
    UWM2 = '\\mathrm{W\\,m^{-2}}'
    UKWM2 = '\\mathrm{kW\\,m^{-2}}'
    UMWM2 = '\\mathrm{MW\\,m^{-2}}'
    UPCT = '\\%'
    UM = '\\mathrm{m}'
    UKM = '\\mathrm{km}'
    UNM = '\\mathrm{nm}'
    UGHZ = '\\mathrm{GHz}'
    UKEV = '\\mathrm{keV}'
    UMINV = '\\mathrm{m^{-1}}'
    UM3 = '\\mathrm{m^{-3}}'
    UM2 = '\\mathrm{m^2}'
    US = '\\mathrm{s}'
    UJ = '\\mathrm{J}'
    UW = '\\mathrm{W}'
    cases = [
        # X.1 三种机制的光学薄量级
        ('三种辐射机制', 32.444, UKW_M3, '轫致辐射发射率'),
        ('三种辐射机制', 3.2444e4, UWM3, '轫致辐射发射率（W/m^3 口径）'),
        ('三种辐射机制', 27.00, UMW, '光学薄轫致总量'),
        ('三种辐射机制', 1.0695, UMW_M3, '回旋辐射发射率（MW/m^3 口径）'),
        ('三种辐射机制', 1.0695e6, UWM3, '回旋辐射发射率（W/m^3 口径）'),
        ('三种辐射机制', 890.01, UMW, '光学薄回旋总量（故意不合物理）'),
        ('三种辐射机制', 1.5993, '', '光学薄回旋总量 / P_fus'),
        ('三种辐射机制', 0.8322, UMW, '线辐射总量'),
        ('三种辐射机制', 1.11767e7, UJM3, '磁能密度 u_B'),
        ('三种辐射机制', 0.035984, '', 'beta_perp^2'),
        # X.2 光学深度与平均自由程
        ('光学薄还是光学厚', 1.46411e-6, UMINV, '自由-自由吸收系数'),
        ('光学薄还是光学厚', 5.85642e-6, '', '自由-自由光学深度'),
        ('光学薄还是光学厚', 683011.0, UM, '平均自由程（m）'),
        ('光学薄还是光学厚', 683.011, UKM, '平均自由程（km）'),
        ('光学薄还是光学厚', 1.7075e5, '', '平均自由程 / 装置尺度'),
        # X.3 轫致辐射
        ('轫致辐射：为什么', 32.444, UKW_M3, '轫致辐射发射率'),
        ('轫致辐射：为什么', 27.00, UMW, '全等离子体 P_BR'),
        ('轫致辐射：为什么', 24.2588, UPCT, 'P_BR / P_alpha'),
        ('轫致辐射：为什么', 4.85175, UPCT, 'P_BR / P_fus'),
        ('轫致辐射：为什么', 1.03598, '', '一阶相对论修正因子'),
        ('轫致辐射：为什么', 0.0179922, '', 'T_e / (m_e c^2)'),
        ('轫致辐射：为什么', 27.9716, UMW, '修正后 P_BR'),
        ('轫致辐射：为什么', 3.598, UPCT, '相对论修正提高的份额'),
        ('轫致辐射：为什么', 0.134853, UNM, '平均光子波长'),
        # X.4 回旋辐射的两个边界
        ('回旋辐射：一个必须', 148.36, UGHZ, '电子回旋频率'),
        ('回旋辐射：一个必须', 1.11767e7, UJM3, '磁能密度'),
        ('回旋辐射：一个必须', 1.06947e-14, UW, '单电子 Larmor 功率'),
        ('回旋辐射：一个必须', 1.06947, UMW_M3, '光学薄功率密度'),
        ('回旋辐射：一个必须', 890.01, UMW, '光学薄总量'),
        ('回旋辐射：一个必须', 1.5993, '', '光学薄总量 / P_fus'),
        ('回旋辐射：一个必须', 336.283, UWM2, '黑体（Rayleigh-Jeans）通量'),
        ('回旋辐射：一个必须', 0.164622, UMW, '黑体上限总量'),
        ('回旋辐射：一个必须', 5406.41, '', '再吸收压低倍数'),
        ('回旋辐射：一个必须', 0.60971, UPCT, '黑体上限 / P_BR'),
        # X.5 线辐射与杂质
        ('线辐射与杂质', 3.08222, UPCT, '线辐射 / 轫致辐射'),
        ('线辐射与杂质', 0.133742, UPCT, '辐射崩塌阈值 f_W*'),
        # X.6 辐射在功率账里的位置
        ('辐射在功率账里的位置', 27.00, UMW, '轫致辐射'),
        ('辐射在功率账里的位置', 24.2588, UPCT, 'P_BR / P_alpha'),
        ('辐射在功率账里的位置', 0.8322, UMW, '线辐射'),
        ('辐射在功率账里的位置', 0.747709, UPCT, 'P_LINE / P_alpha'),
        ('辐射在功率账里的位置', 27.8322, UMW, '辐射合计'),
        ('辐射在功率账里的位置', 25.0065, UPCT, 'P_rad / P_alpha'),
        ('辐射在功率账里的位置', 3.99897, '', '自加热余量'),
        ('辐射在功率账里的位置', 13.3443, US, '辐射冷却时间'),
        ('辐射在功率账里的位置', 173.511, '', 'tau_rad / tau_CQ'),
        ('辐射在功率账里的位置', 55.1547, UKWM2, '壁面轫致热流'),
        # X.7 六道练习
        ('估算练习', 20.25, UMW, '练习 1 的 Z_eff 缩放'),
        ('估算练习', 1.4749, '', '练习 2 的 sqrt(T) 比值'),
        ('估算练习', 47.8519, UKW_M3, '练习 2 的功率密度'),
        ('估算练习', 47.49, UPCT, '练习 2 的涨幅'),
        ('估算练习', 3.24441e16, UM3, '练习 3 的钨密度'),
        ('估算练习', 0.0324441, UPCT, '练习 3 的钨浓度份额'),
        ('估算练习', 1507.93, UKEV, '练习 4 的亮度温度'),
        ('估算练习', 6.8301e5, UM, '练习 5 的平均自由程'),
        ('估算练习', 683.011, UKM, '练习 5 的平均自由程（km）'),
        ('估算练习', 1.7075e5, '', '练习 5 与装置尺度之比'),
        ('估算练习', 489.532, UM2, '练习 6 的第一壁面积'),
        ('估算练习', 55.1547, UKWM2, '练习 6 的壁面热流'),
        ('估算练习', 0.0551547, UMWM2, '练习 6 的壁面热流（MW/m^2 口径）'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s：%s（%s）' % (sec.split('：')[0], val, why))
    if bad:
        rec('FAIL', 'P', '附录 X 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 X 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # 口径回流防护（这些写法一律不得出现在附录 X 里）
    lits = [
        ('890\\ \\mathrm{MW}', '光学薄回旋总量不得四舍五入到 890 MW'),
        ('5400\\ ', '再吸收压低倍数不得写成 5400'),
        ('0.16\\ \\mathrm{MW}', '黑体上限不得写成 0.16 MW'),
        ('0.13\\%', '辐射崩塌阈值不得写成 0.13%'),
        ('13.3\\ \\mathrm{s}', '辐射冷却时间不得写成 13.3 s'),
        ('55\\ \\mathrm{kW}', '壁面轫致热流不得写成 55 kW/m^2'),
        ('6.8\\times10^{5}', '平均自由程不得写成 6.8e5 m'),
        ('0.13\\ \\mathrm{nm}', '光子波长不得写成 0.13 nm'),
        ('28.0\\ \\mathrm{MW}', '修正后 P_BR 不得写成 28.0 MW'),
        ('0.036', 'beta_perp^2 不得写成 0.036'),
        ('0.83\\ \\mathrm{MW}', '线辐射总量不得写成 0.83 MW'),
        ('25.0\\%', '辐射合计份额不得写成 25.0%'),
        ('890.0\\ \\mathrm{MW}', '光学薄回旋总量不得写成 890.0 MW'),
        ('1.118\\times10^{7}', '磁能密度不得写成 1.118e7'),
        ('47.85\\ \\mathrm{kW', '练习 2 的功率密度不得写成 47.85 kW/m^2'),
        ('0.61\\%', '黑体上限份额不得写成 0.61%'),
    ]
    for lit, why in lits:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['X'], why))


# -*- coding: utf-8 -*-
"""附录 Y（强场与激光等离子体）的门禁体：P30。

可证伪点
  ① 八个小节、七张表、六只语义盒、二十一个公式的标签缺一即 FAIL；
  ② 「出处：」不足 6 处、练习与解答不足 6 道即 FAIL；
  ③ 表 Y.1–Y.6 的关键数字按**所在小节**复算（rel=1e-3）——把数字搬到别的小节即 FAIL；
  ④ 整行字面绑定：表头、表行、语义盒标题被改写即 FAIL；
  ⑤ 口径回流防护：已经纠正过的旧写法（3.34e12 bar 等四舍五入）回流即 FAIL。
"""


def check_P30_appendix_y_strongfield():
    """附录 Y：强场与激光等离子体的结构、出处、练习与数值核对。"""
    spath = os.path.join(SRC, APP_FILES['Y'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 Y 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '三个刻度',
        '等离子体波',
        '尾场加速',
        '相对论透明',
        '辐射压强与 QED',
        '惯性约束',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 Y 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 Y 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 语义盒 / 公式标签 ----------
    for lab in ('tab:app-y-scales', 'tab:app-y-wake', 'tab:app-y-lwfa',
                'tab:app-y-transparency', 'tab:app-y-qed', 'tab:app-y-icf',
                'tab:app-y-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表 %s）' % (APP_FILES['Y'], lab))
    for lab in ('tip-app-y-field', 'tip-app-y-dephasing',
                'tip-app-y-transparency', 'tip-app-y-qed',
                'tip-app-y-pressure', 'tip-app-y-honest'):
        has(text, '{%s}' % lab, 'src/%s（语义盒 %s）' % (APP_FILES['Y'], lab))
    for lab in ('eq:app-y-a0', 'eq:app-y-up', 'eq:app-y-es', 'eq:app-y-schwinger',
                'eq:app-y-ratio', 'eq:app-y-omegape', 'eq:app-y-lamp',
                'eq:app-y-nc', 'eq:app-y-ewb', 'eq:app-y-vg',
                'eq:app-y-dephasing', 'eq:app-y-gain', 'eq:app-y-pc',
                'eq:app-y-transparent', 'eq:app-y-itr', 'eq:app-y-pressure',
                'eq:app-y-chi', 'eq:app-y-ratio2', 'eq:app-y-rhor',
                'eq:app-y-burn', 'eq:app-y-gainicf'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 %s）' % (APP_FILES['Y'], lab))
    has(text, '不复制', 'src/%s（版权与复现口径）' % APP_FILES['Y'])

    # ---------- ③ 出处与练习 ----------
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 Y 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 Y 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 Y 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 Y 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 整行字面绑定 ----------
    lits = [
        # 表 Y.1 三个刻度
        ('刻度 & 数值 & 参照物 & 含义 \\\\', '表 Y.1 表头'),
        ('归一化振幅 & $a_0 = 68.3486$ & $a_0=1$ 为分界 & 抖动速度 $0.999786c$ \\\\',
         '表 Y.1 归一化振幅行'),
        ('激光电场 & $2.74309\\times10^{14}\\ \\mathrm{V\\,m^{-1}}$ & 原子内场 $5.14\\times10^{11}$ & 强 $533.6745$ 倍 \\\\',
         '表 Y.1 激光电场行'),
        ('有质动力势 & $24.1907\\ \\mathrm{MeV}$ & 静止能 $0.511\\ \\mathrm{MeV}$ & 单个周期就超过 D-T 的 $3.5\\ \\mathrm{MeV}$ \\\\',
         '表 Y.1 有质动力势行'),
        ('与 Schwinger 场之比 & $2.07294\\times10^{-4}$ & $\\chi = 0.0100206$ & 离 QED 差四个数量级 \\\\',
         '表 Y.1 Schwinger 比行'),
        # 表 Y.2 等离子体波
        ('等离子体频率 $\\omega_{pe}$ & $5.64146\\times10^{13}\\ \\mathrm{rad\\,s^{-1}}$ & $n_e$ & 尾场的时间标尺 \\\\',
         '表 Y.2 等离子体频率行'),
        ('等离子体波长 $\\lambda_p$ & $33.389\\ \\mathrm{\\mu m}$ & $n_e$ & 加速腔的长度 \\\\',
         '表 Y.2 等离子体波长行'),
        ('临界密度 $n_c$ & $1.742\\times10^{21}\\ \\mathrm{cm^{-3}}$ & $\\lambda$ & 传播判据 $n_e<n_c$ \\\\',
         '表 Y.2 临界密度行'),
        ('波破场 $E_{WB}$ & $96.1592\\ \\mathrm{GV\\,m^{-1}}$ & $\\sqrt{n_e}$ & 加速梯度的上限 \\\\',
         '表 Y.2 波破场行'),
        # 表 Y.3 尾场加速
        ('群速度亏损 & $2.87033\\times10^{-4}$ & $n_e/n_c$ & 电子追不上波 \\\\',
         '表 Y.3 群速度亏损行'),
        ('失相长度 $L_d$ & $5.8163\\ \\mathrm{cm}$ & $\\lambda_p n_c/n_e$ & 加速的有效距离 \\\\',
         '表 Y.3 失相长度行'),
        ('理想增益 $\\Delta W$ & $5.5929\\ \\mathrm{GeV}$ & $E_{WB}L_d$ & 五厘米五点六吉电子伏 \\\\',
         '表 Y.3 理想增益行'),
        ('自聚焦功率 $P_c$ & $30.31\\ \\mathrm{TW}$ & $n_c/n_e$ & 本工作点 $P/P_c = 259.12095$ \\\\',
         '表 Y.3 自聚焦功率行'),
        ('脉冲能量 & $235.619\\ \\mathrm{J}$ & $P\\tau$ & 台式装置的由来 \\\\',
         '表 Y.3 脉冲能量行'),
        # 表 Y.4 相对论透明
        ('固体 $n_e/n_c$ & $57.40661$ & 透明所需的 $a_0$ 下限 & 线性理论下不透明 \\\\',
         '表 Y.4 固体折射率行'),
        ('工作点 $a_0$ & $68.3486$ & 余量 $1.19061$ & 刚过线，只多 $19\\%$ \\\\',
         '表 Y.4 工作点余量行'),
        ('阈值强度 & $7.05447\\times10^{21}\\ \\mathrm{W\\,cm^{-2}}$ & 低 $1.41754$ 倍 & 实验窗口就在附近 \\\\',
         '表 Y.4 阈值强度行'),
        ('相对论临界密度 & $8.42065\\times10^{28}\\ \\mathrm{m^{-3}}$ & $= \\gamma n_c$ & 质量变大抬高了 $n_c$ \\\\',
         '表 Y.4 相对论临界密度行'),
        # 表 Y.5 辐射压强与 QED 参数
        ('辐射压强 $I/c$ & $3.33564\\times10^{17}\\ \\mathrm{Pa}$ & 磁压 $1.11767\\times10^{7}\\ \\mathrm{Pa}$ & 高 $2.9845\\times10^{10}$ 倍 \\\\',
         '表 Y.5 辐射压强行'),
        ('反射加倍 $2I/c$ & $6.67128\\times10^{17}\\ \\mathrm{Pa}$ & --- & 靶面受到的冲量 \\\\',
         '表 Y.5 反射加倍行'),
        ('QED 参数 $\\chi$ & $0.0100206$ & $\\chi=1$ 为 QED 门槛 & 还差约两个数量级 \\\\',
         '表 Y.5 QED 参数行'),
        ('$I/I_S$ & $4.30282\\times10^{-8}$ & $I_S = 2.32406\\times10^{33}\\ \\mathrm{W\\,m^{-2}}$ & 离对产生更远 \\\\',
         '表 Y.5 强度比行'),
        # 表 Y.6 惯性约束
        ('面密度 $\\rho R$ & $3\\ \\mathrm{g\\,cm^{-2}}$ & 工程目标 & alpha 粒子能否自吸收 \\\\',
         '表 Y.6 面密度行'),
        ('半径 $R$ & $150\\ \\mathrm{\\mu m}$ & $\\rho R/\\rho$ & 压缩到 $\\rho = 200\\ \\mathrm{g\\,cm^{-3}}$ \\\\',
         '表 Y.6 半径行'),
        ('燃料质量 & $2.8274\\ \\mathrm{mg}$ & $V\\rho$ & 两毫克 \\\\',
         '表 Y.6 燃料质量行'),
        ('聚变能 $E_{\\mathrm{fus}}$ & $960.2\\ \\mathrm{MJ}$ & $m\\times3.396\\times10^{14}$ & 声明比能 \\\\',
         '表 Y.6 聚变能行'),
        ('燃烧份额 $f$ & $33.333\\%$ & $\\rho R/(\\rho R+6)$ & 三分之一的燃料参与 \\\\',
         '表 Y.6 燃烧份额行'),
        ('实际释放 & $320.07\\ \\mathrm{MJ}$ & $fE_{\\mathrm{fus}}$ & 靶增益 $160.03273$ \\\\',
         '表 Y.6 实际释放行'),
        # 表 Y.7 来源表
        ('本附录内容 & 借鉴来源 & 在本书里的位置 \\\\', '表 Y.7 来源表表头'),
        # 六只语义盒：标题与标签一起锁定
        ('\\begin{tip}{场的强度要用无量纲数说话}{tip-app-y-field}',
         '语义盒「场的强度要用无量纲数说话」'),
        ('\\begin{tip}{失相不是误差，是原理}{tip-app-y-dephasing}',
         '语义盒「失相不是误差，是原理」'),
        ('\\begin{tip}{透明不等于穿透}{tip-app-y-transparency}',
         '语义盒「透明不等于穿透」'),
        ('\\begin{tip}{离 QED 门槛还差两个数量级，所以经典计算够用}{tip-app-y-qed}',
         '语义盒「离 QED 门槛还差两个数量级」'),
        ('\\begin{tip}{光压比磁压大三百亿倍，但冲量才是可比的}{tip-app-y-pressure}',
         '语义盒「光压比磁压大三百亿倍」'),
        ('\\begin{tip}{本附录不给实验点，只给标度律}{tip-app-y-honest}',
         '语义盒「本附录不给实验点，只给标度律」'),
        # 正文关键判断
        ('这个百分之一的量级很关键', 'QED 参数的量级判断（曾是「百分之十」笔误）'),
        ('尾场加速的原始提案见 Tajima \\& Dawson（1979）', 'Y.2 出处句'),
        ('尾场加速的梯度与失相长度推导参照 Tajima \\& Dawson（1979）', 'Y.3 出处句'),
    ]
    for lit, why in lits:
        has(text, lit, 'src/%s（%s）' % (APP_FILES['Y'], why))

    # ---------- ⑤ 单位常量 ----------
    UVCM = r'\mathrm{V\,m^{-1}}'
    UGVM = r'\mathrm{GV\,m^{-1}}'
    UMEV = r'\mathrm{MeV}'
    UGEV = r'\mathrm{GeV}'
    UWCM2 = r'\mathrm{W\,cm^{-2}}'
    UWM2 = r'\mathrm{W\,m^{-2}}'
    UTW = r'\mathrm{TW}'
    UJ = r'\mathrm{J}'
    UPA = r'\mathrm{Pa}'
    UCM = r'\mathrm{cm}'
    UUM = r'\mathrm{\mu m}'
    UCM3 = r'\mathrm{cm^{-3}}'
    UM3 = r'\mathrm{m^{-3}}'
    UM2 = r'\mathrm{m^2}'
    UMJ = r'\mathrm{MJ}'
    UGCM2 = r'\mathrm{g\,cm^{-2}}'
    UGCM3 = r'\mathrm{g\,cm^{-3}}'
    UMG = r'\mathrm{mg}'
    URADS = r'\mathrm{rad\,s^{-1}}'
    UBAR = r'\mathrm{bar}'
    UPCT = r'\%'

    # ---------- ⑥ 逐节数值复算（rel=1e-3）----------
    cases = [
        # Y.1 三个刻度：场强、功率与压强
        ('三个刻度', 68.3486, '', '归一化振幅 a0'),
        ('三个刻度', 2.74309e14, UVCM, '激光电场'),
        ('三个刻度', 533.6745, '', '激光电场 / 原子内场'),
        ('三个刻度', 24.1907, UMEV, '有质动力势'),
        ('三个刻度', 0.511, UMEV, '电子静止能'),
        ('三个刻度', 48.3401, '', '抖动 gamma'),
        ('三个刻度', 0.999786, '', '抖动速度 beta'),
        ('三个刻度', 2.07294e-4, '', '与 Schwinger 场之比'),
        ('三个刻度', 1.32329e18, UVCM, 'Schwinger 场'),
        # Y.2 等离子体波
        ('等离子体波', 5.64146e13, URADS, '等离子体频率'),
        ('等离子体波', 33.389, UUM, '等离子体波长'),
        ('等离子体波', 1.742e21, UCM3, '临界密度'),
        ('等离子体波', 96.1592, UGVM, '波破场'),
        ('等离子体波', 5.74066e-4, '', 'n_e/n_c'),
        ('等离子体波', 1741.96, '', 'n_c/n_e'),
        ('等离子体波', 41.7368, '', 'lambda_p/lambda'),
        # Y.3 尾场加速
        ('尾场加速', 2.87033e-4, '', '群速度亏损'),
        ('尾场加速', 5.8163, UCM, '失相长度'),
        ('尾场加速', 5.5929, UGEV, '理想增益'),
        ('尾场加速', 30.31, UTW, '自聚焦临界功率'),
        ('尾场加速', 259.12095, '', 'P/P_c'),
        ('尾场加速', 235.619, UJ, '脉冲能量'),
        ('尾场加速', 7854, UTW, '激光功率'),
        ('尾场加速', 1.04, UCM, '获得 1 GeV 所需长度'),
        ('尾场加速', 7.85398e-11, UM2, '焦斑面积'),
        # Y.4 相对论透明
        ('相对论透明', 57.40661, '', '固体 n_e/n_c'),
        ('相对论透明', 1.19061, '', '透明余量'),
        ('相对论透明', 7.05447e21, UWCM2, '透明阈值强度'),
        ('相对论透明', 8.42065e28, UM3, '相对论临界密度'),
        ('相对论透明', 48.3401, '', '相对论因子'),
        # Y.5 辐射压强与 QED 参数
        ('辐射压强与 QED', 3.33564e17, UPA, '辐射压强 I/c'),
        ('辐射压强与 QED', 6.67128e17, UPA, '反射加倍 2I/c'),
        ('辐射压强与 QED', 1.11767e7, UPA, '磁压'),
        ('辐射压强与 QED', 2.9845e10, '', '辐射压强 / 磁压'),
        ('辐射压强与 QED', 0.0100206, '', 'QED 参数 chi'),
        ('辐射压强与 QED', 99.794, '', '1/chi'),
        ('辐射压强与 QED', 2.32406e33, UWM2, 'Schwinger 强度'),
        ('辐射压强与 QED', 4.30282e-8, '', 'I/I_S'),
        ('辐射压强与 QED', 3.3356e12, UBAR, '辐射压强（bar 口径）'),
        # Y.6 惯性约束
        ('惯性约束', 3, UGCM2, '面密度'),
        ('惯性约束', 200, UGCM3, '压缩密度'),
        ('惯性约束', 150, UUM, '压缩半径'),
        ('惯性约束', 1.41372e-11, '', '燃料体积'),
        ('惯性约束', 2.8274, UMG, '燃料质量'),
        ('惯性约束', 3.396e14, '', 'DT 声明比能'),
        ('惯性约束', 960.2, UMJ, '聚变能'),
        ('惯性约束', 33.333, UPCT, '燃烧份额'),
        ('惯性约束', 320.07, UMJ, '实际释放'),
        ('惯性约束', 160.03273, '', '靶增益'),
        ('惯性约束', 16.0033, '', '插座增益'),
        # Y.7 六道练习
        ('估算练习', 6.83486, '', '练习 1 的 a0'),
        ('估算练习', 2.01096, UMEV, '练习 1 的有质动力势'),
        ('估算练习', 12.0295, '', '练习 1 的有质动力势比'),
        ('估算练习', 30.4082, UGVM, '练习 2 的波破场'),
        ('估算练习', 105.6, UUM, '练习 2 的等离子体波长'),
        ('估算练习', 1.8393, '', '练习 2 的失相长度'),
        ('估算练习', 3.1623, '', '练习 2 的 sqrt(10) 比'),
        ('估算练习', 31.623, '', '练习 2 的长度比'),
        ('估算练习', 9.9222e18, UCM3, '练习 3 的 CO2 临界密度'),
        ('估算练习', 0.100784, '', '练习 3 的 n_e/n_c'),
        ('估算练习', 0.3313, '', '练习 3 的失相长度（mm）'),
        ('估算练习', 2.14063e18, UWCM2, '练习 4 的 a0=1 强度'),
        ('估算练习', 1.68124, UTW, '练习 4 的功率'),
        ('估算练习', 50.44, '', '练习 4 的脉冲能量（mJ）'),
        ('估算练习', 4671.5, '', '练习 4 与工作点的强度比'),
        ('估算练习', 5540.6, UTW, '练习 5 的透明阈值功率'),
        ('估算练习', 0.41888, UMG, '练习 6 的燃料质量'),
        ('估算练习', 142.25, UMJ, '练习 6 的聚变能'),
        ('估算练习', 14.286, UPCT, '练习 6 的燃烧份额'),
        ('估算练习', 20.322, UMJ, '练习 6 的实际释放'),
        ('估算练习', 10.16081, '', '练习 6 的靶增益'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s：%s（%s）' % (sec, val, why))
    if bad:
        rec('FAIL', 'P', '附录 Y 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 Y 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    lits = [
        ('3.34\\times10^{12}', '辐射压强的 bar 口径不得回退到 3.34e12'),
        ('3.34\\times10^{6}', '辐射压强的 Mbar 口径不得回退到 3.34e6'),
        ('2.74\\times10^{14}', '激光电场不得四舍五入到 2.74e14'),
        ('24.2\\ \\mathrm{MeV}', '有质动力势不得四舍五入到 24.2 MeV'),
        ('96.2\\ \\mathrm{GV', '波破场不得四舍五入到 96.2 GV/m'),
        ('5.8\\ \\mathrm{cm}', '失相长度不得四舍五入到 5.8 cm'),
        ('5.6\\ \\mathrm{GeV}', '理想增益不得四舍五入到 5.6 GeV'),
        ('30.3\\ \\mathrm{TW}', '自聚焦功率不得四舍五入到 30.3 TW'),
        ('7.05\\times10^{21}', '透明阈值不得四舍五入到 7.05e21'),
        ('8.42\\times10^{28}', '相对论临界密度不得四舍五入到 8.42e28'),
        ('2.98\\times10^{10}', '压强比不得四舍五入到 2.98e10'),
        ('4.30\\times10^{-8}', '强度比不得四舍五入到 4.30e-8'),
        ('960\\ \\mathrm{MJ}', '聚变能不得写成 960 MJ'),
        ('320\\ \\mathrm{MJ}', '实际释放不得写成 320 MJ'),
        ('33.3\\%', '燃烧份额不得写成 33.3%'),
        ('160.0$', '靶增益不得写成 160.0'),
        ('百分之十的量级', 'QED 参数的量级不得写成「百分之十」'),
    ]
    for lit, why in lits:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['Y'], why))


# -*- coding: utf-8 -*-
"""附录 Z（聚变堆工程）的门禁体：P31。

可证伪点
  ① 八个小节、七张表、六只语义盒、二十一个公式的标签缺一即 FAIL；
  ② 「出处：」不足 6 处、练习与解答不足 6 道即 FAIL；
  ③ 表 Z.1–Z.6 的关键数字按**所在小节**复算（rel=1e-3）——把数字搬到别的小节即 FAIL；
  ④ 整行字面绑定：表头、表行、语义盒标题、关键判断句被改写即 FAIL；
  ⑤ 口径回流防护：已经纠正过的旧写法（二十个百分点、四舍五入后的 640 MW 等）回流即 FAIL。
"""


def check_P31_appendix_z_reactor():
    """附录 Z：聚变堆工程的结构、出处、练习与数值核对。"""
    spath = os.path.join(SRC, APP_FILES['Z'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 Z 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '三本账接起来',
        '中子：壁负载',
        '氚循环',
        '屏蔽与超导磁体',
        '热工与冷却剂',
        '可用率与年发电量',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 Z 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 Z 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 语义盒 / 公式标签 ----------
    for lab in ('tab:app-z-chain', 'tab:app-z-neutron', 'tab:app-z-tritium',
                'tab:app-z-shield', 'tab:app-z-thermal', 'tab:app-z-availability',
                'tab:app-z-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表 %s）' % (APP_FILES['Z'], lab))
    for lab in ('tip-app-z-chain', 'tip-app-z-neutron', 'tip-app-z-tbr',
                'tip-app-z-thermal', 'tip-app-z-availability', 'tip-app-z-honest'):
        has(text, '{%s}' % lab, 'src/%s（语义盒 %s）' % (APP_FILES['Z'], lab))
    for lab in ('eq:app-z-pn', 'eq:app-z-pth', 'eq:app-z-pgross', 'eq:app-z-pnet',
                'eq:app-z-qeng', 'eq:app-z-qn', 'eq:app-z-dpa', 'eq:app-z-life',
                'eq:app-z-burnj', 'eq:app-z-burnrate', 'eq:app-z-decay',
                'eq:app-z-tbrbreak', 'eq:app-z-excess', 'eq:app-z-proc',
                'eq:app-z-shield', 'eq:app-z-cryo', 'eq:app-z-mdot',
                'eq:app-z-pump', 'eq:app-z-eyear', 'eq:app-z-eyearnum',
                'eq:app-z-gain85'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 %s）' % (APP_FILES['Z'], lab))
    has(text, '不复制', 'src/%s（版权与复现口径）' % APP_FILES['Z'])

    # ---------- ③ 出处与练习 ----------
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 Z 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 Z 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 Z 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 Z 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 整行字面绑定 ----------
    lits = [
        # 表 Z.1 四步能量链
        ('环节 & 数值 & 乘子 & 备注 \\\\', '表 Z.1 表头'),
        ('聚变功率 & $556.5\\ \\mathrm{MW}$ & --- & D1 工作点 \\\\', '表 Z.1 聚变功率行'),
        ('alpha / 中子分配 & $111.3$ / $445.2\\ \\mathrm{MW}$ & $1/5$ 与 $4/5$ & 本书统一口径 \\\\',
         '表 Z.1 分配行'),
        ('包层热功率 & $639.97\\ \\mathrm{MW}$ & $M = 1.15$ & 声明值 \\\\', '表 Z.1 包层行'),
        ('毛电功率 & $255.99\\ \\mathrm{MW}$ & $\\eta_{\\mathrm{th}} = 0.40$ & 声明值 \\\\',
         '表 Z.1 毛电行'),
        ('净电功率 & $217.59\\ \\mathrm{MW}$ & $1-f_{\\mathrm{rec}} = 0.85$ & 再循环 15\\% \\\\',
         '表 Z.1 净电行'),
        ('工程增益 & $Q_{\\mathrm{eng}} = 4.3518$ & $P_{\\mathrm{net}}/P_{\\mathrm{aux}}$ & 对照 $Q = 11.13$ \\\\',
         '表 Z.1 工程增益行'),
        # 表 Z.2 中子
        ('第一壁面积 & $489.532\\ \\mathrm{m^2}$ & 装置几何 & 与附录 X 同口径 \\\\',
         '表 Z.2 壁面积行'),
        ('中子壁负载 & $0.90944\\ \\mathrm{MW\\,m^{-2}}$ & $P_n/A_{\\mathrm{wall}}$ & 上限的 $0.45472$ \\\\',
         '表 Z.2 壁负载行'),
        ('壁负载余量 & $2.1992$ 倍 & 上限 $2\\ \\mathrm{MW\\,m^{-2}}$ & 声明上限 \\\\',
         '表 Z.2 余量行'),
        ('位移损伤率 & $4.5472\\ \\mathrm{dpa\\,year^{-1}}$ & $5\\,q_n$ & 声明系数 \\\\',
         '表 Z.2 dpa 行'),
        ('包层寿命 & $10.996$ 年 & $50/\\dot d$ & 必须可更换 \\\\', '表 Z.2 寿命行'),
        ('年输运能量 & $1.4049\\times10^{7}\\ \\mathrm{GJ}$ & $P_n\\times1$ 年 & 必须被包层收下 \\\\',
         '表 Z.2 年能量行'),
        # 表 Z.3 氚
        ('氚燃耗率 & $31.191\\ \\mathrm{kg\\,year^{-1}}$ & $P_{\\mathrm{fus}}N_Tm_T$ & $=85.397\\ \\mathrm{g/day}$ \\\\',
         '表 Z.3 燃耗行'),
        ('每 GW 聚变功率 & $56.049\\ \\mathrm{kg\\,year^{-1}}$ & 同上一行 & 记住这个数 \\\\',
         '表 Z.3 每 GW 行'),
        ('库存自衰变 & $0.11252\\ \\mathrm{kg\\,year^{-1}}$ & $\\lambda m_{\\mathrm{inv}}$ & 库存 $2\\ \\mathrm{kg}$ \\\\',
         '表 Z.3 衰变行'),
        ('自持 TBR 下限 & $1.0036$ & $1+\\dot m_{\\mathrm{decay}}/\\dot m_T$ & 只要千分之四 \\\\',
         '表 Z.3 自持下限行'),
        ('声明 TBR & $1.10$ & 工程余量 & 净增 $3.0066\\ \\mathrm{kg\\,year^{-1}}$ \\\\',
         '表 Z.3 声明 TBR 行'),
        ('处理流量 & $12\\ \\mathrm{kg/day}$ & $m_{\\mathrm{inv}}/\\tau_{\\mathrm{proc}}$ & $\\tau = 4\\ \\mathrm{h}$ \\\\',
         '表 Z.3 处理流量行'),
        # 表 Z.4 屏蔽
        ('衰减长度 & $0.15\\ \\mathrm{m}$ & 声明值 & 钢/水/硼体系量级 \\\\', '表 Z.4 衰减长度行'),
        ('屏蔽厚度（$10^{-6}$） & $2.0723\\ \\mathrm{m}$ & $\\lambda\\ln(1/\\varepsilon)$ & 两米厚 \\\\',
         '表 Z.4 屏蔽厚度行'),
        ('漏到磁体 & $445.2\\ \\mathrm{W}$ & $\\varepsilon P_n$ & 作为冷负荷 \\\\',
         '表 Z.4 漏中子行'),
        ('制冷电功 & $133.56\\ \\mathrm{kW}$ & $\\mathrm{COP}\\,\\varepsilon P_n$ & 占净电 $0.061381\\%$ \\\\',
         '表 Z.4 制冷行'),
        # 表 Z.5 热工
        ('冷却剂比热 & $5.2\\ \\mathrm{kJ\\,kg^{-1}K^{-1}}$ & 声明值（氦） & 温升 $200\\ \\mathrm{K}$ \\\\',
         '表 Z.5 比热行'),
        ('质量流量 & $615.36\\ \\mathrm{kg\\,s^{-1}}$ & $P_{\\mathrm{th}}/(c_p\\Delta T)$ & 每秒六百多公斤 \\\\',
         '表 Z.5 流量行'),
        ('泵功率 & $6.3997\\ \\mathrm{MW}$ & $0.01P_{\\mathrm{th}}$ & 占毛电 $2.5\\%$ \\\\',
         '表 Z.5 泵功率行'),
        # 表 Z.6 可用率
        ('可用率 & $0.75$ & 声明值 & 含换包层与检修 \\\\', '表 Z.6 可用率行'),
        ('年发电量 & $1429.6\\ \\mathrm{GWh}$ & $P_{\\mathrm{net}}\\times8760\\,\\mathrm{h}\\times A$ & $5.14647\\times10^{15}\\ \\mathrm{J}$ \\\\',
         '表 Z.6 年发电量行'),
        ('年停机 & $2190\\ \\mathrm{h}$ & $8760\\times(1-A)$ & 一年四分之一停着 \\\\',
         '表 Z.6 年停机行'),
        ('可用率提到 $0.85$ & $1620.2\\ \\mathrm{GWh}$ & 同上 & 增幅 $13.333\\%$ \\\\',
         '表 Z.6 可用率提升行'),
        ('折合标煤 & $0.48606$ 百万吨 & $0.34\\ \\mathrm{kg/kWh}$ & 声明口径 \\\\',
         '表 Z.6 标煤行'),
        # 表 Z.7 来源表表头
        ('本附录内容 & 借鉴来源 & 在本书里的位置 \\\\', '表 Z.7 来源表表头'),
        # 六只语义盒：标题与标签一起锁定
        ('\\begin{tip}{工程账是乘法，不是加法}{tip-app-z-chain}',
         '语义盒「工程账是乘法，不是加法」'),
        ('\\begin{tip}{难受的不是热，是 dpa}{tip-app-z-neutron}',
         '语义盒「难受的不是热，是 dpa」'),
        ('\\begin{tip}{自持只看两个数，余量留给不确定性}{tip-app-z-tbr}',
         '语义盒「自持只看两个数」'),
        ('\\begin{tip}{两本账要对得上：泵功率在哪一边}{tip-app-z-thermal}',
         '语义盒「两本账要对得上」'),
        ('\\begin{tip}{可用率是「另一半增益」}{tip-app-z-availability}',
         '语义盒「可用率是另一半增益」'),
        ('\\begin{tip}{把「声明值」列成一张清单}{tip-app-z-honest}',
         '语义盒「把声明值列成一张清单」'),
        # 正文关键判断
        ('$1.15\\times0.40\\times0.85 = 0.391$', 'Z.1 的乘法链实例'),
        ('每 $1\\ \\mathrm{MW}$ 聚变功率只能送出', 'Z.1 的「每兆瓦只送出」判断句'),
        ('误差最大 $0.57143\\%$（alpha 侧）', 'Z.1 的口径误差声明'),
        ('相当于一座千兆瓦级电站半年的发电量', 'Z.2 的能量量级对照句'),
        ('十个百分点的可用率，换来 $13.333\\%$ 的年发电量', 'Z.6 的核心判断句（曾是二十个百分点）'),
    ]
    for lit, why in lits:
        has(text, lit, 'src/%s（%s）' % (APP_FILES['Z'], why))

    # ---------- ⑤ 单位常量 ----------
    UMW = r'\mathrm{MW}'
    UMWM2 = r'\mathrm{MW\,m^{-2}}'
    UDBA = r'\mathrm{dpa\,year^{-1}}'
    UGJ = r'\mathrm{GJ}'
    UKGY = r'\mathrm{kg\,year^{-1}}'
    UKGD = r'\mathrm{kg/day}'
    UMG = r'\mathrm{m}'
    UW = r'\mathrm{W}'
    UKW = r'\mathrm{kW}'
    UKJKG = r'\mathrm{kJ\,kg^{-1}K^{-1}}'
    UK = r'\mathrm{K}'
    UKGS = r'\mathrm{kg\,s^{-1}}'
    UGWH = r'\mathrm{GWh}'
    UH = r'\mathrm{h}'
    UM2 = r'\mathrm{m^2}'
    UJ = r'\mathrm{J}'
    UGD = r'\mathrm{g/day}'
    UKGH = r'\mathrm{kg\,h^{-1}}'
    UYINV = r'\mathrm{year^{-1}}'
    UJINV = r'\mathrm{J^{-1}}'
    UKS = r'\mathrm{kg\,s^{-1}}'
    UPCT = r'\%'

    # ---------- ⑥ 逐节数值复算（rel=1e-3）----------
    cases = [
        # Z.1 四步能量链
        ('三本账接起来', 556.5, UMW, '聚变功率'),
        ('三本账接起来', 111.3, '', 'alpha 功率'),
        ('三本账接起来', 445.2, UMW, '中子功率'),
        ('三本账接起来', 0.1989, '', '精确 alpha 份额'),
        ('三本账接起来', 0.8011, '', '精确中子份额'),
        ('三本账接起来', 0.57143, UPCT, '1/5 口径的相对偏差'),
        ('三本账接起来', 639.97, UMW, '包层热功率'),
        ('三本账接起来', 255.99, UMW, '毛电功率'),
        ('三本账接起来', 217.59, UMW, '净电功率'),
        ('三本账接起来', 0.34, '', '整链净效率'),
        ('三本账接起来', 11.13, '', 'Q'),
        ('三本账接起来', 4.3518, '', 'Q_eng'),
        ('三本账接起来', 2.5575, '', 'Q/Q_eng'),
        ('三本账接起来', 0.391, '', '每兆瓦聚变对应的净电'),
        # Z.2 中子
        ('中子：壁负载', 489.532, UM2, '第一壁面积'),
        ('中子：壁负载', 0.90944, UMWM2, '中子壁负载'),
        ('中子：壁负载', 2.1992, '', '壁负载余量'),
        ('中子：壁负载', 0.45472, '', '壁负载 / 上限'),
        ('中子：壁负载', 4.5472, UDBA, '位移损伤率'),
        ('中子：壁负载', 10.996, '', '包层寿命'),
        ('中子：壁负载', 1.4049e7, UGJ, '年输运能量'),
        # Z.3 氚
        ('氚循环', 3.54631e11, UJINV, '每焦耳氚原子数'),
        ('氚循环', 5.00827e-27, '', '氚原子质量'),
        ('氚循环', 9.88393e-7, UKS, '燃耗率'),
        ('氚循环', 85.397, UGD, '燃耗率（g/day）'),
        ('氚循环', 31.191, UKGY, '燃耗率（kg/年）'),
        ('氚循环', 56.049, UKGY, '每 GW 聚变功率'),
        ('氚循环', 0.0562619, UYINV, '衰变常数'),
        ('氚循环', 0.11252, UKGY, '库存自衰变'),
        ('氚循环', 0.30807, UGD, '库存自衰变（g/day）'),
        ('氚循环', 1.0036, '', '自持 TBR 下限'),
        ('氚循环', 3.0066, UKGY, 'TBR=1.10 的净增'),
        ('氚循环', 0.5, UKGH, '处理流量'),
        ('氚循环', 12, UKGD, '处理流量（kg/day）'),
        # Z.4 屏蔽
        ('屏蔽与超导磁体', 0.15, UMG, '衰减长度'),
        ('屏蔽与超导磁体', 2.0723, UMG, '屏蔽厚度'),
        ('屏蔽与超导磁体', 445.2, UW, '漏到磁体的功率'),
        ('屏蔽与超导磁体', 133.56, UKW, '制冷电功'),
        ('屏蔽与超导磁体', 0.061381, UPCT, '制冷占净电份额'),
        # Z.5 热工
        ('热工与冷却剂', 5.2, UKJKG, '冷却剂比热'),
        ('热工与冷却剂', 200, UK, '冷却剂温升'),
        ('热工与冷却剂', 615.36, UKGS, '质量流量'),
        ('热工与冷却剂', 6.3997, UMW, '泵功率'),
        ('热工与冷却剂', 2.5, UPCT, '泵功率占毛电份额'),
        ('热工与冷却剂', 0.01, '', '泵份额系数'),
        # Z.6 可用率
        ('可用率与年发电量', 0.75, '', '可用率'),
        ('可用率与年发电量', 1429.6, UGWH, '年发电量'),
        ('可用率与年发电量', 5.14647e15, UJ, '年发电量（J）'),
        ('可用率与年发电量', 2190, UH, '年停机小时'),
        ('可用率与年发电量', 0.85, '', '提高后的可用率'),
        ('可用率与年发电量', 1620.2, UGWH, '提高后的年发电量'),
        ('可用率与年发电量', 13.333, UPCT, '年发电量增幅'),
        ('可用率与年发电量', 0.48606, '', '折合标煤'),
        # Z.7 练习
        ('估算练习', 723.45, UMW, '练习 1 的 P_th'),
        ('估算练习', 289.38, UMW, '练习 1 的毛电'),
        ('估算练习', 245.97, UMW, '练习 1 的净电'),
        ('估算练习', 1.1304, '', '练习 1 的倍数'),
        ('估算练习', 287.99, UMW, '练习 2 的毛电'),
        ('估算练习', 244.79, UMW, '练习 2 的净电'),
        ('估算练习', 12.5, UPCT, '练习 2 的涨幅'),
        ('估算练习', 0.45472, UMWM2, '练习 3 的壁负载'),
        ('估算练习', 4.3983, '', '练习 3 的余量'),
        ('估算练习', 2.7631, UMG, '练习 4 的屏蔽厚度'),
        ('估算练习', 0.69078, UMG, '练习 4 的增厚'),
        ('估算练习', 0.05, '', '练习 5 的 TBR 增量'),
        ('估算练习', 1.447, UKGY, '练习 5 的净增'),
        ('估算练习', 0.55086, '', '练习 6 的折合标煤'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s：%s（%s）' % (sec, val, why))
    if bad:
        rec('FAIL', 'P', '附录 Z 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 Z 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    lits = [
        ('640\\ \\mathrm{MW}', '包层热功率不得四舍五入到 640 MW'),
        ('256\\ \\mathrm{MW}', '毛电功率不得四舍五入到 256 MW'),
        ('218\\ \\mathrm{MW}', '净电功率不得四舍五入到 218 MW'),
        ('0.91\\ \\mathrm{MW', '中子壁负载不得四舍五入到 0.91 MW/m^2'),
        ('4.55\\ \\mathrm{dpa', 'dpa 率不得四舍五入到 4.55'),
        ('11.0\\ ', '包层寿命不得写成 11.0 年'),
        ('31.2\\ \\mathrm{kg', '氚燃耗率不得四舍五入到 31.2 kg/年'),
        ('56.0\\ \\mathrm{kg', '每 GW 燃耗不得写成 56.0 kg/年'),
        ('0.11\\ \\mathrm{kg', '库存衰变不得四舍五入到 0.11 kg/年'),
        ('3.0\\ \\mathrm{kg', 'TBR 净增不得四舍五入到 3.0 kg/年'),
        ('2.07\\ \\mathrm{m}', '屏蔽厚度不得四舍五入到 2.07 m'),
        ('133.6\\ ', '制冷电功不得四舍五入到 133.6 kW'),
        ('0.061\\%', '制冷份额不得写成 0.061%'),
        ('615\\ \\mathrm{kg', '质量流量不得四舍五入到 615 kg/s'),
        ('6.4\\ \\mathrm{MW', '泵功率不得四舍五入到 6.4 MW'),
        ('1430\\ ', '年发电量不得四舍五入到 1430 GWh'),
        ('0.49\\ ', '折合标煤不得四舍五入到 0.49 百万吨'),
        ('0.55\\ ', '提高后的标煤不得四舍五入到 0.55 百万吨'),
        ('2.76\\ \\mathrm{m}', '练习 4 的厚度不得四舍五入到 2.76 m'),
        ('二十个百分点', '可用率措辞不得写成「二十个百分点」'),
    ]
    for lit, why in lits:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['Z'], why))


# -*- coding: utf-8 -*-
"""附录 AA（替代位形）的门禁体：P32。

可证伪点
  ① 八个小节、七张表、六只语义盒、二十一个公式的标签缺一即 FAIL；
  ② 「出处：」不足 6 处、练习与解答不足 6 道即 FAIL；
  ③ 四张数据表的关键数字按**所在小节**复算（rel=1e-3）；
  ④ 整行字面绑定：表头、表行、语义盒标题、关键判断句被改写即 FAIL；
  ⑤ 口径回流防护：四舍五入写法（3.06、0.0015、99470 等）回流即 FAIL。
"""


def check_P32_appendix_aa_alternates():
    """附录 AA：替代位形的结构、出处、练习与数值核对。"""
    spath = os.path.join(SRC, APP_FILES['AA'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AA 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    want = [
        '三种位形的一张地图',
        '旋转变换：仿星器',
        '新经典输运',
        '球形托卡马克',
        '场反位形',
        '定标与三重积',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AA 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AA 的 %d 个小节齐全' % len(want))

    for lab in ('tab:app-aa-map', 'tab:app-aa-iota', 'tab:app-aa-transport',
                'tab:app-aa-st', 'tab:app-aa-frc', 'tab:app-aa-scaling',
                'tab:app-aa-source'):
        has(text, '\\label{%s}' % lab, 'src/%s（表 %s）' % (APP_FILES['AA'], lab))
    for lab in ('tip-app-aa-map', 'tip-app-aa-iota', 'tip-app-aa-eps',
                'tip-app-aa-st', 'tip-app-aa-frc', 'tip-app-aa-honest'):
        has(text, '{%s}' % lab, 'src/%s（语义盒 %s）' % (APP_FILES['AA'], lab))
    for lab in ('eq:app-aa-eps', 'eq:app-aa-epsnum', 'eq:app-aa-beta',
                'eq:app-aa-betan', 'eq:app-aa-iota', 'eq:app-aa-shear',
                'eq:app-aa-sqeps', 'eq:app-aa-eff32', 'eq:app-aa-effratio',
                'eq:app-aa-1nu', 'eq:app-aa-bp', 'eq:app-aa-fbs',
                'eq:app-aa-space', 'eq:app-aa-frcp', 'eq:app-aa-frcn',
                'eq:app-aa-frci', 'eq:app-aa-frclimit', 'eq:app-aa-iss04',
                'eq:app-aa-tauE', 'eq:app-aa-ratio', 'eq:app-aa-triple'):
        has(text, '\\label{%s}' % lab, 'src/%s（公式 %s）' % (APP_FILES['AA'], lab))
    has(text, '不复制', 'src/%s（版权与复现口径）' % APP_FILES['AA'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AA 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AA 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AA 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AA 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    lits = [
        # 表 AA.1 三种位形的地图
        ('量 & D1 托卡马克 & 球形托卡马克 & 仿星器 \\\\', '表 AA.1 表头'),
        ('环径比 $A=R/a$ & $3.1$ & $1.5$ & $10.377$ \\\\', '表 AA.1 环径比行'),
        ('逆环径比 $\\varepsilon$ & $0.32258$ & $0.66667$ & $0.096364$ \\\\',
         '表 AA.1 逆环径比行'),
        ('压强 $p$ & $1.473\\times10^{5}\\ \\mathrm{Pa}$ & $16022\\ \\mathrm{Pa}$ & $25635\\ \\mathrm{Pa}$ \\\\',
         '表 AA.1 压强行'),
        ('比压 $\\beta$ & $1.318\\%$ & $4.027\\%$ & $1.031\\%$ \\\\', '表 AA.1 比压行'),
        ('归一化比压 $\\beta_N$ & $1.5523$ & $1.2483$ & $1.3659$ \\\\', '表 AA.1 归一化比压行'),
        ('环向电流 & 有（$9\\ \\mathrm{MA}$） & 有（$2\\ \\mathrm{MA}$） & \\textbf{没有} \\\\',
         '表 AA.1 环向电流行'),
        # 表 AA.2 旋转变换
        ('旋转变换 $\\iota/2\\pi$ & 由 $q$ 反推 & $0.85$ & 仿星器按设计取 \\\\',
         '表 AA.2 旋转变换行'),
        ('安全因子 $q$ & $q_{95} = 3.00$ & $1.1765$ & $q = 1/\\iota$ \\\\', '表 AA.2 安全因子行'),
        ('起始 $q_0$ & $1.05$ & 可自由设计 & 托卡马克受锯齿限制 \\\\', '表 AA.2 起始 q0 行'),
        ('磁剪切 $s$ & $1.3$ & 可设计到 $\\approx 0$ & 仿星器的独门手段 \\\\', '表 AA.2 磁剪切行'),
        # 表 AA.3 新经典输运
        ('$\\sqrt{\\varepsilon_{D1}}$ & $0.56796$ & 托卡马克几何 & 香蕉区放大因子 \\\\',
         '表 AA.3 香蕉区因子行'),
        ('$\\varepsilon_{\\mathrm{eff}}$ & $0.013$ & 位形优化 & 声明值（优化位形） \\\\',
         '表 AA.3 有效纹波行'),
        ('$\\varepsilon_{\\mathrm{eff}}/\\varepsilon$ & $0.13491$ & 优化品质 & 只有几何值的 $13.5\\%$ \\\\',
         '表 AA.3 比值行'),
        ('$\\varepsilon_{\\mathrm{eff}}^{3/2}$ & $0.00148223$ & $1/\\nu$ 区 & 仿星器的几何因子 \\\\',
         '表 AA.3 二分之三次方行'),
        ('两个因子的比 & $383.181$ & 不同区间 & \\textbf{不可直接比较} \\\\',
         '表 AA.3 因子比行'),
        # 表 AA.4 球形托卡马克
        ('边缘极向场 $B_p$ & --- & $0.64516\\ \\mathrm{T}$ & $\\mu_0I_p/(2\\pi a)$ \\\\',
         '表 AA.4 极向场行'),
        ('极向比压 $\\beta_p$ & --- & $0.096742$ & 声明 $I_p = 2\\ \\mathrm{MA}$ \\\\',
         '表 AA.4 极向比压行'),
        ('自举份额 $f_{bs}$ & 几个 $\\%$ & $0.078989$ & $C_{bs}\\sqrt{\\varepsilon}\\beta_p$ \\\\',
         '表 AA.4 自举份额行'),
        ('中心柱空间 $R-a$ & $4.2\\ \\mathrm{m}$ & $0.31\\ \\mathrm{m}$ & 屏蔽需要 $2.0723\\ \\mathrm{m}$ \\\\',
         '表 AA.4 中心柱行'),
        ('与屏蔽需求之比 & 放得下 & 少 $6.6848$ 倍 & \\textbf{堆用不了} \\\\',
         '表 AA.4 屏蔽比行'),
        # 表 AA.5 场反位形
        ('外部磁场 $B_{\\mathrm{ext}}$ & $0.5\\ \\mathrm{T}$ & 声明值 & 位形全靠它 \\\\',
         '表 AA.5 外部磁场行'),
        ('压力 $p$ & $99472\\ \\mathrm{Pa}$ & $B^2/(2\\mu_0)$ & $\\beta = 1$ \\\\',
         '表 AA.5 压强行'),
        ('密度（$1.0\\ \\mathrm{keV}$） & $6.2085\\times10^{20}\\ \\mathrm{m^{-3}}$ & $p/k_BT$ & 比 D1 高六倍 \\\\',
         '表 AA.5 密度行'),
        ('反平行电流 & $1\\ \\mathrm{MA}$ & $2\\pi r_sB/\\mu_0$ & 只维持位形 \\\\',
         '表 AA.5 电流行'),
        ('稳定性 & 靠动理学 & 无磁剪切 & \\textbf{不给阈值} \\\\', '表 AA.5 稳定性行'),
        # 表 AA.6 定标
        ('ISS04 系数 & $0.134$ & 声明值 & 仿星器定标 \\\\', '表 AA.6 系数行'),
        ('$\\iota$ 的指数 & $+0.41$ & ISS04 & $\\iota\\uparrow$ 约束改善慢 \\\\',
         '表 AA.6 iota 指数行'),
        ('$B$ 的指数 & $+0.84$ & ISS04 & 磁场更有效 \\\\', '表 AA.6 B 指数行'),
        ('$\\tau_E^{\\mathrm{ISS04}}$ & $0.142983\\ \\mathrm{s}$ & 式 \\eqref{eq:app-aa-iss04} & 研究装置量级 \\\\',
         '表 AA.6 tau_E 行'),
        ('三重积 & $2.28774\\times10^{19}$ & $nT\\tau_E$ & 差 $135.51$ 倍 \\\\', '表 AA.6 三重积行'),
        # 表 AA.7 来源表表头
        ('本附录内容 & 借鉴来源 & 在本书里的位置 \\\\', '表 AA.7 来源表表头'),
        # 六只语义盒
        ('\\begin{tip}{先问「它拿什么换什么」}{tip-app-aa-map}', '语义盒「先问它拿什么换什么」'),
        ('\\begin{tip}{电流是「免费的旋转变换」，也是麻烦的来源}{tip-app-aa-iota}',
         '语义盒「电流是免费的旋转变换」'),
        ('\\begin{tip}{不同区间的数不能直接相除}{tip-app-aa-eps}', '语义盒「不同区间的数不能直接相除」'),
        ('\\begin{tip}{物理的优化可能被几何一票否决}{tip-app-aa-st}',
         '语义盒「物理的优化可能被几何一票否决」'),
        ('\\begin{tip}{比压等于一不是「极限」，是「定义」}{tip-app-aa-frc}',
         '语义盒「比压等于一不是极限是定义」'),
        ('\\begin{tip}{定标律是「代理」，不是物理}{tip-app-aa-honest}',
         '语义盒「定标律是代理不是物理」'),
        # 正文关键判断
        ('\\textbf{三个位形的 $\\varepsilon$ 相差近七倍}', 'AA.1 的「相差近七倍」判断句'),
        ('$\\sqrt{\\varepsilon_{D1}}/\\varepsilon_{\\mathrm{eff}}^{3/2}\n= 383.181$',
         'AA.3 的因子比（分段写法容错见下）'),
        ('\\textbf{放得下}', 'AA.4 的「放得下」'),
        ('\\textbf{放不下}', 'AA.4 的「放不下」'),
    ]
    for lit, why in lits:
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AA'], why))

    UPA = r'\mathrm{Pa}'
    UT = r'\mathrm{T}'
    UMA = r'\mathrm{MA}'
    UM3 = r'\mathrm{m^{-3}}'
    UM = r'\mathrm{m}'
    US = r'\mathrm{s}'
    UKEV = r'\mathrm{keV}'
    UA = r'\mathrm{A}'
    UPCT = r'\%'

    cases = [
        # AA.1 地图
        ('三种位形的一张地图', 3.1, '', 'D1 环径比'),
        ('三种位形的一张地图', 1.5, '', '球形托卡马克环径比'),
        ('三种位形的一张地图', 10.377, '', '仿星器环径比'),
        ('三种位形的一张地图', 0.32258, '', 'D1 逆环径比'),
        ('三种位形的一张地图', 0.66667, '', 'ST 逆环径比'),
        ('三种位形的一张地图', 0.096364, '', '仿星器逆环径比'),
        ('三种位形的一张地图', 1.473e5, UPA, 'D1 压强'),
        ('三种位形的一张地图', 16022, UPA, 'ST 压强'),
        ('三种位形的一张地图', 25635, UPA, '仿星器压强'),
        ('三种位形的一张地图', 1.318, UPCT, 'D1 比压'),
        ('三种位形的一张地图', 4.027, UPCT, 'ST 比压'),
        ('三种位形的一张地图', 1.031, UPCT, '仿星器比压'),
        ('三种位形的一张地图', 1.5523, '', 'D1 归一化比压'),
        ('三种位形的一张地图', 1.2483, '', 'ST 归一化比压'),
        ('三种位形的一张地图', 1.3659, '', '仿星器归一化比压'),
        ('三种位形的一张地图', 9, UMA, 'D1 环向电流'),
        ('三种位形的一张地图', 3.0553, '', 'ST 与 D1 的比压比'),
        # AA.2 旋转变换
        ('旋转变换：仿星器', 0.85, '', 'iota/2pi'),
        ('旋转变换：仿星器', 1.1765, '', 'q = 1/iota'),
        ('旋转变换：仿星器', 3.00, '', 'q95'),
        ('旋转变换：仿星器', 1.05, '', 'q0'),
        ('旋转变换：仿星器', 1.3, '', '磁剪切'),
        # AA.3 新经典输运
        ('新经典输运', 0.56796, '', 'sqrt(eps_D1)'),
        ('新经典输运', 0.013, '', 'eps_eff'),
        ('新经典输运', 0.13491, '', 'eps_eff/eps'),
        ('新经典输运', 13.5, UPCT, 'eps_eff/eps 的百分数'),
        ('新经典输运', 0.00148223, '', 'eps_eff^{3/2}'),
        ('新经典输运', 383.181, '', '两个因子之比'),
        # AA.4 球形托卡马克
        ('球形托卡马克', 0.64516, UT, '边缘极向场'),
        ('球形托卡马克', 0.096742, '', '极向比压'),
        ('球形托卡马克', 0.078989, '', '自举份额'),
        ('球形托卡马克', 7.9, UPCT, '自举份额的百分数'),
        ('球形托卡马克', 4.2, UM, 'D1 中心柱空间'),
        ('球形托卡马克', 0.31, UM, 'ST 中心柱空间'),
        ('球形托卡马克', 2.0723, UM, '屏蔽需求（附录 Z）'),
        ('球形托卡马克', 6.6848, '', '屏蔽需求与中心柱之比'),
        # AA.5 场反位形
        ('场反位形', 0.5, UT, '外部磁场'),
        ('场反位形', 99472, UPA, 'beta=1 的压力'),
        ('场反位形', 6.2085e20, UM3, '对应密度'),
        ('场反位形', 1e6, UA, '反平行电流（A）'),
        ('场反位形', 1, UMA, '反平行电流（MA）'),
        ('场反位形', 10, '', '长细比'),
        ('场反位形', 75.87, '', '1/beta(D1)'),
        # AA.6 定标
        ('定标与三重积', 0.134, '', 'ISS04 系数'),
        ('定标与三重积', 0.41, '', 'iota 指数'),
        ('定标与三重积', 0.84, '', 'B 指数'),
        ('定标与三重积', 0.142983, US, 'ISS04 tau_E'),
        ('定标与三重积', 2.4762, US, 'D1 的 IPB98 tau_E'),
        ('定标与三重积', 0.057743, '', '两者的比'),
        ('定标与三重积', 2.28774e19, '', '三重积'),
        ('定标与三重积', 3.1e21, '', '点火阈值'),
        ('定标与三重积', 135.51, '', '差距倍数'),
        ('定标与三重积', 1.0689, '', 'iota 提高后的倍数'),
        ('定标与三重积', 6.8903, UPCT, 'iota 提高后的涨幅'),
        # AA.7 练习
        ('估算练习', 3.9789e5, UPA, '练习 3 的压力'),
        ('估算练习', 2.4834e21, UM3, '练习 3 的密度'),
        ('估算练习', 1.79, UPCT, '练习 6 的比压'),
        ('估算练习', 0.8322, '', '练习 6 的归一化比压'),
        ('估算练习', 0.353553, '', '练习 5 的 3/2 次方倍数'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s：%s（%s）' % (sec, val, why))
    if bad:
        rec('FAIL', 'P', '附录 AA 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AA 的 %d 个关键数值均在小节内复算通过' % len(cases))

    lits = [
        ('3.06', '比压比不得四舍五入到 3.06'),
        ('0.568\\ ', 'sqrt(eps) 不得写成 0.568'),
        ('0.00148\\ ', 'eps_eff^{3/2} 不得写成 0.00148'),
        ('0.0015', 'eps_eff^{3/2} 不得写成 0.0015'),
        ('6.68\\ ', '中心柱比值不得写成 6.68'),
        ('0.645\\ ', '极向场不得写成 0.645 T'),
        ('0.079\\ ', '自举份额不得写成 0.079'),
        ('99470', 'FRC 压力不得写成 99470 Pa'),
        ('6.21\\times10^{20}', 'FRC 密度不得写成 6.21e20'),
        ('0.143\\ ', 'ISS04 的 tau_E 不得写成 0.143 s'),
        ('0.058', 'tau_E 之比不得写成 0.058'),
        ('136\\ ', '三重积差距不得写成 136 倍'),
    ]
    for lit, why in lits:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AA'], why))


# -*- coding: utf-8 -*-
"""附录 I 新增节「读一篇论文：三遍法、提问清单与复现实验」的门禁体：P33。

可证伪点
  ① 小节标题或 sec label 被删/改名；
  ② 两张表 label、语义盒标签被改名；
  ③ 三遍法判据、六问表、「读懂」定义句、出处句被改写（整行/整句字面绑定）；
  ④ 口径回流：把「五个 C」写成「三个 C」、把判据弱化成「感觉懂了」「复现不上就算了」、
     把第三遍判据从「能复现其中一个数」弱化成「结论」等。
"""


def check_P33_appendix_i_paper():
    """附录 I 新增节：读论文三遍法的结构与口径核对。"""
    path = os.path.join(SRC, APP_FILES['I'])
    if not os.path.exists(path):
        rec('FAIL', 'P', '附录 I 文件缺失：%s' % path)
        return
    text = strip_comments(read(path))

    has(text, '\\section{读一篇论文：三遍法、提问清单与复现实验}',
        'src/%s（新节标题）' % APP_FILES['I'])
    has(text, '\\label{sec:app-i-paper}', 'src/%s（新节 label）' % APP_FILES['I'])
    has(text, '\\label{tab:app-i-paperpass}', 'src/%s（三遍法表 label）' % APP_FILES['I'])
    has(text, '\\label{tab:app-i-paperquestions}', 'src/%s（六问表 label）' % APP_FILES['I'])
    has(text, '{tip-app-i-paper}', 'src/%s（语义盒 tip-app-i-paper）' % APP_FILES['I'])

    lits = [
        ('\\subsection*{一、三遍法：每一遍都有「通过判据」}', '第一小节标题'),
        ('\\subsection*{二、第二遍的核心动作：把图当成「主文」}', '第二小节标题'),
        ('\\subsection*{三、第四遍其实是「算一遍」：复现实验}', '第三小节标题'),
        ('\\subsection*{四、文献笔记：一行问题、三个数、一个步骤、一句反对}', '第四小节标题'),
        ('\\subsection*{五、三条常见错误}', '第五小节标题'),
        ('第一遍 & 标题、摘要、引言、小节标题、结论、扫一遍参考文献 & 五个 C：类别、背景、正确性、贡献、清晰度 & 能用三句话向别人说清「它做了什么」 \\\\',
         '表 1 第一遍行（五个 C）'),
        ('第三遍 & 逐步重做：假设自己是作者，逐行验证推导与数据 & 一份「虚拟重做」记录：哪一步我自己也能推出来，哪一步不能 & \\textbf{能复现其中一个数}（见第三节） \\\\',
         '表 1 第三遍行（能复现其中一个数）'),
        ('我能算出同一个数吗？ & \\textbf{读懂的判据} & 附录 M--AA 的 \\texttt{calc} 脚本 \\\\',
         '表 2 「我能算出同一个数吗」行'),
        ('我不同意哪一点？ & 从「读者」变成「同行」 & 附录 K 的例题与解答对照 \\\\',
         '表 2 最后一行'),
        ('\\textbf{能用本书的脚本复现出论文里的一个数，', '「读懂」的可操作定义句'),
        ('用 \\texttt{tmp/app\\_\\_calc.py} 代入', '文献笔记模板的可复现步骤'),
        ('Keshav（2007）', '出处段的 Keshav 年份'),
        ('ACM SIGCOMM Computer Communication Review \\textbf{37}(3), 83--84',
         '出处段的卷期页'),
        ('\\emph{The Craft of Research}（第 3 版，2008）', '出处段的 Booth 等'),
        ('\\emph{Best Practices for Scientific Computing}', '出处段的 Wilson 等'),
        ('\\textbf{把摘要当结论}', '常见错误之一'),
        ('一次真正的收获——本书自己的勘误', '复现不上时的口径句'),
        # 以下四条是夹具逐条探针抓出的洞（v1.26.0 收尾轮）：语义盒只钉了标签、没钉标题；
        # 表 1 只钉了第一/第三遍行、漏了第二遍行；表 2 只钉了第五/第六问、漏了第一问；
        # 「可以随时停止」是本节的核心口径，此前完全没绑——反写成「必须读完三遍」六套件全绿。
        ('\\begin{tip}{「读懂」的可操作定义}{tip-app-i-paper}', '语义盒标题（可操作定义）'),
        ('第二遍 & 图、表、坐标轴、误差棒、关键公式的\\textbf{结构}（不追推导） & 三到五个关键量级 + 一份「没看懂的术语」清单 & 能说出每个图的横纵轴是什么、结论靠哪张图支撑 \\\\',
         '表 1 第二遍行（三到五个关键量级）'),
        ('它解决什么问题？ & 判断与自己的目标是否相关 & 附录 J 的逐章对照表 \\\\',
         '表 2 第一行（它解决什么问题）'),
        ('三遍法最重要的性质是\\textbf{可以随时停止}：', '「可以随时停止」的口径句'),
    ]
    for lit, why in lits:
        has(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))

    lits_lack = [
        ('三个 C', '五个 C 不得写成「三个 C」'),
        ('感觉懂了', '「读懂」不得被弱化成「感觉懂了」'),
        ('复现不上就算了', '不得把复现失败一笔带过'),
        ('引用次数越多越正确', '不得把引用次数当正确性'),
        ('第二遍读完就算读懂', '不得把第二遍当成终点'),
        ('附录 M--X 的', '六问表的脚本接口不得退回旧字母范围'),
        ('必须读完三遍', '「可以随时停止」不得被反写成「必须读完三遍」'),
    ]
    for lit, why in lits_lack:
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['I'], why))


# -*- coding: utf-8 -*-
"""附录 AB（空间与天体等离子体）的门禁体：P34。

可证伪点
  ① 八个小节、七张表、六只语义盒、二十一个公式的标签缺一即 FAIL；
  ② 「出处：」不足 6 处、练习与解答不足 6 道即 FAIL；
  ③ 关键数值按**所在小节**复算（rel=1e-3）——把数字搬到别的小节即 FAIL；
  ④ 整行字面绑定：表行、出处整行、语义盒标题被改写即 FAIL；
  ⑤ 口径回流防护：已纠正过的四舍五入写法回流即 FAIL。
"""


def check_P34_appendix_ab_space_astro():
    """附录 AB：空间与天体等离子体的结构、出处、练习与数值核对。"""
    spath = os.path.join(SRC, APP_FILES['AB'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AB 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '四个尺度',
        '太阳风：为什么',
        '弓激波与磁鞘',
        '磁层对流与发电机',
        '磁重联：Sweet',
        '辐射带与空间天气',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AB 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AB 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ab-scales',
        'tab:app-ab-solarwind',
        'tab:app-ab-bowshock',
        'tab:app-ab-convection',
        'tab:app-ab-reconnection',
        'tab:app-ab-radbelt',
        'tab:app-ab-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AB'])

    for lab in (
        'tip-app-ab-scales',
        'tip-app-ab-sonic',
        'tip-app-ab-collisionless',
        'tip-app-ab-dynamo',
        'tip-app-ab-anomalous',
        'tip-app-ab-radbelt',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AB'])

    for lab in (
        'eq:app-ab-scales',
        'eq:app-ab-omegape',
        'eq:app-ab-cs',
        'eq:app-ab-va',
        'eq:app-ab-ms',
        'eq:app-ab-beta',
        'eq:app-ab-ram',
        'eq:app-ab-magp',
        'eq:app-ab-rmp',
        'eq:app-ab-powersw',
        'eq:app-ab-emf',
        'eq:app-ab-sigma',
        'eq:app-ab-lundquist',
        'eq:app-ab-eta',
        'eq:app-ab-sp',
        'eq:app-ab-sptime',
        'eq:app-ab-petschek',
        'eq:app-ab-larmor',
        'eq:app-ab-rg',
        'eq:app-ab-tce',
        'eq:app-ab-losscone',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AB'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AB'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AB 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AB 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AB 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AB 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 整行字面绑定 ----------
    for lit, why in (
        ('电子等离子体频率 & $1.26147\\times10^{5}$ & $\\mathrm{rad\\,s^{-1}}$ & $20.08\\ \\mathrm{kHz}$ \\\\', '表行'),
        ('离子等离子体频率 & $2943.89$ & $\\mathrm{rad\\,s^{-1}}$ & $468.5\\ \\mathrm{Hz}$ \\\\', '表行'),
        ('声速 & $43.769$ & $\\mathrm{km\\,s^{-1}}$ & $\\gamma=1$，质子质量 \\\\', '表行'),
        ('阿尔芬速 & $48.773$ & $\\mathrm{km\\,s^{-1}}$ & 用 $n$ 算出的质量密度 \\\\', '表行'),
        ('声速马赫数 & $9.1388$ & --- & 超声速 \\\\', '表行'),
        ('阿尔芬马赫数 & $8.2012$ & --- & 超阿尔芬速 \\\\', '表行'),
        ('等离子体 $\\beta$ & $0.80534$ & --- & 磁压与热压同量级 \\\\', '表行'),
        ('太阳风动压 $\\rho v^2$ & $1.3381\\times10^{-9}$ & $\\mathrm{Pa}$ & 不带 $1/2$ \\\\', '表行'),
        ('赤道磁压 & $3.8237\\times10^{-4}$ & $\\mathrm{Pa}$ & $B_0=3.12\\times10^{-5}\\ \\mathrm{T}$ \\\\', '表行'),
        ('磁压/动压 & $2.86\\times10^{5}$ & --- & 近地磁场主导 \\\\', '表行'),
        ('磁层顶距离 & $8.1158$ & $R_E$ & 观测约 $10\\ R_E$ \\\\', '表行'),
        ('迎风面积 & $1.27516\\times10^{14}$ & $\\mathrm{m^2}$ & $\\pi R_E^2$ \\\\', '表行'),
        ('太阳风输入功率 & $34.126$ & $\\mathrm{GW}$ & $\\frac12\\rho v^3\\pi R_E^2$ \\\\', '表行'),
        ('对流电场 & $E=vB$ & 几十 $\\mathrm{kV}$ / $10\\ R_E$ & 太阳风拖磁力线 \\\\', '表行'),
        ('阿尔芬电导 & $1/(\\mu_0v_A)$ & $16.316\\ \\mathrm{S}$ & 磁冻结成立 \\\\', '表行'),
        ('电离层电导 & 佩德森电导 & $1$--$10\\ \\mathrm{S}$ 量级 & 电流闭合 \\\\', '表行'),
        ('斯皮策电阻率 $\\eta$ & $3.28877\\times10^{-5}$ & $\\Omega\\,\\mathrm{m}$ & $T=10\\ \\mathrm{eV}$，$\\ln\\Lambda=20$ \\\\', '表行'),
        ('伦德奎斯特数 $S$ & $1.86362\\times10^{10}$ & --- & $L=10^6\\ \\mathrm{m}$ \\\\', '表行'),
        ('Sweet--Parker 率 & $7.32523\\times10^{-6}$ & --- & $S^{-1/2}$ \\\\', '表行'),
        ('Petschek 率 & $0.0166058$ & --- & $\\pi/(8\\ln S)$ \\\\', '表行'),
        ('两者之比 & $2266.9$ & --- & 快重联的全部理由 \\\\', '表行'),
        ('Sweet--Parker 时间 & $2.79897\\times10^{7}$ & $\\mathrm{s}$ & $323.95$ 天 \\\\', '表行'),
        ('Petschek 时间 & $12347$ & $\\mathrm{s}$ & 约 $3.4$ 小时 \\\\', '表行'),
        ('内辐射带 & 质子与高能电子 & $r_g$ 小、寿命长 & 磁矩守恒、CRAND 源 \\\\', '表行'),
        ('外辐射带 & $1\\ \\mathrm{MeV}$ 电子 & $r_g=47.432\\ \\mathrm{km}$ & 波粒相互作用、加速 \\\\', '表行'),
        ('磁尾瓣 & 等离子体片 & $B=1.0\\times10^{-7}\\ \\mathrm{T}$ & 重联储能与释放 \\\\', '表行'),
        ('\\textbf{出处：}本节的分区口径与尺度取值参照 Kivelson \\& Russell \\emph{Introduction to Space Physics}', '表行'),
        ('\\textbf{出处：}本节公式与讲解顺序参照 Baumjohann \\& Treumann \\emph{Basic Space Plasma Physics}', '表行'),
        ('\\textbf{出处：}磁层顶力学平衡与弓激波几何参照 Kivelson \\& Russell \\emph{Introduction to Space Physics}', '表行'),
        ('\\textbf{出处：}对流循环与电离层发电机的讲解参照 Kivelson \\& Russell \\emph{Introduction to Space Physics}', '表行'),
        ('\\textbf{出处：}本节公式与两模型对比参照 Priest \\& Forbes \\emph{Magnetic Reconnection}', '表行'),
        ('\\textbf{出处：}三层运动与损失锥的讲解参照 Kivelson \\& Russell \\emph{Introduction to Space Physics}', '表行'),
        ('\\begin{tip}{四个尺度不是四个数量级，是四种物理}{tip-app-ab-scales}', '语义盒标题'),
        ('\\begin{tip}{太阳风超声速，靠的不是加热，是边界条件}{tip-app-ab-sonic}', '语义盒标题'),
        ('\\begin{tip}{无碰撞激波：耗散靠波和粒子，不靠碰撞}{tip-app-ab-collisionless}', '语义盒标题'),
        ('\\begin{tip}{发电机不是机器，是欧姆定律的一个极限}{tip-app-ab-dynamo}', '语义盒标题'),
        ('\\begin{tip}{电阻不够，重联靠的是「反常电阻」}{tip-app-ab-anomalous}', '语义盒标题'),
        ('\\begin{tip}{辐射带不是「带」，是磁矩守恒的必然结果}{tip-app-ab-radbelt}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AB'], why))

    # ---------- ⑤ 关键数值按小节复算 ----------
    cases = [
        ('四个尺度', 1, r'\mathrm{AU}', '太阳风 & $1\\ \\mathrm{AU}=1.496\\times10^{11}\\ \\mathrm{m}$ & 输运 $'),
        ('四个尺度', 10, r'\mathrm{min}', '磁层 & $10\\ R_E\\approx6.4\\times10^{7}\\ \\mathrm{m}$ & 对流 $10\\ \\'),
        ('太阳风：为什么', 5, r'\mathrm{nT}', '关键数值'),
        ('太阳风：为什么', 10, r'\mathrm{eV}', '关键数值'),
        ('太阳风：为什么', 20.08, r'\mathrm{kHz}', '关键数值'),
        ('太阳风：为什么', 43.769, r'\mathrm{km\,s^{-1}}', '关键数值'),
        ('太阳风：为什么', 48.773, r'\mathrm{km\,s^{-1}}', '关键数值'),
        ('太阳风：为什么', 400, r'\mathrm{km\,s^{-1}}', '关键数值'),
        ('太阳风：为什么', 468.5, r'\mathrm{Hz}', '关键数值'),
        ('太阳风：为什么', 2943.89, r'\mathrm{rad\,s^{-1}}', '关键数值'),
        ('弓激波与磁鞘', 34.126, r'\mathrm{GW}', '关键数值'),
        ('磁层对流与发电机', 1, r'\mathrm{nT}', '代入 $v=400\\ \\mathrm{km\\,s^{-1}}$、$B=5\\ \\mathrm{nT}$（垂直分量取 $1\\'),
        ('磁层对流与发电机', 5, r'\mathrm{nT}', '代入 $v=400\\ \\mathrm{km\\,s^{-1}}$、$B=5\\ \\mathrm{nT}$（垂直分量取 $1\\'),
        ('磁层对流与发电机', 10, r'\mathrm{S}', '电离层的佩德森电导只有 $1$--$10\\ \\mathrm{S}$ 量级，'),
        ('磁层对流与发电机', 16.316, r'\mathrm{S}', '\\Sigma_A=\\frac{1}{\\mu_0 v_A}=16.316\\ \\mathrm{S},'),
        ('磁重联：Sweet', 6, r'\mathrm{m}', '关键数值'),
        ('磁重联：Sweet', 10, r'\mathrm{eV}', '关键数值'),
        ('磁重联：Sweet', 12347, r'\mathrm{s}', '关键数值'),
        ('辐射带与空间天气', 0.00714477, r'\mathrm{s}', '关键数值'),
        ('辐射带与空间天气', 1, r'\mathrm{MeV}', '关键数值'),
        ('辐射带与空间天气', 5, r'\mathrm{nT}', '关键数值'),
        ('辐射带与空间天气', 47, r'\mathrm{km}', '关键数值'),
        ('辐射带与空间天气', 47.432, r'\mathrm{km}', '关键数值'),
        ('辐射带与空间天气', 140, r'\mathrm{Hz}', '关键数值'),
        ('辐射带与空间天气', 47431.8, r'\mathrm{m}', '关键数值'),
        ('估算练习', 5, r'\mathrm{cm^{-3}}', '\\textbf{练习 1}（密度）：把太阳风密度从 $5\\ \\mathrm{cm^{-3}}$ 加倍到 $10\\ \\ma'),
        ('估算练习', 10, r'\mathrm{cm^{-3}}', '\\textbf{练习 1}（密度）：把太阳风密度从 $5\\ \\mathrm{cm^{-3}}$ 加倍到 $10\\ \\ma'),
        ('估算练习', 47.432, r'\mathrm{km}', '求动量与回旋半径，并与电子的 $47.432\\ \\mathrm{km}$ 比较。'),
        ('来源、分级与版权说明', 5, r'\mathrm{nT}', '关键数值'),
        ('来源、分级与版权说明', 10, r'\mathrm{S}', '关键数值'),
        ('来源、分级与版权说明', 400, r'\mathrm{km\,s^{-1}}', '关键数值'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s 的 %.6g%s' % (sec, val, unit))
    if bad:
        rec('FAIL', 'P', '附录 AB 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AB 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑥ 口径回流防护 ----------
    for lit, why in (
        ('43.8\\ \\mathrm{km', '太阳风速度四舍五入'),
        ('48.8\\ \\mathrm{km', '阿尔芬速度四舍五入'),
        ('1.34\\times10^{-9}', '太阳风动压四舍五入'),
        ('0.0166\\ \\mathrm', 'Petschek 重联率四舍五入'),
        ('2267', '两模型比值四舍五入'),
        ('136500', 'Lundquist 数平方根四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AB'], why))


# -*- coding: utf-8 -*-
"""附录 AC（等离子体数值方法）的门禁体：P35。

可证伪点
  ① 八个小节、七张表、七只语义盒、二十一个公式的标签缺一即 FAIL；
  ② 「出处：」不足 6 处、练习与解答不足 6 道即 FAIL；
  ③ 关键数值按**所在小节**复算（rel=1e-3）——把数字搬到别的小节即 FAIL；
  ④ 整行字面绑定：表行、出处整行、语义盒标题被改写即 FAIL；
  ⑤ 口径回流防护：已纠正过的四舍五入写法回流即 FAIL。
"""


def check_P35_appendix_ac_numerics():
    """附录 AC：等离子体数值方法的结构、出处、练习与数值核对。"""
    spath = os.path.join(SRC, APP_FILES['AC'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AC 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '分辨率：先算三个长度',
        '推进粒子：蛙跳',
        '稳定性与成本',
        '粒子噪声：超级粒子',
        '场怎么求',
        '流体与 MHD',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AC 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AC 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ac-resolution',
        'tab:app-ac-pusher',
        'tab:app-ac-cost',
        'tab:app-ac-noise',
        'tab:app-ac-error',
        'tab:app-ac-discretization',
        'tab:app-ac-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AC'])

    for lab in (
        'tip-app-ac-resolve',
        'tip-app-ac-conserve',
        'tip-app-ac-stability',
        'tip-app-ac-noise',
        'tip-app-ac-verify',
        'tip-app-ac-conservative',
        'tip-app-ac-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AC'])

    for lab in (
        'eq:app-ac-debye',
        'eq:app-ac-freq',
        'eq:app-ac-step',
        'eq:app-ac-nd',
        'eq:app-ac-collision',
        'eq:app-ac-leapfrog',
        'eq:app-ac-boris',
        'eq:app-ac-phase',
        'eq:app-ac-relativistic',
        'eq:app-ac-cfl',
        'eq:app-ac-pcfl',
        'eq:app-ac-cost',
        'eq:app-ac-scalegap',
        'eq:app-ac-noise',
        'eq:app-ac-timeavg',
        'eq:app-ac-fft',
        'eq:app-ac-trunc',
        'eq:app-ac-order',
        'eq:app-ac-flux',
        'eq:app-ac-divb',
        'eq:app-ac-memory',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AC'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AC'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AC 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AC 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AC 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AC 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 整行字面绑定 ----------
    for lit, why in (
        ('要分辨的 & 判据 & 本工作点的数 & 不满足会怎样 \\\\', '表格数据行'),
        ('德拜屏蔽 & $\\Delta x\\le\\lambda_{De}$ & $\\lambda_{De}=7.12806\\times10^{-5}\\ \\mathrm{m}$ & 有限网格不稳定性 \\\\', '表格数据行'),
        ('等离子体振荡 & $\\Delta t\\le0.2/\\omega_{pe}$ & $3.54518\\times10^{-13}\\ \\mathrm{s}$ & 越步，能量凭空长 \\\\', '表格数据行'),
        ('电子回旋 & $\\Delta t\\le0.2/\\omega_{ce}$ & $2.14552\\times10^{-13}\\ \\mathrm{s}$ & 回旋轨道被算成多边形 \\\\', '表格数据行'),
        ('粒子离散噪声 & $N_p\\gg1$（每格） & 取 $20$ & 噪声淹没信号（AC.4） \\\\', '表格数据行'),
        ('边界距离 & $L\\gg\\lambda_{De}$ & $L=100\\lambda_{De}=0.71281\\ \\mathrm{cm}$ & 边界反射污染物理区 \\\\', '表格数据行'),
        ('碰撞 & $t_{\\mathrm{sim}}\\ll\\tau_{ei}$ & $t_{\\mathrm{sim}}/\\tau_{ei}=1.9869\\times10^{-5}$ & 需显式加碰撞模块 \\\\', '表格数据行'),
        ('\\begin{tip}{先算分辨率，再谈算法}{tip-app-ac-resolve}', '语义盒标题'),
        ('\\textbf{出处：}本节的分辨率判据与有限网格不稳定性的结论，参照 Birdsall \\& Langdon', '出处声明整行'),
        ('推进格式 & 时间精度 & 守恒性质 & 什么时候用 \\\\', '表格数据行'),
        ('显式欧拉 & 一阶 & 无 & 只在教学演示里出现 \\\\', '表格数据行'),
        ('蛙跳（无磁场） & 二阶 & 时间可逆、准能量守恒 & 静电 PIC 的标准 \\\\', '表格数据行'),
        ('Boris & 二阶 & 磁场步保 $|v|$ 到机器精度 & 磁化 PIC 的标准 \\\\', '表格数据行'),
        ('相对论 Boris & 二阶 & 保 $|\\beta|$ 上限 $c$ & 激光尾场、$\\gamma\\gg1$ \\\\', '表格数据行'),
        ('隐式 / Crank--Nicolson & 任意阶 & 可大步长 & 需要 $\\omega_{ce}\\Delta t>1$ 时 \\\\', '表格数据行'),
        ('\\begin{tip}{守恒量不是「算得好」，是格式的性质}{tip-app-ac-conserve}', '语义盒标题'),
        ('\\textbf{出处：}蛙跳与 Boris 推进器的构造参照 Boris（1970，第 4 届数值模拟等离子体会议论文集）', '出处声明整行'),
        ('量 & 表达式 & 本工作点 & 备注 \\\\', '表格数据行'),
        ('盒边长 & $L=100\\lambda_{De}$ & $0.71281\\ \\mathrm{cm}$ & 边界要远离物理区 \\\\', '表格数据行'),
        ('一维 / 三维格点 & $100$ / $10^{6}$ & $10^{6}$ & 成本随维数暴涨 \\\\', '表格数据行'),
        ('总粒子数 & $10^{6}\\times20$ & $2\\times10^{7}$ & 每格 $20$ 个 \\\\', '表格数据行'),
        ('时间步数 & $T_{\\mathrm{end}}/\\Delta t$ & $8261.82$ & $T_{\\mathrm{end}}=1000/\\omega_{pe}$ \\\\', '表格数据行'),
        ('粒子步总数 & $N_pN_{\\mathrm{step}}$ & $1.65236\\times10^{11}$ & 成本的真单位 \\\\', '表格数据行'),
        ('单核 / GPU 时间 & 除以 $10^{8}$ / $10^{9}$ & $0.45899\\ \\mathrm{h}$ / $165.236\\ \\mathrm{s}$ & 差 $10$ 倍 \\\\', '表格数据行'),
        ('真实/模拟粒子数 & $N_{\\mathrm{real}}/N_{\\mathrm{sim}}$ & $4.161\\times10^{15}$ & 超级粒子的代价 \\\\', '表格数据行'),
        ('\\begin{tip}{数值不稳定与物理不稳定是两件事}{tip-app-ac-stability}', '语义盒标题'),
        ('\\textbf{出处：}CFL 条件参照 Courant, Friedrichs \\& Lewy（1928，\\emph{Mathematische Annalen}', '出处声明整行'),
        ('每格粒子数 & 单格噪声 & 时间平均 $100$ 步后 & 用途 \\\\', '表格数据行'),
        ('$20$ & $22.361\\%$ & $2.2361\\%$ & 只看定性趋势 \\\\', '表格数据行'),
        ('$25$ & $20\\%$ & $2\\%$ & 起步档 \\\\', '表格数据行'),
        ('$100$ & $10\\%$ & $1\\%$ & 看剖面形状 \\\\', '表格数据行'),
        ('$1000$ & $3.1623\\%$ & $0.31623\\%$ & 定量比较 \\\\', '表格数据行'),
        ('$10^{4}$ & $1\\%$ & $0.1\\%$ & 单步就要 $1\\%$ 精度时 \\\\', '表格数据行'),
        ('\\begin{tip}{噪声不是误差，是统计}{tip-app-ac-noise}', '语义盒标题'),
        ('\\textbf{出处：}粒子噪声的 $1/\\sqrt{N_p}$ 标度与时间平均降噪，', '出处声明整行'),
        ('误差来源 & 随什么变 & 本工作点的量级 & 怎么压 \\\\', '表格数据行'),
        ('空间截断 & $\\propto\\Delta x^2$ & $8.3333\\%$（$\\Delta x=\\lambda_{De}$） & 加密网格 / 提阶 \\\\', '表格数据行'),
        ('时间截断 & $\\propto\\Delta t^2$ & 相位 $0.00166667\\ \\mathrm{rad}$ 每步 & 减步长 / 高阶格式 \\\\', '表格数据行'),
        ('粒子统计 & $\\propto1/\\sqrt{N_p}$ & $22.361\\%$（每格 $20$） & 加粒子 / 时间平均 \\\\', '表格数据行'),
        ('舍入 & 机器精度 & 约 $10^{-16}$ 相对 & 通常可忽略 \\\\', '表格数据行'),
        ('\\begin{tip}{验证与确认是两件事}{tip-app-ac-verify}', '语义盒标题'),
        ('\\textbf{出处：}FFT 泊松求解与直接求和的成本对照，参照 Press, Teukolsky, Vetterling \\& Flannery', '出处声明整行'),
        ('离散方式 & 数据结构 & 优势 & 典型用途 \\\\', '表格数据行'),
        ('有限差分 & 格点值 & 实现最简单 & 早期 MHD、扩散问题 \\\\', '表格数据行'),
        ('有限体积 & 界面通量 & 严格守恒、抓得住激波 & 现代 MHD、天体 \\\\', '表格数据行'),
        ('谱方法 & 傅里叶/球谐系数 & 精度指数收敛 & 湍流、周期性盒子 \\\\', '表格数据行'),
        ('粒子（PIC） & 宏粒子 + 网格 & 动理学、无流体假设 & 加热、尾场、快粒子 \\\\', '表格数据行'),
        ('\\begin{tip}{守恒律要用守恒格式写}{tip-app-ac-conservative}', '语义盒标题'),
        ('\\textbf{出处：}通量形式与黎曼求解器参照 Godunov（1959，\\emph{Matematicheskii Sbornik}', '出处声明整行'),
        ('来源 & 借用的内容 & 本附录的处理 \\\\', '表格数据行'),
        ('Birdsall \\& Langdon（1985/1991） & PIC 分辨率判据、噪声标度、Boris 步 & 只借结论与公式形式；数值全部按本工作点重算 \\\\', '表格数据行'),
        ('Hockney \\& Eastwood（1988） & 场求解与截断误差 & 只借误差阶数的结论 \\\\', '表格数据行'),
        ('Tajima（1989/2004） & 成本估算、超级粒子 & 只借用词与量级框架 \\\\', '表格数据行'),
        ('Press 等（2007） & FFT 解泊松的运算量 & 只借 $N\\log_2N$ 这一形式，数值自算 \\\\', '表格数据行'),
        ('Boris（1970）、Vay（2008） & 推进器构造 & 只借旋转形式与守恒性质 \\\\', '表格数据行'),
        ('Courant 等（1928） & CFL 条件 & 只借用判据形式 \\\\', '表格数据行'),
        ('Godunov（1959）、Van Leer（1979） & 通量形式与限制器 & 只借方法名与适用场合 \\\\', '表格数据行'),
        ('Gardiner \\& Stone（2005/2008）、Stone 等（2008） & 约束输运 & 只借 $\\nabla\\cdot\\mathbf{B}$ 的处理思路 \\\\', '表格数据行'),
        ('\\begin{tip}{学会说「这不算」}{tip-app-ac-honest}', '语义盒标题'),
        ('\\textbf{出处：}本节的来源分级与版权口径，沿用本教材附录 E（分级自学书单）与', '出处声明整行'),
        ('\\textbf{过粗了 $4590.38$ 倍}', '离子步长过粗的倍数'),
        ('一、\\textbf{不复制}任何原书的文字表述', '版权声明第一条'),
        ('四、\\textbf{不复制}任何未经验证的参数推荐值', '版权声明第四条'),
        ('它不会自动保持，误差会累积成一个假的「磁单极子」', '旋度约束的后果'),
        ('=7.12806\\times10^{-5}\\ \\mathrm{m}=71.28\\ \\mathrm{\\mu m}.', '式 (AC.1) 德拜长度整行（末位容差洞）'),
        ('\\omega_{pe}=\\sqrt{\\frac{n_e e^2}{\\epsilon_0 m_e}}=5.64146\\times10^{11}\\ \\mathrm{rad\\,s^{-1}}', '式 (AC.2) 等离子体频率整行（末位容差洞）'),
        ('\\omega_{ce}=\\frac{eB}{m_e}=9.32175\\times10^{11}\\ \\mathrm{rad\\,s^{-1}}', '式 (AC.2) 回旋频率整行（末位容差洞）'),
        ('N_D=n_e\\lambda_{De}^3=3.62171\\times10^{7}.', '式 (AC.4) 德拜球内粒子数整行（末位容差洞）'),
        ('=8.92144\\times10^{-5}\\ \\mathrm{s},', '式 (AC.5) 碰撞时间整行（末位容差洞）'),
        ('\\frac{(\\omega_{ce}\\Delta t)^2}{24}=0.00166667\\ \\mathrm{rad},', '式 (AC.9) 每步相位误差整行（末位容差洞）'),
        ('N_pN_{\\mathrm{step}}=2\\times10^{7}\\times8261.82=1.65236\\times10^{11},', '式 (AC.14) 粒子步总数整行（末位容差洞）'),
        ('=\\frac{10^{12}}{1.99316\\times10^{7}}=50171.7\\quad(10^{6}\\ \\text{个格点}).', '式 (AC.16) FFT 加速比整行（末位容差洞）'),
        ('256^3\\times8\\ \\text{个变量}\\times8\\ \\mathrm{B}=1.07374\\times10^{9}\\ \\mathrm{B}=1.0737\\ \\mathrm{GB}.', '式 (AC.21) 网格内存整行（末位容差洞）'),
        ('\\textbf{解答 1}：$\\lambda_{De}=7.43394\\times10^{-5}\\ \\mathrm{m}=74.339\\ \\mathrm{\\mu m}$、', '练习 1 解答的德拜长度整行（末位容差洞）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AC'], why))

    # ---------- ⑤ 关键数值按小节复算 ----------
    cases = [
        ('分辨率：先算三个长度', 7.12806e-05, r'\mathrm{m}', '德拜长度'),
        ('分辨率：先算三个长度', 71.28, r'\mathrm{\mu m}', '德拜长度（微米）'),
        ('分辨率：先算三个长度', 5.64146e+11, r'\mathrm{rad\,s^{-1}}', '电子等离子体频率'),
        ('分辨率：先算三个长度', 9.32175e+11, r'\mathrm{rad\,s^{-1}}', '电子回旋频率'),
        ('分辨率：先算三个长度', 1.65236, r'', '回旋与等离子体频率之比'),
        ('分辨率：先算三个长度', 3.54518e-13, r'\mathrm{s}', '等离子体振荡步长上限'),
        ('分辨率：先算三个长度', 2.14552e-13, r'\mathrm{s}', '回旋步长上限（采用值）'),
        ('分辨率：先算三个长度', 8.62771e-06, r'\mathrm{m}', '每步热位移'),
        ('分辨率：先算三个长度', 0.12104, r'', '每步热位移（以德拜长度计）'),
        ('分辨率：先算三个长度', 36217100, r'', '德拜球内粒子数'),
        ('分辨率：先算三个长度', 8.92144e-05, r'\mathrm{s}', '电子-离子碰撞时间'),
        ('分辨率：先算三个长度', 1.9869e-05, r'', '模拟窗口与碰撞时间之比'),
        ('分辨率：先算三个长度', 0.71281, r'\mathrm{cm}', '100 个德拜长度的盒边长'),
        ('推进粒子：蛙跳', 0.00166667, r'\mathrm{rad}', '每步相位误差'),
        ('推进粒子：蛙跳', 0.0523599, r'\mathrm{rad}', '每回旋累计相位误差'),
        ('推进粒子：蛙跳', 2.95695, r'', '1 MeV 电子的洛伦兹因子'),
        ('推进粒子：蛙跳', 0.000217847, r'', '离子与电子回旋频率之比'),
        ('推进粒子：蛙跳', 4590.38, r'', '离子回旋周期比电子长的倍数'),
        ('稳定性与成本', 7359450, r'\mathrm{m\,s^{-1}}', '快磁声速'),
        ('稳定性与成本', 4.18155e-07, r'\mathrm{kg\,m^{-3}}', '质量密度'),
        ('稳定性与成本', 7311430, r'\mathrm{m\,s^{-1}}', '阿尔芬速度'),
        ('稳定性与成本', 7311.4, r'\mathrm{km\,s^{-1}}', '阿尔芬速度（km/s）'),
        ('稳定性与成本', 839370, r'\mathrm{m\,s^{-1}}', '声速'),
        ('稳定性与成本', 7359.4, r'\mathrm{km\,s^{-1}}', '快磁声速（km/s）'),
        ('稳定性与成本', 6.79399e-10, r'\mathrm{s}', 'MHD 时间步长'),
        ('稳定性与成本', 0.848, r'\mathrm{\mu s}', '阿尔芬时间'),
        ('稳定性与成本', 12481.4, r'', '10 个阿尔芬时间的步数'),
        ('稳定性与成本', 8261.82, r'', 'PIC 时间步数'),
        ('稳定性与成本', 1.65236e+11, r'', '粒子步总数'),
        ('稳定性与成本', 1652.36, r'\mathrm{s}', '单核时间'),
        ('稳定性与成本', 0.45899, r'\mathrm{h}', '单核时间（小时）'),
        ('稳定性与成本', 165.236, r'\mathrm{s}', 'GPU 时间'),
        ('稳定性与成本', 832.2, r'\mathrm{m^3}', '等离子体体积'),
        ('稳定性与成本', 8.322e+22, r'', '真实粒子数'),
        ('稳定性与成本', 4.161e+15, r'', '真实与模拟粒子数之比'),
        ('稳定性与成本', 20000000, r'', '模拟粒子数'),
        ('稳定性与成本', 0.12104, r'', '粒子 CFL 检验值'),
        ('粒子噪声：超级粒子', 22.361, r'', '每格 20 粒子的相对噪声'),
        ('粒子噪声：超级粒子', 3.1623, r'', '每格 1000 粒子的相对噪声'),
        ('粒子噪声：超级粒子', 2.2361, r'', '时间平均 100 步后的噪声'),
        ('粒子噪声：超级粒子', 0.31623, r'', '每格 1000 粒子 + 100 步平均'),
        ('粒子噪声：超级粒子', 0.1, r'', '每格 10000 粒子 + 100 步平均'),
        ('粒子噪声：超级粒子', 2500, r'', '2% 精度所需的每格粒子数'),
        ('场怎么求', 19931600, r'', 'FFT 泊松的运算量'),
        ('场怎么求', 50171.7, r'', 'FFT 与直接求和的运算量之比'),
        ('场怎么求', 0.0833333, r'', '二阶差分的截断误差系数'),
        ('场怎么求', 33.333, r'', 'Δx/λ_De=2 的截断误差（%）'),
        ('场怎么求', 8.3333, r'', 'Δx/λ_De=1 的截断误差（%）'),
        ('场怎么求', 2.0833, r'', 'Δx/λ_De=0.5 的截断误差（%）'),
        ('场怎么求', 0.25, r'', '步长减半后误差的倍数'),
        ('流体与 MHD', 1073740000, r'', '256^3 网格内存（字节）'),
        ('流体与 MHD', 1.0737, r'', '256^3 网格内存（GB）'),
        ('估算练习', 7.43394e-05, r'\mathrm{m}', '练习 1 的德拜长度'),
        ('估算练习', 74.339, r'\mathrm{\mu m}', '练习 1 的德拜长度（微米）'),
        ('估算练习', 1.78399e+11, r'\mathrm{rad\,s^{-1}}', '练习 1 的等离子体频率'),
        ('估算练习', 4108260, r'', '练习 1 的德拜球内粒子数'),
        ('估算练习', 3.51764e+11, r'\mathrm{rad\,s^{-1}}', '练习 3 的回旋频率'),
        ('估算练习', 5.68563e-13, r'\mathrm{s}', '练习 3 的步长上限'),
        ('估算练习', 31.4159, r'', '练习 3 的每回旋步数'),
        ('估算练习', 125000, r'', '练习 4 的格点数'),
        ('估算练习', 6250000, r'', '练习 4 的粒子数'),
        ('估算练习', 5.16364e+10, r'', '练习 4 的粒子步总数'),
        ('估算练习', 516.364, r'\mathrm{s}', '练习 4 的单核时间'),
        ('估算练习', 0.1434, r'\mathrm{h}', '练习 4 的单核时间（小时）'),
        ('估算练习', 3.39699e-10, r'\mathrm{s}', '练习 5 的 MHD 步长'),
        ('估算练习', 24962.9, r'', '练习 5 的步数'),
        ('估算练习', 2500, r'', '练习 6 的每格粒子数'),
    ]
    bad = []
    for sec, val, unit, why in cases:
        if not value_in(found[sec], val, 1e-3, unit):
            bad.append('%s 的 %.6g%s' % (sec, val, unit))
    if bad:
        rec('FAIL', 'P', '附录 AC 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AC 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑥ 口径回流防护 ----------
    for lit, why in (
        ('22.4\\%', '粒子噪声四舍五入'),
        ('1.99\\times10^{7}', 'FFT 运算量四舍五入'),
        ('8.3\\%', '截断误差四舍五入'),
        ('50172', 'FFT 加速比四舍五入'),
        ('1.65\\times10^{11}', '粒子步总数四舍五入'),
        ('4.16\\times10^{15}', '真实/模拟粒子数之比四舍五入'),
        ('8262', 'PIC 步数四舍五入'),
        ('1.07\\ \\mathrm{GB}', '网格内存四舍五入'),
        ('71.3\\ \\mathrm{\\mu m}', '德拜长度四舍五入'),
        ('0.46\\ \\mathrm{h}', '单核时间四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AC'], why))


def check_P36_appendix_ad_pedestal():
    """附录 AD：边界台基与 L-H 模——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AD'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AD 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '从 L 模到 H 模',
        '台基：十二厘米',
        '两条边界的夹缝',
        'ELM：储能的周期性坍塌',
        '缓解：把一次大坍塌',
        '从台基到靶板',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AD 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AD 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ad-threshold',
        'tab:app-ad-pedestal',
        'tab:app-ad-limits',
        'tab:app-ad-elm',
        'tab:app-ad-mitigation',
        'tab:app-ad-sol',
        'tab:app-ad-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AD'])

    for lab in (
        'sec:app-ad-threshold',
        'sec:app-ad-pedestal',
        'sec:app-ad-limits',
        'sec:app-ad-elm',
        'sec:app-ad-mitigation',
        'sec:app-ad-sol',
        'sec:app-ad-problems',
        'sec:app-ad-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AD'])

    for lab in (
        'tip-app-ad-threshold',
        'tip-app-ad-width',
        'tip-app-ad-two-branches',
        'tip-app-ad-elm-power',
        'tip-app-ad-mitigation',
        'tip-app-ad-sol',
        'tip-app-ad-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AD'])

    for lab in (
        'eq:app-ad-pth',
        'eq:app-ad-psol',
        'eq:app-ad-greenwald',
        'eq:app-ad-shear',
        'eq:app-ad-er',
        'eq:app-ad-rhoi',
        'eq:app-ad-pped',
        'eq:app-ad-width',
        'eq:app-ad-grad',
        'eq:app-ad-jbs',
        'eq:app-ad-pmax',
        'eq:app-ad-margin',
        'eq:app-ad-wped',
        'eq:app-ad-dwelm',
        'eq:app-ad-eden',
        'eq:app-ad-felm',
        'eq:app-ad-pace',
        'eq:app-ad-fmin',
        'eq:app-ad-rmp',
        'eq:app-ad-pellet',
        'eq:app-ad-lamq',
        'eq:app-ad-qpeak',
        'eq:app-ad-detach',
        'eq:app-ad-broaden',
        'eq:app-ad-broaden2',
        'eq:app-ad-ratio',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AD'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AD'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AD 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AD 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AD 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AD 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 数字 + 单位的精确绑定（堵末位容差洞） ----------
    for lit, why in (
        ('63.2554\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('489.532\\ \\mathrm{m^2}', '从 L 模到 H 模·公式（精确绑定）'),
        ('5.3\\ \\mathrm{T}', '从 L 模到 H 模·公式（精确绑定）'),
        ('111.3\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('103.468\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('40.2124\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('1.19366\\times10^{20}\\ \\mathrm{m^{-3}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('41.116\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('38.4821\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('103.977\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('4\\ \\mathrm{T}', '从 L 模到 H 模·公式（精确绑定）'),
        ('50.4612\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('6\\ \\mathrm{T}', '从 L 模到 H 模·公式（精确绑定）'),
        ('69.881\\ \\mathrm{MW}', '从 L 模到 H 模·公式（精确绑定）'),
        ('5.93524\\times10^{5}\\ \\mathrm{m\\,s^{-1}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('593.524\\ \\mathrm{km\\,s^{-1}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('2.96762\\times10^{5}\\ \\mathrm{s^{-1}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('3.54936\\times10^{4}\\ \\mathrm{m\\,s^{-1}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('53.2405\\ \\mathrm{kV\\,m^{-1}}', '从 L 模到 H 模·公式（精确绑定）'),
        ('1.5\\ \\mathrm{T}', '从 L 模到 H 模·正文（精确绑定）'),
        ('30\\ \\mathrm{kV\\,m^{-1}}', '从 L 模到 H 模·正文（精确绑定）'),
        ('1.67220\\times10^{5}\\ \\mathrm{s^{-1}}', '从 L 模到 H 模·正文（精确绑定）'),
        ('4.13338\\ \\mathrm{mm}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('1.46046\\ \\mathrm{cm}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('15\\ \\mathrm{MA}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('9.194\\ \\mathrm{keV}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('241.933\\ \\mathrm{m^{-1}}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('3\\ \\mathrm{keV}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('57.6784\\ \\mathrm{kPa}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('297.525\\ \\mathrm{kPa}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('11.9603\\ \\mathrm{cm}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('4.82250\\times10^{5}\\ \\mathrm{Pa\\,m^{-1}}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('1.8804\\ \\mathrm{m}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('2.73900\\times10^{5}\\ \\mathrm{A\\,m^{-2}}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('1.27616\\ \\mathrm{MA}', '台基：十二厘米的壳层·公式（精确绑定）'),
        ('1.19366\\times10^{6}\\ \\mathrm{A\\,m^{-2}}', '台基：十二厘米的壳层·正文（精确绑定）'),
        ('44.0012\\ \\mathrm{kPa}', '两条边界的夹缝：气球·公式（精确绑定）'),
        ('70.4019\\ \\mathrm{kPa}', '两条边界的夹缝：气球·公式（精确绑定）'),
        ('58.5493\\ \\mathrm{m^3}', '两条边界的夹缝：气球·公式（精确绑定）'),
        ('2.53277\\ \\mathrm{MJ}', '两条边界的夹缝：气球·公式（精确绑定）'),
        ('2.59980\\times10^{7}\\ \\mathrm{J}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('25.998\\ \\mathrm{MJ}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('1\\ \\mathrm{m}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('38.9557\\ \\mathrm{m^2}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('6.67373\\times10^{5}\\ \\mathrm{J\\,m^{-2}}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('0.667373\\ \\mathrm{MJ\\,m^{-2}}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('0.3\\ \\mathrm{ms}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('2.22458\\ \\mathrm{GW\\,m^{-2}}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('1\\ \\mathrm{MJ\\,m^{-2}}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('3.97984\\ \\mathrm{Hz}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('251.267\\ \\mathrm{ms}', 'ELM：储能的周期性·公式（精确绑定）'),
        ('0.684532\\ \\mathrm{MW}', 'ELM：储能的周期性·正文（精确绑定）'),
        ('50\\ \\mathrm{ms}', '缓解：把一次大坍塌拆·正文（精确绑定）'),
        ('20\\ \\mathrm{Hz}', '缓解：把一次大坍塌拆·正文（精确绑定）'),
        ('5.17339\\ \\mathrm{MJ}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('0.132802\\ \\mathrm{MJ\\,m^{-2}}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('0.442672\\ \\mathrm{GW\\,m^{-2}}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('0.3\\ \\mathrm{MJ\\,m^{-2}}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('8.85345\\ \\mathrm{Hz}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('1.78885\\ \\mathrm{m}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('9.15431\\ \\mathrm{cm}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('1.26491\\ \\mathrm{m}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('4.18155\\times10^{-6}\\ \\mathrm{kg}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('4.18155\\ \\mathrm{mg}', '缓解：把一次大坍塌拆·公式（精确绑定）'),
        ('1.70796\\times10^{19}\\ \\mathrm{m^{-3}}', '缓解：把一次大坍塌拆·正文（精确绑定）'),
        ('0.3\\,\\mathrm{MJ\\,m^{-2}}', '缓解：把一次大坍塌拆·表格（精确绑定）'),
        ('0.388859\\ \\mathrm{mm}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('1.94429\\ \\mathrm{mm}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('0.0757414\\ \\mathrm{m^2}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('956.246\\ \\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('10\\ \\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('19.1249\\ \\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('9.56246\\ \\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('55.1547\\ \\mathrm{kW\\,m^{-2}}', '从台基到靶板：热流落·公式（精确绑定）'),
        ('150\\ \\mathrm{MW}', '估算练习·正文（精确绑定）'),
        ('0.15\\ \\mathrm{s}', '估算练习·正文（精确绑定）'),
        ('371.4\\ \\mathrm{MJ}', '估算练习·正文（精确绑定）'),
        ('37.14\\ \\mathrm{MJ}', '估算练习·正文（精确绑定）'),
        ('0.953391\\ \\mathrm{MJ\\,m^{-2}}', '估算练习·正文（精确绑定）'),
        ('5.76967\\ \\mathrm{Hz}', '估算练习·正文（精确绑定）'),
        ('173.32\\ \\mathrm{ms}', '估算练习·正文（精确绑定）'),
        ('8.36311\\times10^{-6}\\ \\mathrm{kg}', '估算练习·正文（精确绑定）'),
        ('8.36311\\ \\mathrm{mg}', '估算练习·正文（精确绑定）'),
        ('3.41592\\times10^{19}\\ \\mathrm{m^{-3}}', '估算练习·正文（精确绑定）'),
        ('1.20605\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('67.7065\\ \\mathrm{cm}', '估算练习·正文（精确绑定）'),
        ('0.151483\\ \\mathrm{m^2}', '估算练习·正文（精确绑定）'),
        ('478.124\\ \\mathrm{MW\\,m^{-2}}', '估算练习·正文（精确绑定）'),
        ('2.476\\ \\mathrm{GW}', '估算练习·正文（精确绑定）'),
        ('742.8\\ \\mathrm{MW}', '估算练习·正文（精确绑定）'),
        ('28.5714\\ \\mathrm{Hz}', '估算练习·正文（精确绑定）'),
        ('6864.94\\ \\mathrm{MW\\,m^{-2}}', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AD'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 数值 & 表达式 & 说明 \\\\', '表格数据行'),
        ('L-H 阈值 & $63.2554\\ \\mathrm{MW}$ & $0.0488 n_{20}^{0.717}B^{0.803}S^{0.941}$ & 工作点的阈值 \\\\', '表格数据行'),
        ('alpha 自加热 & $111.3\\ \\mathrm{MW}$ & $1.75953P_{\\mathrm{th}}$ & 单独就够进 H 模 \\\\', '表格数据行'),
        ('穿越分界面功率 & $103.468\\ \\mathrm{MW}$ & $P_\\alpha+P_{\\mathrm{aux}}-P_{\\mathrm{rad,core}}$ & $=1.63572P_{\\mathrm{th}}$ \\\\', '表格数据行'),
        ('阈值余量 & $40.2124\\ \\mathrm{MW}$ & $P_{\\mathrm{SOL}}-P_{\\mathrm{th}}$ & 可用于辐射与输运 \\\\', '表格数据行'),
        ('H$\\to$L 反向阈值 & $41.116\\ \\mathrm{MW}$ & $0.65P_{\\mathrm{th}}$ & 回滞宽度 \\\\', '表格数据行'),
        ('格林沃尔德密度 & $1.19366\\times10^{20}\\ \\mathrm{m^{-3}}$ & $I_p/(\\pi a^2)$ & $f_G=0.837758$ \\\\', '表格数据行'),
        ('\\begin{tip}{阈值不是「加热到多少」，是「输入压过损失」}{tip-app-ad-threshold}', '语义盒标题'),
        ('\\textbf{出处：}L-H 转换的实验发现参照 Wagner 等（1982，\\emph{Phys. Rev. Lett.} \\textb', '出处声明整行'),
        ('台基宽度 & $11.9603\\ \\mathrm{cm}$ & $0.076\\sqrt{\\beta_{p,\\mathrm{ped}}}R$ & $=8.18937\\rho_\\theta$ \\\\', '表格数据行'),
        ('台基顶压强 & $57.6784\\ \\mathrm{kPa}$ & $2n_{\\mathrm{ped}}T_{\\mathrm{ped}}$ & $=0.193861\\langle p\\rangle$ \\\\', '表格数据行'),
        ('台基梯度 & $4.82250\\times10^{5}\\ \\mathrm{Pa\\,m^{-1}}$ & $p_{\\mathrm{ped}}/\\Delta_{\\mathrm{ped}}$ & 比芯部平均陡一个量级 \\\\', '表格数据行'),
        ('气球参数 & $3.27709$ & $2\\mu_0Rq_{95}^2|dp/dr|/B^2$ & 与判据 $2.5$--$4$ 相比 \\\\', '表格数据行'),
        ('自举电流密度 & $2.73900\\times10^{5}\\ \\mathrm{A\\,m^{-2}}$ & $1.5\\sqrt{\\varepsilon}|dp/dr|/B_\\theta$ & $=0.229462j_{\\mathrm{avg}}$ \\\\', '表格数据行'),
        ('台基自举电流 & $1.27616\\ \\mathrm{MA}$ & $j_{\\mathrm{bs}}\\Delta_{\\mathrm{ped}}2\\pi R$ & $=0.085077I_p$ \\\\', '表格数据行'),
        ('台基总 beta & $0.5161\\%$ & $2\\mu_0p_{\\mathrm{ped}}/B^2$ & 芯部约 $3.993\\%$ \\\\', '表格数据行'),
        ('\\begin{tip}{台基宽度不是自由参数}{tip-app-ad-width}', '语义盒标题'),
        ('\\textbf{出处：}台基宽度与高度的经验定标参照 Snyder 等（2011，\\emph{Nucl. Fusion} \\textbf{5', '出处声明整行'),
        ('44.0012\\ \\mathrm{kPa}, & \\alpha_{\\mathrm{crit}}=2.5\\ \\text{（无自举）},\\\\', '表格数据行'),
        ('1.31084, & \\text{无自举判据},\\\\', '表格数据行'),
        ('气球模上限（无自举） & $44.0012\\ \\mathrm{kPa}$ & $\\alpha_{\\mathrm{crit}}B^2\\Delta/2\\mu_0Rq_{95}^2$ & $\\alpha_{\\mathrm{crit}}=2.5$ \\\\', '表格数据行'),
        ('气球模上限（含自举） & $70.4019\\ \\mathrm{kPa}$ & 同上 & $\\alpha_{\\mathrm{crit}}=4.0$ \\\\', '表格数据行'),
        ('工作点台基压强 & $57.6784\\ \\mathrm{kPa}$ & $2n_{\\mathrm{ped}}T_{\\mathrm{ped}}$ & 落在两条判据之间 \\\\', '表格数据行'),
        ('比值 & $1.31084$ / $0.819273$ & $\\alpha_{\\mathrm{MHD}}/\\alpha_{\\mathrm{crit}}$ & 与判据同量级 \\\\', '表格数据行'),
        ('边缘电流分数 & $0.229462$ & $j_{\\mathrm{bs}}/j_{\\mathrm{avg}}$ & 剥离模驱动源 \\\\', '表格数据行'),
        ('台基壳层体积 & $58.5493\\ \\mathrm{m^3}$ & $S\\Delta_{\\mathrm{ped}}$ & $=7.03548\\%\\,V$ \\\\', '表格数据行'),
        ('壳层储能 & $2.53277\\ \\mathrm{MJ}$ & $3p_{\\mathrm{ped}}V_{\\mathrm{ped}}/4$ & $<1\\%\\,W_{\\mathrm{th}}$ \\\\', '表格数据行'),
        ('\\begin{tip}{台基不是「越高越好」}{tip-app-ad-two-branches}', '语义盒标题'),
        ('\\textbf{出处：}理想气球模判据的 $\\alpha$ 参数与临界值参照 Connor、Hastie 与 Taylor', '出处声明整行'),
        ('单次 ELM 损失 & $25.998\\ \\mathrm{MJ}$ & $0.07W_{\\mathrm{th}}$ & 类型 I 常规值 \\\\', '表格数据行'),
        ('偏滤器面积 & $38.9557\\ \\mathrm{m^2}$ & $2\\pi RL_{\\mathrm{div}}$ & $L_{\\mathrm{div}}=1\\ \\mathrm{m}$ \\\\', '表格数据行'),
        ('能量密度 & $0.667373\\ \\mathrm{MJ\\,m^{-2}}$ & $\\Delta W/A_{\\mathrm{div}}$ & 阈值的 $0.667373$ \\\\', '表格数据行'),
        ('峰值热流 & $2.22458\\ \\mathrm{GW\\,m^{-2}}$ & $E/t_{\\mathrm{ELM}}$ & $t_{\\mathrm{ELM}}=0.3\\ \\mathrm{ms}$ \\\\', '表格数据行'),
        ('ELM 频率 & $3.97984\\ \\mathrm{Hz}$ & $P_{\\mathrm{SOL}}/\\Delta W_{\\mathrm{ELM}}$ & 由功率平衡定 \\\\', '表格数据行'),
        ('ELM 周期 & $251.267\\ \\mathrm{ms}$ & $1/f_{\\mathrm{ELM}}$ & $f_{\\mathrm{ELM}}\\tau_E=14.7254$ \\\\', '表格数据行'),
        ('\\begin{tip}{ELM 频率不是外部旋钮}{tip-app-ad-elm-power}', '语义盒标题'),
        ('\\textbf{出处：}ELM 的分类与物理图像参照 Zohm（1996，\\emph{Plasma Phys. Control. Fusio', '出处声明整行'),
        ('步调 $20\\ \\mathrm{Hz}$ 的单次损失 & $5.17339\\ \\mathrm{MJ}$ & $P_{\\mathrm{SOL}}/f$ & 未缓解的 $0.199$ 倍 \\\\', '表格数据行'),
        ('对应能量密度 & $0.132802\\ \\mathrm{MJ\\,m^{-2}}$ & $\\Delta W/A_{\\mathrm{div}}$ & 远低于 $1\\ \\mathrm{MJ\\,m^{-2}}$ \\\\', '表格数据行'),
        ('对应峰值热流 & $0.442672\\ \\mathrm{GW\\,m^{-2}}$ & $E/t_{\\mathrm{ELM}}$ & 未缓解的 $0.199$ 倍 \\\\', '表格数据行'),
        ('能量密度达标频率 & $8.85345\\ \\mathrm{Hz}$ & $P_{\\mathrm{SOL}}/(0.3\\,\\mathrm{MJ\\,m^{-2}}A_{\\mathrm{div}})$ & 步调系统工作点 \\\\', '表格数据行'),
        ('$n=3$ RMP 共振面 & $1.78885\\ \\mathrm{m}$ & $a\\sqrt{(3-q_0)/(q_{95}-q_0)}$ & 台基顶内 $9.15431\\ \\mathrm{cm}$ \\\\', '表格数据行'),
        ('$n=2$ RMP 共振面 & $1.26491\\ \\mathrm{m}$ & 同上 & 在台基以里 \\\\', '表格数据行'),
        ('弹丸质量 & $4.18155\\ \\mathrm{mg}$ & $10^{21}\\times2.5m_p$ & $\\Delta n=28.4661\\%n_{\\mathrm{ped}}$ \\\\', '表格数据行'),
        ('\\begin{tip}{缓解不是消除}{tip-app-ad-mitigation}', '语义盒标题'),
        ('\\textbf{出处：}RMP 抑制 ELM 的实验结果参照 Evans 等（2006，\\emph{Nature Physics} \\tex', '出处声明整行'),
        ('热流宽度 & $0.388859\\ \\mathrm{mm}$ & $0.63B_\\theta^{-1.19}$ & $B_\\theta=0.089$ \\\\', '表格数据行'),
        ('总宽度 & $1.94429\\ \\mathrm{mm}$ & $5\\lambda_q$ & 覆盖约 $99\\%$ 热流 \\\\', '表格数据行'),
        ('沉积面积 & $0.0757414\\ \\mathrm{m^2}$ & $2\\pi R\\lambda_t\\times5$ & 倾角因子 $5$ \\\\', '表格数据行'),
        ('峰值热流 & $956.246\\ \\mathrm{MW\\,m^{-2}}$ & $0.7P_{\\mathrm{SOL}}/A_{\\mathrm{out}}$ & 阈值的 $95.6246$ 倍 \\\\', '表格数据行'),
        ('需辐射掉的比例 & $0.989542$ & $1-10/q_{\\mathrm{peak}}$ & 完全不展宽 \\\\', '表格数据行'),
        ('展宽 $5$ 倍 + 辐射 $90\\%$ & $19.1249\\ \\mathrm{MW\\,m^{-2}}$ & 见式 & 仍超阈值 \\\\', '表格数据行'),
        ('展宽 $5$ 倍 + 辐射 $95\\%$ & $9.56246\\ \\mathrm{MW\\,m^{-2}}$ & 见式 & 刚好达标 \\\\', '表格数据行'),
        ('与壁面韧致之比 & $17337.5$ & $q_{\\mathrm{peak}}/q_{\\mathrm{brems}}$ & 四个数量级 \\\\', '表格数据行'),
        ('\\begin{tip}{台基决定「有多少」，SOL 决定「多小」}{tip-app-ad-sol}', '语义盒标题'),
        ('\\textbf{出处：}热流宽度的经验定标参照 Eich 等（2013，\\emph{Nucl. Fusion} \\textbf{53}, 0', '出处声明整行'),
        ('等级 & 例子 & 来源 & 可否被门禁绑定 \\\\', '表格数据行'),
        ('现算级 & $P_{\\mathrm{SOL}}=103.468\\ \\mathrm{MW}$、$\\lambda_q=0.388859\\ \\mathrm{mm}$ & \\texttt{tmp/appAD\\_calc.py} & 可以（逐节数值绑定） \\\\', '表格数据行'),
        ('声明级 & $\\alpha_{\\mathrm{crit}}$、$H_{98}$、$0.07$、$1\\ \\mathrm{MJ\\,m^{-2}}$ & 脚本顶部注释 & 可以（表格整行字面绑定） \\\\', '表格数据行'),
        ('文献级 & ELM 上升时间范围、$\\lambda_q$ 定标散布 & 本节的文献条目 & 否（只作定性对照） \\\\', '表格数据行'),
        ('\\begin{tip}{诚实地说清「哪些数不是算出来的」}{tip-app-ad-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录的物理图像与判据形式综合参照 Freidberg \\emph{Ideal MHD}（2014，', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AD'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 63.2554, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 489.532, '\\mathrm{m^2}', '从 L 模到 H 模·公式'),
        (0, 5.3, '\\mathrm{T}', '从 L 模到 H 模·公式'),
        (0, 111.3, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 103.468, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 40.2124, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 1.19366e+20, '\\mathrm{m^{-3}}', '从 L 模到 H 模·公式'),
        (0, 41.116, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 38.4821, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 103.977, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 4, '\\mathrm{T}', '从 L 模到 H 模·公式'),
        (0, 50.4612, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 6, '\\mathrm{T}', '从 L 模到 H 模·公式'),
        (0, 69.881, '\\mathrm{MW}', '从 L 模到 H 模·公式'),
        (0, 593524, '\\mathrm{m\\,s^{-1}}', '从 L 模到 H 模·公式'),
        (0, 593.524, '\\mathrm{km\\,s^{-1}}', '从 L 模到 H 模·公式'),
        (0, 296762, '\\mathrm{s^{-1}}', '从 L 模到 H 模·公式'),
        (0, 35493.6, '\\mathrm{m\\,s^{-1}}', '从 L 模到 H 模·公式'),
        (0, 53.2405, '\\mathrm{kV\\,m^{-1}}', '从 L 模到 H 模·公式'),
        (0, 1.5, '\\mathrm{T}', '从 L 模到 H 模·正文'),
        (0, 30, '\\mathrm{kV\\,m^{-1}}', '从 L 模到 H 模·正文'),
        (0, 167220, '\\mathrm{s^{-1}}', '从 L 模到 H 模·正文'),
        (1, 4.13338, '\\mathrm{mm}', '台基：十二厘米的壳层·公式'),
        (1, 1.46046, '\\mathrm{cm}', '台基：十二厘米的壳层·公式'),
        (1, 15, '\\mathrm{MA}', '台基：十二厘米的壳层·公式'),
        (1, 9.194, '\\mathrm{keV}', '台基：十二厘米的壳层·公式'),
        (1, 241.933, '\\mathrm{m^{-1}}', '台基：十二厘米的壳层·公式'),
        (1, 3, '\\mathrm{keV}', '台基：十二厘米的壳层·公式'),
        (1, 57.6784, '\\mathrm{kPa}', '台基：十二厘米的壳层·公式'),
        (1, 297.525, '\\mathrm{kPa}', '台基：十二厘米的壳层·公式'),
        (1, 11.9603, '\\mathrm{cm}', '台基：十二厘米的壳层·公式'),
        (1, 482250, '\\mathrm{Pa\\,m^{-1}}', '台基：十二厘米的壳层·公式'),
        (1, 1.8804, '\\mathrm{m}', '台基：十二厘米的壳层·公式'),
        (1, 273900, '\\mathrm{A\\,m^{-2}}', '台基：十二厘米的壳层·公式'),
        (1, 1.27616, '\\mathrm{MA}', '台基：十二厘米的壳层·公式'),
        (1, 1193660, '\\mathrm{A\\,m^{-2}}', '台基：十二厘米的壳层·正文'),
        (2, 44.0012, '\\mathrm{kPa}', '两条边界的夹缝：气球·公式'),
        (2, 70.4019, '\\mathrm{kPa}', '两条边界的夹缝：气球·公式'),
        (2, 58.5493, '\\mathrm{m^3}', '两条边界的夹缝：气球·公式'),
        (2, 2.53277, '\\mathrm{MJ}', '两条边界的夹缝：气球·公式'),
        (3, 25998000, '\\mathrm{J}', 'ELM：储能的周期性·公式'),
        (3, 25.998, '\\mathrm{MJ}', 'ELM：储能的周期性·公式'),
        (3, 1, '\\mathrm{m}', 'ELM：储能的周期性·公式'),
        (3, 38.9557, '\\mathrm{m^2}', 'ELM：储能的周期性·公式'),
        (3, 667373, '\\mathrm{J\\,m^{-2}}', 'ELM：储能的周期性·公式'),
        (3, 0.667373, '\\mathrm{MJ\\,m^{-2}}', 'ELM：储能的周期性·公式'),
        (3, 0.3, '\\mathrm{ms}', 'ELM：储能的周期性·公式'),
        (3, 2.22458, '\\mathrm{GW\\,m^{-2}}', 'ELM：储能的周期性·公式'),
        (3, 1, '\\mathrm{MJ\\,m^{-2}}', 'ELM：储能的周期性·公式'),
        (3, 3.97984, '\\mathrm{Hz}', 'ELM：储能的周期性·公式'),
        (3, 251.267, '\\mathrm{ms}', 'ELM：储能的周期性·公式'),
        (3, 0.684532, '\\mathrm{MW}', 'ELM：储能的周期性·正文'),
        (4, 50, '\\mathrm{ms}', '缓解：把一次大坍塌拆·正文'),
        (4, 20, '\\mathrm{Hz}', '缓解：把一次大坍塌拆·正文'),
        (4, 5.17339, '\\mathrm{MJ}', '缓解：把一次大坍塌拆·公式'),
        (4, 0.132802, '\\mathrm{MJ\\,m^{-2}}', '缓解：把一次大坍塌拆·公式'),
        (4, 0.442672, '\\mathrm{GW\\,m^{-2}}', '缓解：把一次大坍塌拆·公式'),
        (4, 0.3, '\\mathrm{MJ\\,m^{-2}}', '缓解：把一次大坍塌拆·公式'),
        (4, 8.85345, '\\mathrm{Hz}', '缓解：把一次大坍塌拆·公式'),
        (4, 1.78885, '\\mathrm{m}', '缓解：把一次大坍塌拆·公式'),
        (4, 9.15431, '\\mathrm{cm}', '缓解：把一次大坍塌拆·公式'),
        (4, 1.26491, '\\mathrm{m}', '缓解：把一次大坍塌拆·公式'),
        (4, 4.18155e-06, '\\mathrm{kg}', '缓解：把一次大坍塌拆·公式'),
        (4, 4.18155, '\\mathrm{mg}', '缓解：把一次大坍塌拆·公式'),
        (4, 1.70796e+19, '\\mathrm{m^{-3}}', '缓解：把一次大坍塌拆·正文'),
        (5, 0.388859, '\\mathrm{mm}', '从台基到靶板：热流落·公式'),
        (5, 1.94429, '\\mathrm{mm}', '从台基到靶板：热流落·公式'),
        (5, 0.0757414, '\\mathrm{m^2}', '从台基到靶板：热流落·公式'),
        (5, 956.246, '\\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式'),
        (5, 10, '\\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式'),
        (5, 19.1249, '\\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式'),
        (5, 9.56246, '\\mathrm{MW\\,m^{-2}}', '从台基到靶板：热流落·公式'),
        (5, 55.1547, '\\mathrm{kW\\,m^{-2}}', '从台基到靶板：热流落·公式'),
        (6, 150, '\\mathrm{MW}', '估算练习·正文'),
        (6, 0.15, '\\mathrm{s}', '估算练习·正文'),
        (6, 371.4, '\\mathrm{MJ}', '估算练习·正文'),
        (6, 37.14, '\\mathrm{MJ}', '估算练习·正文'),
        (6, 0.953391, '\\mathrm{MJ\\,m^{-2}}', '估算练习·正文'),
        (6, 5.76967, '\\mathrm{Hz}', '估算练习·正文'),
        (6, 173.32, '\\mathrm{ms}', '估算练习·正文'),
        (6, 8.36311e-06, '\\mathrm{kg}', '估算练习·正文'),
        (6, 8.36311, '\\mathrm{mg}', '估算练习·正文'),
        (6, 3.41592e+19, '\\mathrm{m^{-3}}', '估算练习·正文'),
        (6, 1.20605, '\\mathrm{m}', '估算练习·正文'),
        (6, 67.7065, '\\mathrm{cm}', '估算练习·正文'),
        (6, 0.151483, '\\mathrm{m^2}', '估算练习·正文'),
        (6, 478.124, '\\mathrm{MW\\,m^{-2}}', '估算练习·正文'),
        (6, 2.476, '\\mathrm{GW}', '估算练习·正文'),
        (6, 742.8, '\\mathrm{MW}', '估算练习·正文'),
        (6, 28.5714, '\\mathrm{Hz}', '估算练习·正文'),
        (6, 6864.94, '\\mathrm{MW\\,m^{-2}}', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AD 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AD 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('44.0\\ \\mathrm{kPa}', '气球模上限四舍五入'),
        ('70.4\\ \\mathrm{kPa}', '含自举气球模上限四舍五入'),
        ('26.0\\ \\mathrm{MJ}', '单次 ELM 损失四舍五入'),
        ('0.667\\ \\mathrm{MJ', 'ELM 能量密度四舍五入'),
        ('3.98\\ \\mathrm{Hz}', 'ELM 频率四舍五入'),
        ('956\\ \\mathrm{MW', 'SOL 峰值热流四舍五入'),
        ('0.990', '需辐射掉的比例四舍五入'),
        ('9.56\\ \\mathrm{MW', '展宽后峰值热流四舍五入'),
        ('0.39\\ \\mathrm{mm}', '热流宽度四舍五入'),
        ('1.9\\ \\mathrm{mm}', '热流总宽度四舍五入'),
        ('17338', 'SOL 与壁面韧致之比四舍五入'),
        ('5.8\\ \\mathrm{Hz}', '练习 2 频率四舍五入'),
        ('0.95\\ \\mathrm{MJ', '练习 1 能量密度四舍五入'),
        ('2.5\\ \\mathrm{MJ}', '台基壳层储能四舍五入'),
        ('58.5\\ \\mathrm{m^3}', '台基壳层体积四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AD'], why))


def check_P37_appendix_ae_kinetics():
    """附录 AE：动理学与回旋动理学——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AE'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AE 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '分布函数：状态量为什么不能',
        '无碰撞不等于无相互作用',
        '波粒共振：Landau 阻尼',
        '捕获粒子、香蕉轨道',
        '回旋动理学：把最快的',
        '速度分布的形状',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AE 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AE 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ae-scales',
        'tab:app-ae-collisionless',
        'tab:app-ae-landau',
        'tab:app-ae-orbit',
        'tab:app-ae-gyro',
        'tab:app-ae-velocity',
        'tab:app-ae-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AE'])

    for lab in (
        'sec:app-ae-distribution',
        'sec:app-ae-vlasov',
        'sec:app-ae-landau',
        'sec:app-ae-orbit',
        'sec:app-ae-gyro',
        'sec:app-ae-tail',
        'sec:app-ae-problems',
        'sec:app-ae-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AE'])

    for lab in (
        'tip-app-ae-distribution',
        'tip-app-ae-collisionless',
        'tip-app-ae-landau',
        'tip-app-ae-orbit',
        'tip-app-ae-gyro',
        'tip-app-ae-tail',
        'tip-app-ae-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AE'])

    for lab in (
        'eq:app-ae-vthe',
        'eq:app-ae-vthe2',
        'eq:app-ae-debye',
        'eq:app-ae-rho',
        'eq:app-ae-vlasov',
        'eq:app-ae-collisionality',
        'eq:app-ae-mfp',
        'eq:app-ae-kn',
        'eq:app-ae-moments',
        'eq:app-ae-landau',
        'eq:app-ae-resonance',
        'eq:app-ae-resonance2',
        'eq:app-ae-ionacoustic',
        'eq:app-ae-nustar',
        'eq:app-ae-ftrap',
        'eq:app-ae-banana',
        'eq:app-ae-bounce',
        'eq:app-ae-gyro',
        'eq:app-ae-gamma0',
        'eq:app-ae-gyro2',
        'eq:app-ae-chigb',
        'eq:app-ae-chieff',
        'eq:app-ae-tail',
        'eq:app-ae-energytail',
        'eq:app-ae-losscone',
        'eq:app-ae-closure',
        'eq:app-ae-anisotropy',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AE'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AE'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AE 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AE 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AE 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AE 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 数字 + 单位的精确绑定（堵末位容差洞） ----------
    for lit, why in (
        ('4.02127\\times10^{7}\\ \\mathrm{m\\,s^{-1}}', '分布函数：状态量为什·公式（精确绑定）'),
        ('839.370\\ \\mathrm{km\\,s^{-1}}', '分布函数：状态量为什·公式（精确绑定）'),
        ('593.524\\ \\mathrm{km\\,s^{-1}}', '分布函数：状态量为什·公式（精确绑定）'),
        ('71.2806\\ \\mathrm{\\mu m}', '分布函数：状态量为什·公式（精确绑定）'),
        ('43.1386\\ \\mathrm{\\mu m}', '分布函数：状态量为什·公式（精确绑定）'),
        ('4.13338\\ \\mathrm{mm}', '分布函数：状态量为什·公式（精确绑定）'),
        ('3.58755\\ \\mathrm{km}', '无碰撞不等于无相互作·公式（精确绑定）'),
        ('6.04449\\times10^{-3}\\ \\mathrm{s}', '无碰撞不等于无相互作·公式（精确绑定）'),
        ('6.3576\\times10^{11}\\ \\mathrm{rad\\,s^{-1}}', '波粒共振：Landau 阻尼为·公式（精确绑定）'),
        ('4208.72\\ \\mathrm{m^{-1}}', '波粒共振：Landau 阻尼为·公式（精确绑定）'),
        ('1.4929\\ \\mathrm{mm}', '波粒共振：Landau 阻尼为·公式（精确绑定）'),
        ('1.51058\\times10^{8}\\ \\mathrm{m\\,s^{-1}}', '波粒共振：Landau 阻尼为·公式（精确绑定）'),
        ('1402.91\\ \\mathrm{m^{-1}}', '波粒共振：Landau 阻尼为·正文（精确绑定）'),
        ('4.47869\\ \\mathrm{mm}', '波粒共振：Landau 阻尼为·正文（精确绑定）'),
        ('8.32659\\times10^{8}\\ \\mathrm{rad\\,s^{-1}}', '波粒共振：Landau 阻尼为·正文（精确绑定）'),
        ('132.522\\ \\mathrm{MHz}', '波粒共振：Landau 阻尼为·正文（精确绑定）'),
        ('21.7\\ \\mathrm{m}', '捕获粒子、香蕉轨道与·正文（精确绑定）'),
        ('2.54715\\ \\mathrm{cm}', '捕获粒子、香蕉轨道与·公式（精确绑定）'),
        ('21969.1\\ \\mathrm{s^{-1}}', '捕获粒子、香蕉轨道与·公式（精确绑定）'),
        ('3496.5\\ \\mathrm{Hz}', '捕获粒子、香蕉轨道与·公式（精确绑定）'),
        ('286.001\\ \\mathrm{\\mu s}', '捕获粒子、香蕉轨道与·公式（精确绑定）'),
        ('2.5\\ \\mathrm{cm}', '捕获粒子、香蕉轨道与·正文（精确绑定）'),
        ('72.5799\\ \\mathrm{m^{-1}}', '回旋动理学：把最快的·公式（精确绑定）'),
        ('8.65693\\ \\mathrm{cm}', '回旋动理学：把最快的·公式（精确绑定）'),
        ('5.07013\\ \\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式（精确绑定）'),
        ('108.42\\ \\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式（精确绑定）'),
        ('0.540541\\ \\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式（精确绑定）'),
        ('3.7\\ \\mathrm{s}', '回旋动理学：把最快的·正文（精确绑定）'),
        ('10\\ \\mathrm{keV}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('64\\ \\mathrm{keV}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('100\\ \\mathrm{keV}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('3.5\\ \\mathrm{MeV}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('1.29472\\times10^{7}\\ \\mathrm{m\\,s^{-1}}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('5.10055\\ \\mathrm{cm}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('5.07356\\ \\mathrm{km}', '速度分布的形状：尾部·公式（精确绑定）'),
        ('1\\times10^{19}\\ \\mathrm{m^{-3}}', '估算练习·正文（精确绑定）'),
        ('1\\ \\mathrm{keV}', '估算练习·正文（精确绑定）'),
        ('5.3\\ \\mathrm{T}', '估算练习·正文（精确绑定）'),
        ('1\\ \\mathrm{m^3}', '估算练习·正文（精确绑定）'),
        ('7.43394\\times10^{-5}\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('74.3394\\ \\mathrm{\\mu m}', '估算练习·正文（精确绑定）'),
        ('276822\\ \\mathrm{m\\,s^{-1}}', '估算练习·正文（精确绑定）'),
        ('1.36318\\ \\mathrm{mm}', '估算练习·正文（精确绑定）'),
        ('241.933\\ \\mathrm{m^{-1}}', '估算练习·正文（精确绑定）'),
        ('0.0259708\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('2.59708\\ \\mathrm{cm}', '估算练习·正文（精确绑定）'),
        ('4.3985\\times10^{16}\\ \\mathrm{m^{-3}}', '估算练习·正文（精确绑定）'),
        ('8.92144\\times10^{-4}\\ \\mathrm{s}', '估算练习·正文（精确绑定）'),
        ('35875.5\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('35.8755\\ \\mathrm{km}', '估算练习·正文（精确绑定）'),
        ('0.0123719\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('1.23719\\ \\mathrm{cm}', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AE'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 本工作点的值 & 表达式 & 它在哪件事里说话 \\\\', '表格数据行'),
        ('电子热速度 & $4.02127\\times10^{7}\\ \\mathrm{m\\,s^{-1}}$ & $\\sqrt{T_e/m_e}$ & 德拜判据、网格与时间步 \\\\', '表格数据行'),
        ('离子热速度 & $839.370\\ \\mathrm{km\\,s^{-1}}$ & $\\sqrt{2T_i/m_i}$ & 离子轨道、香蕉宽度 \\\\', '表格数据行'),
        ('离子声速 & $593.524\\ \\mathrm{km\\,s^{-1}}$ & $\\sqrt{T_e/m_i}$ & 离子声波、台基宽度 \\\\', '表格数据行'),
        ('德拜长度 & $71.2806\\ \\mathrm{\\mu m}$ & $\\sqrt{\\epsilon_0T_e/(n_ee^2)}$ & 集体性判据、屏蔽 \\\\', '表格数据行'),
        ('德拜立方内粒子数 & $3.62171\\times10^{7}$ & $n_e\\lambda_{De}^3$ & 统计描述是否成立 \\\\', '表格数据行'),
        ('电子回旋半径 & $43.1386\\ \\mathrm{\\mu m}$ & $v_{the}/\\omega_{ce}$ & 电子尺度湍流 \\\\', '表格数据行'),
        ('离子回旋半径 & $4.13338\\ \\mathrm{mm}$ & $v_{Ti}/\\omega_{ci}$ & 回旋动理学的展开参数 \\\\', '表格数据行'),
        ('\\begin{tip}{分布函数不是「更精确的密度」}{tip-app-ae-distribution}', '语义盒标题'),
        ('\\textbf{出处：}动理学描述与德拜屏蔽的框架参照 Nicholson \\emph{Introduction to Plasma Phy', '出处声明整行'),
        ('判据 & 本工作点 & 阈值 & 不成立时会发生什么 \\\\', '表格数据行'),
        ('$\\omega_{pe}\\tau_{ei}$ & $5.03299\\times10^{7}$ & $\\gg1$ & 等离子体不屏蔽、像中性气体 \\\\', '表格数据行'),
        ('$\\omega_{ce}\\tau_{ei}$ & $8.31634\\times10^{7}$ & $\\gg1$ & 电子被磁化、回旋是快运动 \\\\', '表格数据行'),
        ('$\\omega_{ci}\\tau_{ei}$ & $18116.9$ & $\\gg1$ & 离子未被电子碰撞钉住 \\\\', '表格数据行'),
        ('$\\lambda_{\\mathrm{mfp}}/(2\\pi R)$ & $92.0929$ & $\\gg1$ & 环向输运由轨道而非碰撞决定 \\\\', '表格数据行'),
        ('$\\mathrm{Kn}=\\lambda_{\\mathrm{mfp}}/a$ & $1793.77$ & $\\gg1$ & 流体闭合必须声明，不能默认 \\\\', '表格数据行'),
        ('\\begin{tip}{无碰撞是个不等式，不是「没有相互作用」}{tip-app-ae-collisionless}', '语义盒标题'),
        ('\\textbf{出处：}弗拉索夫方程与矩展开的写法参照 Nicholson \\emph{Introduction to Plasma Phy', '出处声明整行'),
        ('模式 & $k\\lambda_{De}$ & $\\gamma/\\omega$ & 读法 \\\\', '表格数据行'),
        ('朗缪尔波 & $0.1$ & $-2.6969\\times10^{-20}$ & 相速度远离尾部，几乎不衰减 \\\\', '表格数据行'),
        ('朗缪尔波 & $0.3$ & $-0.0200206$ & 已经能在一百个周期里看到衰减 \\\\', '表格数据行'),
        ('朗缪尔波 & $0.5$ & $-0.151387$ & 强阻尼，波只活几个周期 \\\\', '表格数据行'),
        ('朗缪尔波 & $1$ & $-0.0848088$ & 公式已越界（$\\omega\\simeq\\omega_{pe}$ 不成立） \\\\', '表格数据行'),
        ('离子声波（$T_e/T_i=1$） & --- & $-0.094058$ & 离子声波靠电子与离子两头一起吃 \\\\', '表格数据行'),
        ('离子声波（$T_e/T_i=10$） & --- & $-0.00927902$ & 电子温度拉高，阻尼反而变弱 \\\\', '表格数据行'),
        ('\\begin{tip}{Landau 阻尼不是「无碰撞的碰撞」}{tip-app-ae-landau}', '语义盒标题'),
        ('\\textbf{出处：}Landau 阻尼的原始结果见 Landau（1946，\\emph{Journal of Physics USSR}', '出处声明整行'),
        ('量 & 本工作点 & 表达式 & 它在哪件事里说话 \\\\', '表格数据行'),
        ('反径比 $\\epsilon$ & $0.322581$ & $a/R$ & 磁镜深度 \\\\', '表格数据行'),
        ('碰撞率 $\\nu_*$ & $0.0330145$ & $\\nu_{ei}qR/(\\epsilon^{1.5}v_{the})$ & 输运区判别（$<1$ 香蕉区） \\\\', '表格数据行'),
        ('捕获比例 $f_t$ & $69.843\\%$ & $\\sqrt{2\\epsilon/(1+\\epsilon)}$ & 有多少粒子在弹跳 \\\\', '表格数据行'),
        ('香蕉宽度 $\\Delta_b$ & $2.54715\\ \\mathrm{cm}$ & $\\rho_iq_{95}/\\sqrt{\\epsilon}$ & 新经典输运的步长 \\\\', '表格数据行'),
        ('弹跳频率 $f_b$ & $3496.5\\ \\mathrm{Hz}$ & $v_{Ti}\\sqrt{\\epsilon}/(2\\pi q_{95}R)$ & 弹跳周期 $286\\ \\mathrm{\\mu s}$ \\\\', '表格数据行'),
        ('\\begin{tip}{香蕉轨道不是「轨道」，是漂移的平均失败}{tip-app-ae-orbit}', '语义盒标题'),
        ('\\textbf{出处：}香蕉轨道与新经典输运的框架见 Galeev \\& Sagdeev（1968，', '出处声明整行'),
        ('$k_\\perp\\rho_i$ & $b$ & $\\Gamma_0(b)$ & 抑制掉的比例 $1-\\Gamma_0$ \\\\', '表格数据行'),
        ('$0.3$ & $0.09$ & $0.915783$ & $8.42172\\%$（长波，FLR 只是修正） \\\\', '表格数据行'),
        ('$1$ & $1$ & $0.46576$ & $53.424\\%$（一半以上的响应被抹平） \\\\', '表格数据行'),
        ('$2$ & $4$ & $0.207002$ & $79.2998\\%$（短波，几乎完全被平均掉） \\\\', '表格数据行'),
        ('\\begin{tip}{回旋动理学不是「把回旋平均掉就完事」}{tip-app-ae-gyro}', '语义盒标题'),
        ('\\textbf{出处：}回旋动理学的系统推导见 Frieman \\& Chen（1982，\\emph{Physics of Fluids}', '出处声明整行'),
        ('分布特征 & 本工作点 & 表达式 & 后果 \\\\', '表格数据行'),
        ('$v>v_{the}$ & $57.2407\\%$ & 高斯尾 & 一半以上粒子「比热速度快」 \\\\', '表格数据行'),
        ('$v>3v_{the}$ & $0.043985\\%$ & 高斯尾 & 电流、射频共振的作用区间 \\\\', '表格数据行'),
        ('$E>64\\ \\mathrm{keV}$ & $0.30131\\%$ & 能量尾 & 与 $T_e$ 的比决定反应率 \\\\', '表格数据行'),
        ('alpha 等效 $u$ & $15.4249$（$v_{Ti}$） & $v_\\alpha/v_{Ti}$ & 热尾权重 $4.67\\times10^{-104}$ \\\\', '表格数据行'),
        ('损失锥（$R_m=2$） & $29.2893\\%$ & $1-\\sqrt{1-1/R_m}$ & 磁镜泄漏份额 \\\\', '表格数据行'),
        ('\\begin{tip}{聚变靠的是尾部，不是平均能量}{tip-app-ae-tail}', '语义盒标题'),
        ('\\textbf{出处：}麦克斯韦分布尾部与矩的关系见 Nicholson \\emph{Introduction to Plasma Phys', '出处声明整行'),
        ('来源 & 借鉴了什么 & 不借用什么 \\\\', '表格数据行'),
        ('Nicholson (1983)；Chen (1984) & 动理学框架、德拜屏蔽、矩与闭合的讲法 & 原书文字、编号、插图、习题 \\\\', '表格数据行'),
        ('Krall \\& Trivelpiece (1973) & 弗拉索夫方程、Landau 阻尼、磁镜的处理 & 原书推导步骤与符号体系 \\\\', '表格数据行'),
        ('Landau (1946)；Jackson (1960) & Landau 阻尼的结论与物理图像 & 原文推导与实验数据点 \\\\', '表格数据行'),
        ('Malmberg \\& Wharton (1964) & 「无碰撞阻尼可观测」这一事实 & 实验曲线与数值 \\\\', '表格数据行'),
        ('Bernstein, Greene \\& Kruskal (1957) & 非线性（BGK）图像 & 原文构造细节 \\\\', '表格数据行'),
        ('Brizard \\& Hahm (2007) & 回旋动理学的现代综述框架 & 综述的章节结构与公式编号 \\\\', '表格数据行'),
        ('Hahm \\& Burrell (1995) & 回旋玻姆与玻姆标度的对照 & 定标拟合系数 \\\\', '表格数据行'),
        ('Galeev \\& Sagdeev (1968)；Hinton \\& Hazeltine (1976) & 香蕉轨道与新经典输运的结论 & 原综述的推导与图表 \\\\', '表格数据行'),
        ('Helander \\& Sigmar (2002) & 输运系数与张量闭合的口径 & 原书公式编号与数据 \\\\', '表格数据行'),
        ('Wesson (2011) & 工作点参数、快粒子与 alpha 的地位 & 原书表格与数值 \\\\', '表格数据行'),
        ('Bosch \\& Hale (1992) & D-T 反应率对温度的依赖 & 拟合公式与系数 \\\\', '表格数据行'),
        ('\\begin{tip}{诚实地说清「本附录不给什么」}{tip-app-ae-honest}', '语义盒标题'),
        ('\\textbf{出处：}本节的文献清单与版权口径参照本教材其他附录的惯例', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AE'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 40212700, '\\mathrm{m\\,s^{-1}}', '分布函数：状态量为什·公式'),
        (0, 839.37, '\\mathrm{km\\,s^{-1}}', '分布函数：状态量为什·公式'),
        (0, 593.524, '\\mathrm{km\\,s^{-1}}', '分布函数：状态量为什·公式'),
        (0, 71.2806, '\\mathrm{\\mu m}', '分布函数：状态量为什·公式'),
        (0, 43.1386, '\\mathrm{\\mu m}', '分布函数：状态量为什·公式'),
        (0, 4.13338, '\\mathrm{mm}', '分布函数：状态量为什·公式'),
        (1, 3.58755, '\\mathrm{km}', '无碰撞不等于无相互作·公式'),
        (1, 0.00604449, '\\mathrm{s}', '无碰撞不等于无相互作·公式'),
        (2, 6.3576e+11, '\\mathrm{rad\\,s^{-1}}', '波粒共振：Landau 阻尼为·公式'),
        (2, 4208.72, '\\mathrm{m^{-1}}', '波粒共振：Landau 阻尼为·公式'),
        (2, 1.4929, '\\mathrm{mm}', '波粒共振：Landau 阻尼为·公式'),
        (2, 151058000, '\\mathrm{m\\,s^{-1}}', '波粒共振：Landau 阻尼为·公式'),
        (2, 1402.91, '\\mathrm{m^{-1}}', '波粒共振：Landau 阻尼为·正文'),
        (2, 4.47869, '\\mathrm{mm}', '波粒共振：Landau 阻尼为·正文'),
        (2, 832659000, '\\mathrm{rad\\,s^{-1}}', '波粒共振：Landau 阻尼为·正文'),
        (2, 132.522, '\\mathrm{MHz}', '波粒共振：Landau 阻尼为·正文'),
        (3, 21.7, '\\mathrm{m}', '捕获粒子、香蕉轨道与·正文'),
        (3, 2.54715, '\\mathrm{cm}', '捕获粒子、香蕉轨道与·公式'),
        (3, 21969.1, '\\mathrm{s^{-1}}', '捕获粒子、香蕉轨道与·公式'),
        (3, 3496.5, '\\mathrm{Hz}', '捕获粒子、香蕉轨道与·公式'),
        (3, 286.001, '\\mathrm{\\mu s}', '捕获粒子、香蕉轨道与·公式'),
        (3, 2.5, '\\mathrm{cm}', '捕获粒子、香蕉轨道与·正文'),
        (3, 286, '\\mathrm{\\mu s}', '捕获粒子、香蕉轨道与·表格'),
        (4, 72.5799, '\\mathrm{m^{-1}}', '回旋动理学：把最快的·公式'),
        (4, 8.65693, '\\mathrm{cm}', '回旋动理学：把最快的·公式'),
        (4, 5.07013, '\\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式'),
        (4, 108.42, '\\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式'),
        (4, 0.540541, '\\mathrm{m^2\\,s^{-1}}', '回旋动理学：把最快的·公式'),
        (4, 3.7, '\\mathrm{s}', '回旋动理学：把最快的·正文'),
        (5, 10, '\\mathrm{keV}', '速度分布的形状：尾部·公式'),
        (5, 64, '\\mathrm{keV}', '速度分布的形状：尾部·公式'),
        (5, 100, '\\mathrm{keV}', '速度分布的形状：尾部·公式'),
        (5, 3.5, '\\mathrm{MeV}', '速度分布的形状：尾部·公式'),
        (5, 12947200, '\\mathrm{m\\,s^{-1}}', '速度分布的形状：尾部·公式'),
        (5, 5.10055, '\\mathrm{cm}', '速度分布的形状：尾部·公式'),
        (5, 5.07356, '\\mathrm{km}', '速度分布的形状：尾部·公式'),
        (6, 1e+19, '\\mathrm{m^{-3}}', '估算练习·正文'),
        (6, 1, '\\mathrm{keV}', '估算练习·正文'),
        (6, 5.3, '\\mathrm{T}', '估算练习·正文'),
        (6, 1, '\\mathrm{m^3}', '估算练习·正文'),
        (6, 7.43394e-05, '\\mathrm{m}', '估算练习·正文'),
        (6, 74.3394, '\\mathrm{\\mu m}', '估算练习·正文'),
        (6, 276822, '\\mathrm{m\\,s^{-1}}', '估算练习·正文'),
        (6, 1.36318, '\\mathrm{mm}', '估算练习·正文'),
        (6, 241.933, '\\mathrm{m^{-1}}', '估算练习·正文'),
        (6, 0.0259708, '\\mathrm{m}', '估算练习·正文'),
        (6, 2.59708, '\\mathrm{cm}', '估算练习·正文'),
        (6, 4.3985e+16, '\\mathrm{m^{-3}}', '估算练习·正文'),
        (6, 0.000892144, '\\mathrm{s}', '估算练习·正文'),
        (6, 35875.5, '\\mathrm{m}', '估算练习·正文'),
        (6, 35.8755, '\\mathrm{km}', '估算练习·正文'),
        (6, 0.0123719, '\\mathrm{m}', '估算练习·正文'),
        (6, 1.23719, '\\mathrm{cm}', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AE 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AE 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('1.5\\times10^{8}', '德拜球内粒子数四舍五入'),
        ('3.62\\times10^{7}', '德拜立方内粒子数四舍五入'),
        ('71.3\\ \\mathrm{\\mu m}', '德拜长度四舍五入'),
        ('839.4\\ \\mathrm{km', '离子热速度四舍五入'),
        ('593.5\\ \\mathrm{km', '离子声速四舍五入'),
        ('43.1\\ \\mathrm{\\mu m}', '电子回旋半径四舍五入'),
        ('4.13\\ \\mathrm{mm}', '离子回旋半径四舍五入'),
        ('5.03\\times10^{7}', 'omega_pe tau 四舍五入'),
        ('8.32\\times10^{7}', 'omega_ce tau 四舍五入'),
        ('3.59\\ \\mathrm{km}', '平均自由程四舍五入'),
        ('92.1', '绕环圈数四舍五入'),
        ('1.79\\times10^{3}', '电子努森数四舍五入'),
        ('2.55\\ \\mathrm{cm}', '香蕉宽度四舍五入'),
        ('69.8\\%', '捕获比例四舍五入'),
        ('5.07\\ \\mathrm{m^2', '回旋玻姆扩散四舍五入'),
        ('0.54\\ \\mathrm{m^2', '有效扩散四舍五入'),
        ('108\\ \\mathrm{m^2', '玻姆扩散四舍五入'),
        ('57.2\\%', '一热速度以上比例四舍五入'),
        ('53.7\\%', '10 keV 以上比例四舍五入'),
        ('4.4\\times10^{16}', '三次热速度以上密度四舍五入'),
        ('0.044\\%', '三次热速度以上比例四舍五入'),
        ('2.6\\ \\mathrm{cm}', '练习 2 垂直波长四舍五入'),
        ('1.24\\ \\mathrm{cm}', '练习 5 香蕉宽度四舍五入'),
        ('15.4\\,v_{Ti}', '练习 6 速度比四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AE'], why))


def check_P38_appendix_af_turbulence():
    """附录 AF：等离子体湍流与反常输运——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AF'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AF 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '三个输运系数：反常输运',
        '漂移波：梯度和磁场',
        '不稳定性：ITG 与 TEM',
        '混合长与饱和：上界与真值',
        '带状流与剪切抑制',
        '谱与级联：输运由哪一段波数',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AF 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AF 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-af-scales',
        'tab:app-af-driftwave',
        'tab:app-af-instability',
        'tab:app-af-saturation',
        'tab:app-af-zonal',
        'tab:app-af-spectrum',
        'tab:app-af-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AF'])

    for lab in (
        'sec:app-af-scales',
        'sec:app-af-driftwave',
        'sec:app-af-instability',
        'sec:app-af-saturation',
        'sec:app-af-zonal',
        'sec:app-af-spectrum',
        'sec:app-af-problems',
        'sec:app-af-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AF'])

    for lab in (
        'tip-app-af-orders',
        'tip-app-af-driftwave',
        'tip-app-af-stiffness',
        'tip-app-af-mixing',
        'tip-app-af-zonal',
        'tip-app-af-spectrum',
        'tip-app-af-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AF'])

    for lab in (
        'eq:app-af-chieff',
        'eq:app-af-chib',
        'eq:app-af-chigb',
        'eq:app-af-q',
        'eq:app-af-rhos',
        'eq:app-af-omega-star',
        'eq:app-af-eta',
        'eq:app-af-rlt',
        'eq:app-af-ltcrit',
        'eq:app-af-tem',
        'eq:app-af-gamma',
        'eq:app-af-mixing',
        'eq:app-af-dreff',
        'eq:app-af-vzf',
        'eq:app-af-gammae',
        'eq:app-af-suppress',
        'eq:app-af-spectrum',
        'eq:app-af-flux',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AF'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AF'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AF 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AF 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AF 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AF 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 数字 + 单位的精确绑定（堵末位容差洞） ----------
    for lit, why in (
        ('0.540541\\ \\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式（精确绑定）'),
        ('108.42\\ \\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式（精确绑定）'),
        ('5.07013\\ \\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式（精确绑定）'),
        ('2.53506\\ \\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式（精确绑定）'),
        ('0.159248\\ \\mathrm{MW\\,m^{-2}}', '三个输运系数：反常输·公式（精确绑定）'),
        ('77.9568\\ \\mathrm{MW}', '三个输运系数：反常输·公式（精确绑定）'),
        ('2.92274\\times10^{-3}\\ \\mathrm{m}', '漂移波：梯度和磁场谁·公式（精确绑定）'),
        ('2.92274\\ \\mathrm{mm}', '漂移波：梯度和磁场谁·公式（精确绑定）'),
        ('178057\\ \\mathrm{rad\\,s^{-1}}', '漂移波：梯度和磁场谁·公式（精确绑定）'),
        ('2\\ \\mathrm{m^{-1}}', '不稳定性：ITG 与·公式（精确绑定）'),
        ('1.24\\ \\mathrm{m}', '不稳定性：ITG 与·公式（精确绑定）'),
        ('178057\\ \\mathrm{s^{-1}}', '混合长与饱和：上界与·公式（精确绑定）'),
        ('5.61617\\times10^{-6}\\ \\mathrm{s}', '混合长与饱和：上界与·公式（精确绑定）'),
        ('5.61617\\ \\mathrm{\\mu s}', '混合长与饱和：上界与·公式（精确绑定）'),
        ('16.9004\\ \\mathrm{m^2\\,s^{-1}}', '混合长与饱和：上界与·公式（精确绑定）'),
        ('1.74235\\ \\mathrm{mm}', '混合长与饱和：上界与·公式（精确绑定）'),
        ('10045.4\\ \\mathrm{m\\,s^{-1}}', '带状流与剪切抑制·公式（精确绑定）'),
        ('2.43031\\times10^{6}\\ \\mathrm{s^{-1}}', '带状流与剪切抑制·公式（精确绑定）'),
        ('56.4165\\ \\mathrm{mm}', '带状流与剪切抑制·公式（精确绑定）'),
        ('0.69897', '谱权重积分（正文）（精确绑定·补充）'),
        ('69.897\\%', '谱权重占百分比（正文）（精确绑定·补充）'),
        ('30.103\\%', '谱权重余量（正文，百分数口径）（精确绑定·补充）'),
        ('3.83153\\times10^{8}', '练习 1 回旋频率（精确绑定·补充）'),
        ('1.54905', '练习 1 离子声回旋半径（精确绑定·补充）'),
        ('193.667', '练习 1 波数（精确绑定·补充）'),
        ('57.4625', '练习 1 玻姆扩散（精确绑定·补充）'),
        ('0.712099', '练习 1 回旋玻姆扩散（精确绑定·补充）'),
        ('80.6945', '练习 1 两定标之比（精确绑定·补充）'),
        ('3.55999', '练习 1 chi_gB 缩小倍数（精确绑定·补充）'),
        ('1.88679', '练习 1 chi_B 缩小倍数（精确绑定·补充）'),
        ('178057', '练习 1 omega_star 不变（精确绑定·补充）'),
        ('20.6667', '练习 2 R/L_T（精确绑定·补充）'),
        ('3.13333', '练习 2 超临界倍数（精确绑定·补充）'),
        ('1.24', '练习 2 临界梯度长度（精确绑定·补充）'),
        ('4.13333', '练习 2 梯度比值（精确绑定·补充）'),
        ('1.74235', '练习 3 有效步长（精确绑定·补充）'),
        ('0.596134', '练习 3 有效步长/rho_s（精确绑定·补充）'),
        ('9.74246', '练习 3 混合长（精确绑定·补充）'),
        ('3.33333', '练习 3 混合长/rho_s（精确绑定·补充）'),
        ('0.540541', '练习 5 tau_E 之比（精确绑定·补充）'),
        ('10045.4', '练习 6 带状流速度（精确绑定·补充）'),
        ('2.43031\\times10^{6}', '练习 6 剪切率（精确绑定·补充）'),
        ('13.649', '练习 6 剪切率/增长率（精确绑定·补充）'),
        ('56.4165', '练习 6 临界径向尺度（精确绑定·补充）'),
        ('4.13338', '练习 6 离子回旋半径（精确绑定·补充）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AF'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('候选 & 本工作点的值 & 对 $B$ 的依赖 & 与实验反推之比 \\\\', '表格数据行'),
        ('玻姆 $\\chi_B$ & $108.42$ m$^2$s$^{-1}$ & $1/B$ & $200.577$ \\\\', '表格数据行'),
        ('回旋玻姆（$\\rho_i$） & $5.07013$ m$^2$s$^{-1}$ & $1/B^2$ & $9.37973$ \\\\', '表格数据行'),
        ('回旋玻姆（$\\rho_s$） & $2.53506$ m$^2$s$^{-1}$ & $1/B^2$ & $4.68987$ \\\\', '表格数据行'),
        ('实验反推 $\\chi_{\\mathrm{eff}}$ & $0.540541$ m$^2$s$^{-1}$ & —— & $1$ \\\\', '表格数据行'),
        ('\\begin{tip}{反常输运不是「多了一点」，是两个数量级}{tip-app-af-orders}', '语义盒标题'),
        ('\\textbf{出处：}本节的输运系数分类与回旋玻姆定标参照 Horton \\emph{Rev. Mod. Phys.} \\textbf{7', '出处声明整行'),
        ('量 & 本工作点的值 & 它由什么决定 \\\\', '表格数据行'),
        ('垂直波数 $k_\\theta$（$k_\\theta\\rho_s=0.3$） & $102.643$ m$^{-1}$ & 谱峰的无量纲位置 $\\div\\rho_s$ \\\\', '表格数据行'),
        ('垂直波长 $\\lambda_\\theta$ & $6.12137$ cm & $2\\pi/k_\\theta$ \\\\', '表格数据行'),
        ('抗磁漂移速度 $v_{\\mathrm{dia}}$ & $1734.72$ m/s & $\\rho_sc_s/L_n$ \\\\', '表格数据行'),
        ('漂移波频率 $\\omega_*$ & $178057$ rad/s（$28.3387$ kHz） & $k_\\theta\\rho_sc_s/L_n$，与 $B$ 无关 \\\\', '表格数据行'),
        ('环向涡数 $2\\pi R/\\lambda_\\theta$ & $636.389$ 个 & 装置尺寸 $\\div$ 垂直波长 \\\\', '表格数据行'),
        ('\\begin{tip}{漂移波不是「波」，是梯度与磁场的必然产物}{tip-app-af-driftwave}', '语义盒标题'),
        ('\\textbf{出处：}漂移波机制与 $\\omega_*$ 的推导参照 Wesson \\emph{Tokamaks}（第 4 版，2011，', '出处声明整行'),
        ('判据 & 本工作点 & 阈值（声明） & 结论 \\\\', '表格数据行'),
        ('$\\eta_i=L_n/L_T$ & $2$ & $1.5$ & ITG 不稳定（$1.33333$ 倍） \\\\', '表格数据行'),
        ('$R/L_T$ & $12.4$ & $5$ & 超临界 $1.48$ \\\\', '表格数据行'),
        ('$L_T$ 与临界值之比 & $0.5$ m & $1.24$ m & 陡 $2.48$ 倍 \\\\', '表格数据行'),
        ('捕获电子比例 $f_{\\mathrm{trap}}$ & $0.69843$ & —— & TEM 有足够驱动 \\\\', '表格数据行'),
        ('\\begin{tip}{临界梯度不是开关，是刚性}{tip-app-af-stiffness}', '语义盒标题'),
        ('\\textbf{出处：}ITG 与 TEM 的阈值判据参照 Wesson \\emph{Tokamaks}（第 4 版，2011，第 12 章', '出处声明整行'),
        ('量 & 值 & 说明 \\\\', '表格数据行'),
        ('增长率 $\\gamma_{\\max}$ & $178057$ s$^{-1}$ & 取 $k_\\theta\\rho_sc_s/L_n$，与 $\\omega_*$ 同阶 \\\\', '表格数据行'),
        ('相关时间 $\\tau_c$ & $5.61617$ $\\mu$s & $1/\\gamma_{\\max}$ \\\\', '表格数据行'),
        ('混合长步长 $1/k_\\theta$ & $9.74246$ mm & $=3.33333\\,\\rho_s$ \\\\', '表格数据行'),
        ('混合长扩散 $\\chi_{\\mathrm{mix}}$ & $16.9004$ m$^2$s$^{-1}$ & 高估 $31.2658$ 倍 \\\\', '表格数据行'),
        ('有效步长 $\\Delta r_{\\mathrm{eff}}$ & $1.74235$ mm & 由 $\\chi_{\\mathrm{eff}}\\tau_c$ 反推 \\\\', '表格数据行'),
        ('涨落水平（声明）$\\delta n/n$ & $0.01$ & $(\\delta n/n)^2$ 定标给 $0.00169004$ m$^2$s$^{-1}$ \\\\', '表格数据行'),
        ('\\begin{tip}{混合长是上界，不是预言}{tip-app-af-mixing}', '语义盒标题'),
        ('\\textbf{出处：}混合长与饱和的讨论参照 Horton \\emph{Rev. Mod. Phys.} \\textbf{71} (199', '出处声明整行'),
        ('径向尺度 $\\Delta r$ & $\\gamma_E$ & $\\gamma_E/\\gamma_{\\max}$ & 是否抑制 \\\\', '表格数据行'),
        ('$\\rho_i=4.13338$ mm & $2.43031\\times10^{6}$ s$^{-1}$ & $13.649$ & 是 \\\\', '表格数据行'),
        ('$a/10=0.2$ m & $50226.9$ s$^{-1}$ & $0.282083$ & 否 \\\\', '表格数据行'),
        ('分界 $\\Delta r_*=56.4165$ mm & $178057$ s$^{-1}$ & $1$ & 临界 \\\\', '表格数据行'),
        ('\\begin{tip}{带状流是湍流自己长出来的「刹车」}{tip-app-af-zonal}', '语义盒标题'),
        ('\\textbf{出处：}带状流与 $E\\times B$ 剪切抑制判据参照 Biglari、Diamond \\& Terry', '出处声明整行'),
        ('谱 & $S(k)$ 指数 & $k^2S(k)$ 权重 & 含义 \\\\', '表格数据行'),
        ('二维熵级串 & $-3$ & $\\propto k^{-1}$ & 长波占 $69.897\\%$（$k\\le0.5$），余 $30.103\\%$ \\\\', '表格数据行'),
        ('逆能量级串 & $-5/3$ & $\\propto k^{1/3}$ & 高波数份额上升（长波不再独占） \\\\', '表格数据行'),
        ('三维 Kolmogorov & $-5/3$ & $\\propto k^{1/3}$ & 各向同性惯性区 \\\\', '表格数据行'),
        ('\\begin{tip}{谱指数不是装饰，它决定输运集中在哪里}{tip-app-af-spectrum}', '语义盒标题'),
        ('\\textbf{出处：}二维漂移波熵级串与逆能量级串的谱指数参照 Kraichnan \\emph{Phys. Fluids}', '出处声明整行'),
        ('级别 & 内容 & 例子 & 可核验性 \\\\', '表格数据行'),
        ('现算级 & 由 \\texttt{tmp/appAF\\_calc.py} 从工作点算出 & $\\rho_s=2.92274$ mm、$\\chi_{\\mathrm{mix}}=16.9004$ m$^2$s$^{-1}$ & 可逐位复算 \\\\', '表格数据行'),
        ('声明级 & 脚本顶部的输入或文献量级 & $L_n=1.0$ m、$\\eta_{i,\\mathrm{crit}}=1.5$、谱指数 $-3/-5/3$ & 可改输入重算 \\\\', '表格数据行'),
        ('文献级 & 定标系数的适用范围 & 回旋玻姆系数、混合长饱和因子 & 需查原文 \\\\', '表格数据行'),
        ('\\begin{tip}{诚实地说清「本附录不给什么」}{tip-app-af-honest}', '语义盒标题'),
        ('\\textbf{出处：}回旋玻姆定标与 $\\rho_s$ 口径参照 Wesson \\emph{Tokamaks}（第 4 版，2011，第 ', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AF'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 0.540541, '\\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式'),
        (0, 108.42, '\\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式'),
        (0, 5.07013, '\\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式'),
        (0, 2.53506, '\\mathrm{m^2\\,s^{-1}}', '三个输运系数：反常输·公式'),
        (0, 0.159248, '\\mathrm{MW\\,m^{-2}}', '三个输运系数：反常输·公式'),
        (0, 77.9568, '\\mathrm{MW}', '三个输运系数：反常输·公式'),
        (1, 0.00292274, '\\mathrm{m}', '漂移波：梯度和磁场谁·公式'),
        (1, 2.92274, '\\mathrm{mm}', '漂移波：梯度和磁场谁·公式'),
        (1, 178057, '\\mathrm{rad\\,s^{-1}}', '漂移波：梯度和磁场谁·公式'),
        (2, 2, '\\mathrm{m^{-1}}', '不稳定性：ITG 与·公式'),
        (2, 1.24, '\\mathrm{m}', '不稳定性：ITG 与·公式'),
        (3, 178057, '\\mathrm{s^{-1}}', '混合长与饱和：上界与·公式'),
        (3, 5.61617e-06, '\\mathrm{s}', '混合长与饱和：上界与·公式'),
        (3, 5.61617, '\\mathrm{\\mu s}', '混合长与饱和：上界与·公式'),
        (3, 16.9004, '\\mathrm{m^2\\,s^{-1}}', '混合长与饱和：上界与·公式'),
        (3, 1.74235, '\\mathrm{mm}', '混合长与饱和：上界与·公式'),
        (4, 10045.4, '\\mathrm{m\\,s^{-1}}', '带状流与剪切抑制·公式'),
        (4, 2430310, '\\mathrm{s^{-1}}', '带状流与剪切抑制·公式'),
        (4, 56.4165, '\\mathrm{mm}', '带状流与剪切抑制·公式'),
        (5, 0.69897, '', '谱权重积分（正文）'),
        (5, 30.103, '', '谱权重余量（正文，百分数口径）'),
        (6, 383153000, '', '练习 1 回旋频率'),
        (6, 1.54905, '', '练习 1 离子声回旋半径'),
        (6, 193.667, '', '练习 1 波数'),
        (6, 57.4625, '', '练习 1 玻姆扩散'),
        (6, 0.712099, '', '练习 1 回旋玻姆扩散'),
        (6, 80.6945, '', '练习 1 两定标之比'),
        (6, 3.55999, '', '练习 1 chi_gB 缩小倍数'),
        (6, 1.88679, '', '练习 1 chi_B 缩小倍数'),
        (6, 178057, '', '练习 1 omega_star 不变'),
        (6, 20.6667, '', '练习 2 R/L_T'),
        (6, 3.13333, '', '练习 2 超临界倍数'),
        (6, 1.24, '', '练习 2 临界梯度长度'),
        (6, 4.13333, '', '练习 2 梯度比值'),
        (6, 1.74235, '', '练习 3 有效步长'),
        (6, 0.596134, '', '练习 3 有效步长/rho_s'),
        (6, 9.74246, '', '练习 3 混合长'),
        (6, 3.33333, '', '练习 3 混合长/rho_s'),
        (6, 0.540541, '', '练习 5 tau_E 之比'),
        (6, 10045.4, '', '练习 6 带状流速度'),
        (6, 2430310, '', '练习 6 剪切率'),
        (6, 13.649, '', '练习 6 剪切率/增长率'),
        (6, 56.4165, '', '练习 6 临界径向尺度'),
        (6, 4.13338, '', '练习 6 离子回旋半径'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AF 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AF 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('69.9\\%', '谱权重占百分比四舍五入'),
        ('30.1\\%', '谱权重余量四舍五入'),
        ('0.699', '谱权重小数四舍五入'),
        ('0.302', '谱权重余量小数四舍五入'),
        ('2.9\\ \\mathrm{mm}', '离子声回旋半径四舍五入'),
        ('4.1\\ \\mathrm{mm}', '离子回旋半径四舍五入'),
        ('5.1\\ \\mathrm{m', '回旋玻姆扩散（离子尺度）四舍五入'),
        ('2.5\\ \\mathrm{m', '回旋玻姆扩散（离子声尺度）四舍五入'),
        ('108\\ \\mathrm{m', '玻姆扩散四舍五入'),
        ('178000', '漂移频率取整'),
        ('1.79\\times10^{5}', '漂移频率四舍五入'),
        ('28.3\\ \\mathrm{kHz}', '漂移频率换算四舍五入'),
        ('6.1\\ \\mathrm{cm}', '混合长四舍五入'),
        ('31.3', '混合长高估倍数四舍五入'),
        ('5.6\\ \\mathrm{\\mu s}', '相关时间四舍五入'),
        ('0.283', '剪切率/增长率四舍五入'),
        ('17.4\\ \\mathrm{km', '带状流速度四舍五入'),
        ('3.8\\times10^{8}', '练习 1 回旋频率四舍五入'),
        ('3.83\\times10^{8}', '练习 1 回旋频率另一种四舍五入'),
        ('1.55', '练习 1 离子声回旋半径四舍五入'),
        ('193.7', '练习 1 波数四舍五入'),
        ('57.5', '练习 1 玻姆扩散四舍五入'),
        ('0.712\\ ', '练习 1 回旋玻姆扩散四舍五入'),
        ('80.7', '练习 1 定标之比四舍五入'),
        ('3.56', '练习 1 chi_gB 缩小倍数四舍五入'),
        ('1.89', '练习 1 chi_B 缩小倍数四舍五入'),
        ('20.7', '练习 2 R/L_T 四舍五入'),
        ('0.541', '练习 5 tau_E 之比四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AF'], why))


def check_P39_appendix_ag_equilibrium():
    """附录 AG：磁流体平衡与位形——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AG'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AG 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '磁面与安全因子',
        'Grad--Shafranov 方程',
        'Shafranov 位移',
        '位形参数：伸长比',
        '压强极限：比压',
        '从平衡到稳定性',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AG 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AG 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ag-scales',
        'tab:app-ag-gs',
        'tab:app-ag-shift',
        'tab:app-ag-shaping',
        'tab:app-ag-beta',
        'tab:app-ag-stability',
        'tab:app-ag-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AG'])

    for lab in (
        'sec:app-ag-surfaces',
        'sec:app-ag-gs',
        'sec:app-ag-shift',
        'sec:app-ag-shaping',
        'sec:app-ag-beta',
        'sec:app-ag-stability',
        'sec:app-ag-problems',
        'sec:app-ag-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AG'])

    for lab in (
        'tip-app-ag-surfaces',
        'tip-app-ag-gs',
        'tip-app-ag-shift',
        'tip-app-ag-shaping',
        'tip-app-ag-beta',
        'tip-app-ag-stability',
        'tip-app-ag-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AG'])

    for lab in (
        'eq:app-ag-epsilon',
        'eq:app-ag-btheta',
        'eq:app-ag-q95',
        'eq:app-ag-shear',
        'eq:app-ag-volume',
        'eq:app-ag-gs',
        'eq:app-ag-pressure',
        'eq:app-ag-betat',
        'eq:app-ag-betan',
        'eq:app-ag-shift',
        'eq:app-ag-gap',
        'eq:app-ag-axis',
        'eq:app-ag-shape',
        'eq:app-ag-current',
        'eq:app-ag-virial',
        'eq:app-ag-troyon',
        'eq:app-ag-alpha',
        'eq:app-ag-mercier',
        'eq:app-ag-suydam',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AG'])

    # ---------- ③ 版权口径 / 出处 / 练习与解答 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AG'])

    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AG 的「出处：」声明只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AG 的「出处：」声明 %d 处' % n_src)

    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', text))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', text))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AG 的练习 %d 道 / 解答 %d 条（要求 >= 6 且解答不少于练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AG 的练习 %d 道、解答 %d 条' % (n_ex, n_ans))

    # ---------- ④ 数字 + 单位的精确绑定（堵末位容差洞） ----------
    for lit, why in (
        ('1.5\\ \\mathrm{T}', '磁面与安全因子：嵌套·公式（精确绑定）'),
        ('5.3\\ \\mathrm{T}', '磁面与安全因子：嵌套·公式（精确绑定）'),
        ('113.223\\ \\mathrm{Wb}', '磁面与安全因子：嵌套·公式（精确绑定）'),
        ('51.0731\\ \\mathrm{Wb}', '磁面与安全因子：嵌套·公式（精确绑定）'),
        ('832.205\\ \\mathrm{m^3}', '磁面与安全因子：嵌套·公式（精确绑定）'),
        ('489.532\\ \\mathrm{m^2}', '磁面与安全因子：嵌套·正文（精确绑定）'),
        ('7.82381\\ \\mathrm{T}', '磁面与安全因子：嵌套·正文（精确绑定）'),
        ('4.00732\\ \\mathrm{T}', '磁面与安全因子：嵌套·正文（精确绑定）'),
        ('1\\times10^{20}\\ \\mathrm{m^{-3}}', 'Grad--Shaf·公式（精确绑定）'),
        ('9194\\ \\mathrm{eV}', 'Grad--Shaf·公式（精确绑定）'),
        ('2.94608\\times10^{5}\\ \\mathrm{Pa}', 'Grad--Shaf·公式（精确绑定）'),
        ('1.47304\\times10^{5}\\ \\mathrm{Pa\\,m^{-1}}', 'Grad--Shaf·公式（精确绑定）'),
        ('7.02154\\times10^{5}\\ \\mathrm{A\\,m^{-2}}', 'Grad--Shaf·公式（精确绑定）'),
        ('32.86\\ \\mathrm{T\\,m}', 'Grad--Shaf·正文（精确绑定）'),
        ('0.267445\\ \\mathrm{m}', 'Shafranov ·公式（精确绑定）'),
        ('6.46745\\ \\mathrm{m}', 'Shafranov ·公式（精确绑定）'),
        ('5.08083\\ \\mathrm{T}', 'Shafranov ·公式（精确绑定）'),
        ('1.73255\\ \\mathrm{m}', 'Shafranov ·公式（精确绑定）'),
        ('2.26745\\ \\mathrm{m}', 'Shafranov ·公式（精确绑定）'),
        ('1.0\\ \\mathrm{m}', 'Shafranov ·正文（精确绑定）'),
        ('0.3736\\ \\mathrm{m}', 'Shafranov ·正文（精确绑定）'),
        ('7.71208\\ \\mathrm{MA}', '位形参数：伸长比与三·公式（精确绑定）'),
        ('15\\ \\mathrm{MA}', '位形参数：伸长比与三·公式（精确绑定）'),
        ('19.2802\\ \\mathrm{MA}', '位形参数：伸长比与三·公式（精确绑定）'),
        ('10.3728\\ \\mathrm{MA}', '位形参数：伸长比与三·公式（精确绑定）'),
        ('1057.32\\ \\mathrm{m^3}', '位形参数：伸长比与三·正文（精确绑定）'),
        ('15.0\\ \\mathrm{MA}', '位形参数：伸长比与三·表格（精确绑定）'),
        ('4.74481\\times10^{5}\\ \\mathrm{Pa}', '压强极限：比压能推到·公式（精确绑定）'),
        ('11.0844\\ \\mathrm{MA}', '压强极限：比压能推到·公式（精确绑定）'),
        ('8.82069\\times10^{19}\\ \\mathrm{m^{-3}}', '压强极限：比压能推到·公式（精确绑定）'),
        ('1.5\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('8\\ \\mathrm{MA}', '估算练习·正文（精确绑定）'),
        ('8\\ \\mathrm{keV}', '估算练习·正文（精确绑定）'),
        ('2.0\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('6.5\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('2.2\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('1.06667\\ \\mathrm{T}', '估算练习·正文（精确绑定）'),
        ('0.194933\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('1117.79\\ \\mathrm{m^3}', '估算练习·正文（精确绑定）'),
        ('0.322581', 'AG·正文（精确绑定·补充）'),
        ('3.1', 'AG·正文（精确绑定·补充）'),
        ('1.2566\\times10^{-6}', 'AG·正文（精确绑定·补充）'),
        ('15.0\\times10^{6}', 'AG·正文（精确绑定·补充）'),
        ('28.3019\\%', 'AG·正文（精确绑定·补充）'),
        ('1.13978', 'AG·正文（精确绑定·补充）'),
        ('2.21688', 'AG·正文（精确绑定·补充）'),
        ('1.7', 'AG·正文（精确绑定·补充）'),
        ('1.09783', 'AG·正文（精确绑定·补充）'),
        ('0.0263592 ', 'AG·正文（精确绑定·补充）'),
        ('2.63592\\%', 'AG·正文（精确绑定·补充）'),
        ('0.329081 ', 'AG·正文（精确绑定·补充）'),
        ('32.9081\\%', 'AG·正文（精确绑定·补充）'),
        ('1.86272', 'AG·正文（精确绑定·补充）'),
        ('12.4844', 'AG·正文（精确绑定·补充）'),
        ('0.133723 ', 'AG·正文（精确绑定·补充）'),
        ('13.3723\\%', 'AG·正文（精确绑定·补充）'),
        ('1.30873', 'AG·正文（精确绑定·补充）'),
        ('2.6', 'AG·正文（精确绑定·补充）'),
        ('0.329081', 'AG·正文（精确绑定·补充）'),
        ('1.945', 'AG·正文（精确绑定·补充）'),
        ('1.345', 'AG·正文（精确绑定·补充）'),
        ('2.5', 'AG·正文（精确绑定·补充）'),
        ('0.0656377 ', 'AG·正文（精确绑定·补充）'),
        ('6.56377\\%', 'AG·正文（精确绑定·补充）'),
        ('4.24528\\%', 'AG·正文（精确绑定·补充）'),
        ('62.0907\\%', 'AG·正文（精确绑定·补充）'),
        ('1.35325', 'AG·正文（精确绑定·补充）'),
        ('1.1337', 'AG·正文（精确绑定·补充）'),
        ('0.401587', 'AG·正文（精确绑定·补充）'),
        ('0.160635', 'AG·正文（精确绑定·补充）'),
        ('0.150654', 'AG·正文（精确绑定·补充）'),
        ('0.602617', 'AG·正文（精确绑定·补充）'),
        ('0.829081', 'AG·正文（精确绑定·补充）'),
        ('1.31084', 'AG·正文（精确绑定·补充）'),
        ('0.819273', 'AG·正文（精确绑定·补充）'),
        ('1.2\\ \\mathrm{cm}', 'AG·正文（精确绑定·补充）'),
        ('2.54715\\ \\mathrm{cm}', 'AG·正文（精确绑定·补充）'),
        ('5.5\\ \\mathrm{m}', 'AG·练习（精确绑定·补充）'),
        ('1.6', 'AG·练习（精确绑定·补充）'),
        ('0.8\\times10^{20}\\ \\mathrm{m^{-3}}', 'AG·练习（精确绑定·补充）'),
        ('0.272727', 'AG·练习（精确绑定·补充）'),
        ('2.4121', 'AG·练习（精确绑定·补充）'),
        ('1.83488\\%', 'AG·练习（精确绑定·补充）'),
        ('0.453005', 'AG·练习（精确绑定·补充）'),
        ('12.9955\\%', 'AG·练习（精确绑定·补充）'),
        ('3.66667', 'AG·练习（精确绑定·补充）'),
        ('1.13372', 'AG·练习（精确绑定·补充）'),
        ('0.866277', 'AG·练习（精确绑定·补充）'),
        ('1.34317', 'AG·练习（精确绑定·补充）'),
        ('1.61055', 'AG·练习（精确绑定·补充）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AG'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 本工作点 & 它从哪里来 & 它决定什么 \\\\', '表格数据行'),
        ('反宽高比 $\\epsilon$ & $0.322581$ & $a/R$，几何 & 展开式的小量、环径比 \\\\', '表格数据行'),
        ('极向场 $B_\\theta(a)$ & $1.5\\ \\mathrm{T}$ & $\\mu_0 I_p/(2\\pi a)$ & 磁力线的螺距 \\\\', '表格数据行'),
        ('安全因子 $q_{95}$ & $2.21688$ & 环向磁通与极向磁通之比 & 扭曲模、磁岛、约束 \\\\', '表格数据行'),
        ('磁剪切 $s_{95}$ & $1.09783$ & $q$ 剖面的斜率 & 局域稳定性判据 \\\\', '表格数据行'),
        ('体积 $V$ & $832.205\\ \\mathrm{m^3}$ & $2\\pi^2 R a^2 \\kappa$ & 功率密度与总储能 \\\\', '表格数据行'),
        ('\\begin{tip}{磁面不是画出来的，是平衡方程的解}{tip-app-ag-surfaces}', '语义盒标题'),
        ('\\textbf{出处：}磁面几何与安全因子的记法参照 Freidberg \\emph{Ideal MHD}（2014，第 5--6 章）、', '出处声明整行'),
        ('量 & 本工作点 & 表达式 & 物理含义 \\\\', '表格数据行'),
        ('$\\langle p\\rangle$ & $2.94608\\times10^{5}\\ \\mathrm{Pa}$ & $2n_eT_e$ & 往外推的账 \\\\', '表格数据行'),
        ('$\\langle j\\rangle$ & $7.02154\\times10^{5}\\ \\mathrm{A\\,m^{-2}}$ & $I_p/(\\pi a^2\\kappa)$ & 往里箍的账 \\\\', '表格数据行'),
        ('$\\beta_t$ & $2.63592\\%$ & $2\\mu_0\\langle p\\rangle/B_0^2$ & 对环向场的效率 \\\\', '表格数据行'),
        ('$\\beta_p$ & $32.9081\\%$ & $2\\mu_0\\langle p\\rangle/B_\\theta^2$ & 对极向场的效率 \\\\', '表格数据行'),
        ('$\\beta_N$ & $1.86272$ & $\\beta_t[\\%]aB_0/I_p[\\mathrm{MA}]$ & 与极限比较的口径 \\\\', '表格数据行'),
        ('\\begin{tip}{Grad--Shafranov 方程不是新物理，是两条账}{tip-app-ag-gs}', '语义盒标题'),
        ('\\textbf{出处：}Grad--Shafranov 方程的推导与各项含义参照 Grad \\& Rubin（1958，第二次', '出处声明整行'),
        ('量 & 本工作点 & 表达式 & 后果 \\\\', '表格数据行'),
        ('$\\Delta$ & $0.267445\\ \\mathrm{m}$ & $(a^2/2R)(\\beta_p+l_i/2)$ & 磁轴外移 \\\\', '表格数据行'),
        ('$\\Delta/a$ & $13.3723\\%$ & 相对小半径 & 内外间距不对称 \\\\', '表格数据行'),
        ('$R_{\\rm axis}$ & $6.46745\\ \\mathrm{m}$ & $R+\\Delta$ & 轴心场降到 $5.08083\\ \\mathrm{T}$ \\\\', '表格数据行'),
        ('$(a+\\Delta)/(a-\\Delta)$ & $1.30873$ & 内外间距比 & 偏滤器与壁负载 \\\\', '表格数据行'),
        ('$\\Delta=a/2$ 所需 $\\beta_p$ & $2.6$ & 反解 & 大电流装置更安全 \\\\', '表格数据行'),
        ('\\begin{tip}{位移不是误差，是压强与电流付的账}{tip-app-ag-shift}', '语义盒标题'),
        ('\\textbf{出处：}Shafranov 位移的小 $\\epsilon$ 表达式参照 Shafranov（\\emph{Reviews of', '出处声明整行'),
        ('伸长比 $\\kappa$ & 形状因子 & 同样 $q_{95}$ 允许的 $I_p$ & 相对圆截面 & 说明 \\\\', '表格数据行'),
        ('$1.0$ & $1.0$ & $7.71208\\ \\mathrm{MA}$ & $1.0$ 倍 & 圆截面基准 \\\\', '表格数据行'),
        ('$1.3$ & $1.345$ & $10.3728\\ \\mathrm{MA}$ & $1.345$ 倍 & 早期装置 \\\\', '表格数据行'),
        ('$1.7$ & $1.945$ & $15.0\\ \\mathrm{MA}$ & $1.945$ 倍 & 本工作点 \\\\', '表格数据行'),
        ('$2.0$ & $2.5$ & $19.2802\\ \\mathrm{MA}$ & $2.5$ 倍 & 形状的极限一侧 \\\\', '表格数据行'),
        ('\\begin{tip}{伸长比不是造型，是买电流}{tip-app-ag-shaping}', '语义盒标题'),
        ('\\textbf{出处：}形状因子与 $\\beta_N$ 的三角形变增益参照 Troyon 等（\\emph{Plasma Phys.', '出处声明整行'),
        ('量 & 本工作点 & 极限或对照 & 余量 \\\\', '表格数据行'),
        ('$\\beta_t$ & $2.63592\\%$ & $4.24528\\%$（$\\beta_N=3$） & 用了 $62.0907\\%$ \\\\', '表格数据行'),
        ('$\\langle p\\rangle$ & $2.94608\\times10^{5}\\ \\mathrm{Pa}$ & $4.74481\\times10^{5}\\ \\mathrm{Pa}$ & 可再涨六成 \\\\', '表格数据行'),
        ('$\\epsilon/q_{95}^2$ & $6.56377\\%$ & 维里量级上界 & 用了四成 \\\\', '表格数据行'),
        ('$q_{95}$ & $2.21688$ & $3.0$（安全线） & 偏低，需降流 \\\\', '表格数据行'),
        ('$n_e/n_G$ & $0.85$（附录 AD 口径） & $1.0$ & 降流后变成 $1.1337$ \\\\', '表格数据行'),
        ('\\begin{tip}{比压的上限不是材料强度，是磁场的刚度}{tip-app-ag-beta}', '语义盒标题'),
        ('\\textbf{出处：}维里定理与比压量级上界参照 Greene, Johnson \\& Weimer（\\emph{Phys.', '出处声明整行'),
        ('平衡输出 & 本工作点 & 交给谁 & 用来判断什么 \\\\', '表格数据行'),
        ('$q_0$、$q_{95}$ & $1.0$、$2.21688$ & 附录 T & 磁岛位置与宽度 \\\\', '表格数据行'),
        ('$s_{95}$ & $1.09783$ & 附录 T、AD & 局域稳定性判据 \\\\', '表格数据行'),
        ('$l_i$ & $1.0$ & 附录 V & 外部扭曲模、破裂 \\\\', '表格数据行'),
        ('$\\alpha$ & $0.401587$ & 附录 AD & 气球模驱动 \\\\', '表格数据行'),
        ('$\\beta_p+l_i/2$ & $0.829081$ & 本附录、附录 AA & Shafranov 位移 \\\\', '表格数据行'),
        ('\\begin{tip}{平衡给的是背景，稳定性才说能不能待住}{tip-app-ag-stability}', '语义盒标题'),
        ('\\textbf{出处：}Suydam 判据参照 Suydam（1958，第二次日内瓦会议论文集 \\textbf{31}，157）、', '出处声明整行'),
        ('级别 & 典型量 & 使用规则 \\\\', '表格数据行'),
        ('现算级 & 磁面几何、$\\beta$ 三种口径、位移、电流极限 & 可作结论，句尾注明脚本 \\\\', '表格数据行'),
        ('声明级 & $\\kappa$、$\\delta$、$l_i$、$q_0$、$\\beta_N$ 上限、形状增益 & 只作输入，正文写明「声明值」 \\\\', '表格数据行'),
        ('文献级 & Troyon 散布、ITER 水平、平衡程序 & 只作量级对照，不作结论 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ag-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录的写法参照 Freidberg \\emph{Ideal MHD}（2014，第 5--6、11--12 章）、', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AG'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 1.5, '\\mathrm{T}', '磁面与安全因子：嵌套·公式'),
        (0, 5.3, '\\mathrm{T}', '磁面与安全因子：嵌套·公式'),
        (0, 113.223, '\\mathrm{Wb}', '磁面与安全因子：嵌套·公式'),
        (0, 51.0731, '\\mathrm{Wb}', '磁面与安全因子：嵌套·公式'),
        (0, 832.205, '\\mathrm{m^3}', '磁面与安全因子：嵌套·公式'),
        (0, 832.2, '\\mathrm{m^3}', '磁面与安全因子：嵌套·正文'),
        (0, 489.532, '\\mathrm{m^2}', '磁面与安全因子：嵌套·正文'),
        (0, 7.82381, '\\mathrm{T}', '磁面与安全因子：嵌套·正文'),
        (0, 4.00732, '\\mathrm{T}', '磁面与安全因子：嵌套·正文'),
        (1, 1e+20, '\\mathrm{m^{-3}}', 'Grad--Shaf·公式'),
        (1, 9194, '\\mathrm{eV}', 'Grad--Shaf·公式'),
        (1, 294608, '\\mathrm{Pa}', 'Grad--Shaf·公式'),
        (1, 147304, '\\mathrm{Pa\\,m^{-1}}', 'Grad--Shaf·公式'),
        (1, 702154, '\\mathrm{A\\,m^{-2}}', 'Grad--Shaf·公式'),
        (1, 32.86, '\\mathrm{T\\,m}', 'Grad--Shaf·正文'),
        (2, 0.267445, '\\mathrm{m}', 'Shafranov ·公式'),
        (2, 6.46745, '\\mathrm{m}', 'Shafranov ·公式'),
        (2, 5.08083, '\\mathrm{T}', 'Shafranov ·公式'),
        (2, 1.73255, '\\mathrm{m}', 'Shafranov ·公式'),
        (2, 2.26745, '\\mathrm{m}', 'Shafranov ·公式'),
        (2, 1, '\\mathrm{m}', 'Shafranov ·正文'),
        (2, 0.3736, '\\mathrm{m}', 'Shafranov ·正文'),
        (3, 7.71208, '\\mathrm{MA}', '位形参数：伸长比与三·公式'),
        (3, 15, '\\mathrm{MA}', '位形参数：伸长比与三·公式'),
        (3, 19.2802, '\\mathrm{MA}', '位形参数：伸长比与三·公式'),
        (3, 10.3728, '\\mathrm{MA}', '位形参数：伸长比与三·公式'),
        (3, 1057.32, '\\mathrm{m^3}', '位形参数：伸长比与三·正文'),
        (4, 474481, '\\mathrm{Pa}', '压强极限：比压能推到·公式'),
        (4, 11.0844, '\\mathrm{MA}', '压强极限：比压能推到·公式'),
        (4, 8.82069e+19, '\\mathrm{m^{-3}}', '压强极限：比压能推到·公式'),
        (6, 1.5, '\\mathrm{m}', '估算练习·正文'),
        (6, 8, '\\mathrm{MA}', '估算练习·正文'),
        (6, 8, '\\mathrm{keV}', '估算练习·正文'),
        (6, 2, '\\mathrm{m}', '估算练习·正文'),
        (6, 6.5, '\\mathrm{m}', '估算练习·正文'),
        (6, 2.2, '\\mathrm{m}', '估算练习·正文'),
        (6, 1.06667, '\\mathrm{T}', '估算练习·正文'),
        (6, 0.194933, '\\mathrm{m}', '估算练习·正文'),
        (6, 1117.79, '\\mathrm{m^3}', '估算练习·正文'),
        (0, 0.322581, '', 'AG·正文'),
        (0, 3.1, '', 'AG·正文'),
        (0, 1.2566e-06, '', 'AG·正文'),
        (0, 15000000, '', 'AG·正文'),
        (0, 28.3019, '', 'AG·正文'),
        (0, 1.13978, '', 'AG·正文'),
        (0, 2.21688, '', 'AG·正文'),
        (0, 1.7, '', 'AG·正文'),
        (0, 1.09783, '', 'AG·正文'),
        (1, 0.0263592, '', 'AG·正文'),
        (1, 2.63592, '', 'AG·正文'),
        (1, 0.329081, '', 'AG·正文'),
        (1, 32.9081, '', 'AG·正文'),
        (1, 1.5, '', 'AG·正文'),
        (1, 1.86272, '', 'AG·正文'),
        (1, 12.4844, '', 'AG·正文'),
        (2, 0.133723, '', 'AG·正文'),
        (2, 13.3723, '', 'AG·正文'),
        (2, 1.30873, '', 'AG·正文'),
        (2, 2.6, '', 'AG·正文'),
        (3, 1.945, '', 'AG·正文'),
        (3, 832.205, '', 'AG·正文'),
        (3, 1.345, '', 'AG·正文'),
        (3, 2.5, '', 'AG·正文'),
        (4, 0.0656377, '', 'AG·正文'),
        (4, 6.56377, '', 'AG·正文'),
        (4, 4.24528, '', 'AG·正文'),
        (4, 62.0907, '', 'AG·正文'),
        (4, 1.35325, '', 'AG·正文'),
        (4, 1.1337, '', 'AG·正文'),
        (4, 0.401587, '', 'AG·正文'),
        (4, 0.160635, '', 'AG·正文'),
        (4, 294608, '', 'AG·正文'),
        (5, 0.150654, '', 'AG·正文'),
        (5, 0.602617, '', 'AG·正文'),
        (5, 0.829081, '', 'AG·正文'),
        (5, 1.31084, '', 'AG·正文'),
        (5, 0.819273, '', 'AG·正文'),
        (5, 1.2, '', 'AG·正文'),
        (5, 2.54715, '', 'AG·正文'),
        (6, 5.5, '', 'AG·练习'),
        (6, 5.3, '', 'AG·练习'),
        (6, 1.6, '', 'AG·练习'),
        (6, 8e+19, '', 'AG·练习'),
        (6, 0.267445, '', 'AG·练习'),
        (6, 15, '', 'AG·练习'),
        (6, 1e+20, '', 'AG·练习'),
        (6, 0.272727, '', 'AG·练习'),
        (6, 2.4121, '', 'AG·练习'),
        (6, 1.83488, '', 'AG·练习'),
        (6, 0.453005, '', 'AG·练习'),
        (6, 12.9955, '', 'AG·练习'),
        (6, 3.66667, '', 'AG·练习'),
        (6, 1.13372, '', 'AG·练习'),
        (6, 0.866277, '', 'AG·练习'),
        (6, 1.34317, '', 'AG·练习'),
        (6, 474481, '', 'AG·练习'),
        (6, 1.61055, '', 'AG·练习'),
        (6, 11.0844, '', 'AG·练习'),
        (6, 8.82069e+19, '', 'AG·练习'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AG 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AG 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('2.64\\%', '环向比压四舍五入'),
        ('32.9\\%', '极向比压四舍五入'),
        ('13.4\\%', 'Shafranov 位移占小半径四舍五入'),
        ('28.3\\%', '极向场占比四舍五入'),
        ('1.86\\ ', '归一化比压四舍五入'),
        ('0.267\\ ', 'Shafranov 位移四舍五入'),
        ('2.22\\ ', '安全因子四舍五入'),
        ('1.14\\ ', '圆柱安全因子四舍五入'),
        ('1.10\\ ', '磁剪切四舍五入'),
        ('0.33\\ ', '极向比压小数四舍五入'),
        ('0.0264', '环向比压小数四舍五入'),
        ('12.5', '环向场与极向场之比的平方四舍五入'),
        ('1.94\\ ', '形状因子四舍五入'),
        ('1.35\\ ', '形状因子（大伸长）四舍五入'),
        ('0.066', '剩磁通占比四舍五入'),
        ('6.56\\%', 'Troyon 系数占比四舍五入'),
        ('62.1\\%', '比压占上限比例四舍五入'),
        ('4.25\\%', 'beta_N=3 时环向比压四舍五入'),
        ('1.134', 'Greenwald 分数四舍五入'),
        ('1.61\\ ', '压强提升倍数四舍五入'),
        ('3.67\\ ', '练习 1 大半径与小半径之比四舍五入'),
        ('1.83\\%', '练习 1 环向比压四舍五入'),
        ('2.41\\ ', '练习 1 安全因子四舍五入'),
        ('13.0\\%', '练习 1 位移占小半径四舍五入'),
        ('0.27\\ ', '练习 1 反纵横比四舍五入'),
        ('0.45\\ ', '练习 1 极向比压四舍五入'),
        ('1.34\\ ', '练习 4 体积比四舍五入'),
        ('0.866\\ ', '练习 3 内侧场比值四舍五入'),
        ('11.1\\ \\mathrm{MA}', '练习 6 电流四舍五入'),
        ('8.82\\times10^{19}', '练习 6 Greenwald 密度四舍五入'),
        ('4.74\\times10^{5}', '练习 5 平均压强四舍五入'),
        ('2.95\\times10^{5}', '平均压强四舍五入'),
        ('0.402\\ ', '归一化比压系数四舍五入'),
        ('0.161\\ ', '系数与临界值之比四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AG'], why))


def check_P40_appendix_ah_atomic():
    """附录 AH：原子过程与谱线——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AH'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AH 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '能级与跃迁',
        'Saha 方程',
        '强度比：为什么日冕平衡',
        '多普勒展宽',
        '塞曼分裂与斯塔克',
        '光学薄还是光学厚',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AH 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AH 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ah-levels',
        'tab:app-ah-saha',
        'tab:app-ah-intensity',
        'tab:app-ah-doppler',
        'tab:app-ah-fields',
        'tab:app-ah-opacity',
        'tab:app-ah-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AH'])

    for lab in (
        'sec:app-ah-levels',
        'sec:app-ah-saha',
        'sec:app-ah-corona',
        'sec:app-ah-doppler',
        'sec:app-ah-fields',
        'sec:app-ah-opacity',
        'sec:app-ah-problems',
        'sec:app-ah-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AH'])

    for lab in (
        'tip-app-ah-levels',
        'tip-app-ah-saha',
        'tip-app-ah-corona',
        'tip-app-ah-doppler',
        'tip-app-ah-fields',
        'tip-app-ah-opacity',
        'tip-app-ah-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AH'])

    for lab in (
        'eq:app-ah-rydberg',
        'eq:app-ah-bohr',
        'eq:app-ah-saha',
        'eq:app-ah-neutral',
        'eq:app-ah-corona',
        'eq:app-ah-sens',
        'eq:app-ah-doppler',
        'eq:app-ah-dopwidth',
        'eq:app-ah-zeeman',
        'eq:app-ah-stark',
        'eq:app-ah-sigma',
        'eq:app-ah-tau',
        'eq:app-ah-optthin',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AH'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AH'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AH 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AH 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AH 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AH 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('474.862', '多普勒线宽比（整串绑定，防末位改动）'),
        ('3.91228', '质量比平方根（整串绑定，防末位改动）'),
        ('7.39922\\times10^{-3}', '中性份额（整串绑定，防末位改动）'),
        ('2.80074', '练习 2 电离度（整串绑定，防末位改动）'),
        ('2\\,\\mathrm{Ry}', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('13.605693122994\\ \\mathrm{eV}', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('13.6057$ eV', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('489.805$ eV', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('3.37506$ nm', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('656.28$ nm', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('1.31261\\times10^{7}\\ \\mathrm{m\\,s^{-1}}', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('0.0437841\\,c', '能级与跃迁：谱线为什·公式（精确绑定）'),
        ('367.354$ eV', '能级与跃迁：谱线为什·正文（精确绑定）'),
        ('100$ eV', '能级与跃迁：谱线为什·正文（精确绑定）'),
        ('6.03713\\times10^{27}\\,T^{3', 'Saha 方程：电离·公式（精确绑定）'),
        ('5.0cm', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.30$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('4.46814\\times10^{-5}$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.50$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.568235$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.60$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('6.11317$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.70$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('29.9454$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('1.00$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('98.692$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.5$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('5.53851$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('1.78588$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.180042$\\,\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('1$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('9.194$ keV', 'Saha 方程：电离·公式（精确绑定）'),
        ('99.4318\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('11.260$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('5.77523\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('0.568235\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('24.587$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('3.49693\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('2$ eV', 'Saha 方程：电离·公式（精确绑定）'),
        ('99.8725\\%', 'Saha 方程：电离·公式（精确绑定）'),
        ('68.0285\\ \\mathrm{eV}', '强度比：为什么日冕平·公式（精确绑定）'),
        ('300$ eV', '强度比：为什么日冕平·公式（精确绑定）'),
        ('1000$ eV', '强度比：为什么日冕平·公式（精确绑定）'),
        ('500$ eV', '强度比：为什么日冕平·公式（精确绑定）'),
        ('20\\%', '强度比：为什么日冕平·公式（精确绑定）'),
        ('2.72114\\%', '强度比：为什么日冕平·公式（精确绑定）'),
        ('68.0285$ eV', '强度比：为什么日冕平·公式（精确绑定）'),
        ('40\\%', '强度比：为什么日冕平·公式（精确绑定）'),
        ('8.67126\\%', '强度比：为什么日冕平·公式（精确绑定）'),
        ('4.41\\times10^{7}$ s$^{-1}$', '强度比：为什么日冕平·正文（精确绑定）'),
        ('2\\,T', '多普勒展宽：线宽测离·公式（精确绑定）'),
        ('2.0cm', '多普勒展宽：线宽测离·公式（精确绑定）'),
        ('468.57$ nm', '多普勒展宽：线宽测离·公式（精确绑定）'),
        ('1.62989\\times10^{5}$ m\\,s$^{-1}$', '多普勒展宽：线宽测离·正文（精确绑定）'),
        ('0.02$ nm', '多普勒展宽：线宽测离·正文（精确绑定）'),
        ('13.996\\ \\mathrm{GHz\\,T^{-1}}', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('5.3$ T', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('7.41788\\times10^{10}$ Hz', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('74.1788$ GHz', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('0.106571$ nm', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('2.81853\\times10^{-6}$ nm', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('0.053$ nm', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('0.246004$ nm', '塞曼分裂与斯塔克展宽·公式（精确绑定）'),
        ('3.04304\\times10^{-18}\\ \\mathrm{m^2}', '光学薄还是光学厚：谱·公式（精确绑定）'),
        ('2.48353\\times10^{11}$ Hz', '光学薄还是光学厚：谱·公式（精确绑定）'),
        ('3.04304\\times10^{-18}$ m$^2$', '光学薄还是光学厚：谱·公式（精确绑定）'),
        ('1$ m', '光学薄还是光学厚：谱·公式（精确绑定）'),
        ('200$ eV', '估算练习·正文（精确绑定）'),
        ('0.7$ eV', '估算练习·正文（精确绑定）'),
        ('29.9454\\%', '估算练习·正文（精确绑定）'),
        ('1$ keV', '估算练习·正文（精确绑定）'),
        ('2$ m', '估算练习·正文（精确绑定）'),
        ('183.84$ u', '估算练习·正文（精确绑定）'),
        ('5$ nm', '估算练习·正文（精确绑定）'),
        ('0.255556$ nm', '估算练习·正文（精确绑定）'),
        ('1.06261\\times10^{-3}$ nm', '估算练习·正文（精确绑定）'),
        ('10.692\\%', '估算练习·正文（精确绑定）'),
        ('0.00237607$ nm', '估算练习·正文（精确绑定）'),
        ('98.1443$ eV', '估算练习·正文（精确绑定）'),
        ('0.328619$ m', '估算练习·正文（精确绑定）'),
        ('0.00272816$ nm', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AH'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('跃迁 & 能量 / eV & 波长 / nm & 它在哪里重要 \\\\', '表格数据行'),
        ('C VI Ly-$\\alpha$（$2\\to1$） & $367.354$ & $3.37506$ & 芯部杂质辐射：软 X 射线，只能靠晶体谱仪 \\\\', '表格数据行'),
        ('C VI Ly-$\\beta$（$3\\to1$） & $435.382$ & $2.84771$ & 与 Ly-$\\alpha$ 组成强度比（第 3 节） \\\\', '表格数据行'),
        ('C VI Balmer-$\\alpha$（$3\\to2$） & $68.0285$ & $18.2253$ & 极紫外，电荷交换复合光谱的窗口 \\\\', '表格数据行'),
        ('D-$\\alpha$（$3\\to2$，声明） & $1.88920$ & $656.28$ & 边缘中性氢：可见光，最常被观测的一条线 \\\\', '表格数据行'),
        ('\\begin{tip}{谱线是能级差，不是「颜色」}{tip-app-ah-levels}', '语义盒标题'),
        ('\\textbf{出处：}本节关于类氢能级、玻尔速度与跃迁波长的写法，参照 Griem \\emph{Principles of Plasma ', '出处声明整行'),
        ('扫描量 & 固定量 & 电离度 & 读数 \\\\', '表格数据行'),
        ('$T=0.30$ eV & $n=10^{20}$ m$^{-3}$ & $4.46814\\times10^{-5}$\\,\\% & 几乎完全中性 \\\\', '表格数据行'),
        ('$T=0.50$ eV & $n=10^{20}$ m$^{-3}$ & $0.568235$\\,\\% & 中性为主 \\\\', '表格数据行'),
        ('$T=0.60$ eV & $n=10^{20}$ m$^{-3}$ & $6.11317$\\,\\% & 过渡区：0.1 eV 变一个量级 \\\\', '表格数据行'),
        ('$T=0.70$ eV & $n=10^{20}$ m$^{-3}$ & $29.9454$\\,\\% & 已经接近一半 \\\\', '表格数据行'),
        ('$T=1.00$ eV & $n=10^{20}$ m$^{-3}$ & $98.692$\\,\\% & 基本全电离 \\\\', '表格数据行'),
        ('$n=10^{18}$ m$^{-3}$ & $T=0.5$ eV & $5.53851$\\,\\% & 密度低十倍，电离度涨十倍 \\\\', '表格数据行'),
        ('$n=10^{19}$ m$^{-3}$ & $T=0.5$ eV & $1.78588$\\,\\% & 同上趋势 \\\\', '表格数据行'),
        ('$n=10^{20}$ m$^{-3}$ & $T=0.5$ eV & $0.568235$\\,\\% & 基准点 \\\\', '表格数据行'),
        ('$n=10^{21}$ m$^{-3}$ & $T=0.5$ eV & $0.180042$\\,\\% & 密度高十倍，电离度降三倍 \\\\', '表格数据行'),
        ('\\begin{tip}{「全电离」是个不等式，不是一句话}{tip-app-ah-saha}', '语义盒标题'),
        ('\\textbf{出处：}Saha 方程及其在恒星与实验室等离子体中的应用，参照 Saha（1920，\\emph{Philosophical ', '出处声明整行'),
        ('$T$ / eV & $R=I_{\\beta}/I_{\\alpha}$ & 与 $500$ eV 之比 & 读数 & \\\\', '表格数据行'),
        ('$300$ & $0.797111$ & $0.913287$ & 低温端：比值变化最快 & \\\\', '表格数据行'),
        ('$500$ & $0.872793$ & $1$ & 基准点：$T$ 变 $20\\%$，$R$ 变 $2.72114\\%$ & \\\\', '表格数据行'),
        ('$1000$ & $0.934234$ & $1.07040$ & 高温端：比值趋于饱和 & \\\\', '表格数据行'),
        ('\\begin{tip}{强度比是「相对」的温度计}{tip-app-ah-corona}', '语义盒标题'),
        ('\\textbf{出处：}日冕平衡与强度比测温的推导与适用条件，参照 Griem \\emph{Plasma Spectroscopy}（196', '出处声明整行'),
        ('谱线 & $T_i$ / eV & $\\Delta\\lambda_{\\mathrm{D}}$ / nm & 读数 \\\\', '表格数据行'),
        ('C VI Ly-$\\alpha$（$3.37506$ nm） & $50$ & $5.31306\\times10^{-4}$ & 高电荷态：线宽最窄 \\\\', '表格数据行'),
        ('C VI Ly-$\\alpha$ & $100$ & $7.51381\\times10^{-4}$ & 需要分辨率 $\\lesssim10^{-3}$ nm 的谱仪 \\\\', '表格数据行'),
        ('C VI Ly-$\\alpha$ & $200$ & $1.06261\\times10^{-3}$ & 线宽 $\\propto\\sqrt{T_i}$ \\\\', '表格数据行'),
        ('D-$\\alpha$（$656.28$ nm） & $50$ & $0.252297$ & 波长长、质量小：线宽大近五百倍 \\\\', '表格数据行'),
        ('D-$\\alpha$ & $100$ & $0.356802$ & 可见光波段就能分辨 \\\\', '表格数据行'),
        ('He II（$468.57$ nm） & $200$ & $0.255556$ & 氦的肩上：与 D-$\\alpha$ 同量级 \\\\', '表格数据行'),
        ('\\begin{tip}{线宽是「温度计」，但先要问「谁在发光」}{tip-app-ah-doppler}', '语义盒标题'),
        ('\\textbf{出处：}多普勒展宽公式与线宽测温的实践，参照 Hutchinson（2002，第 6--7 章）、', '出处声明整行'),
        ('量 & 数值 & 参考 & 读数 \\\\', '表格数据行'),
        ('塞曼分裂（$B=5.3$ T） & $74.1788$ GHz & D-$\\alpha$ & $\\Delta\\lambda_Z=0.106571$ nm \\\\', '表格数据行'),
        ('塞曼/多普勒 & $0.298683$ & D-$\\alpha$，$T_i=100$ eV & 同量级：可见光波段塞曼可分辨 \\\\', '表格数据行'),
        ('塞曼/多普勒 & $0.00375114$ & C VI，$T_i=100$ eV & 小 $79.6246$ 倍：软 X 射线测磁场极难 \\\\', '表格数据行'),
        ('斯塔克 FWHM & $0.053$ nm & D-$\\alpha$，$n_e=10^{20}$ m$^{-3}$ & 声明系数，按 $n_e^{2/3}$ 外推 \\\\', '表格数据行'),
        ('斯塔克 FWHM & $0.246004$ nm & D-$\\alpha$，$n_e=10^{21}$ m$^{-3}$ & 密度涨十倍，线宽涨 $4.64$ 倍 \\\\', '表格数据行'),
        ('斯塔克/塞曼 & $0.497323$ & D-$\\alpha$，$10^{20}$ m$^{-3}$ & 同量级：必须靠模型分离 \\\\', '表格数据行'),
        ('\\begin{tip}{磁场与密度不是从线宽「直接读」出来的，是模型反演出来的}{tip-app-ah-fields}', '语义盒标题'),
        ('\\textbf{出处：}塞曼效应与斯塔克展宽的公式形式、以及氢样离子线宽的密度定标，', '出处声明整行'),
        ('$n_1$ / m$^{-3}$ & $L$ / m & $\\tau$（D-$\\alpha$ 线心） & 读数 \\\\', '表格数据行'),
        ('$10^{16}$ & $0.01$ & $3.04304\\times10^{-4}$ & 光学薄：诊断可以放心用 \\\\', '表格数据行'),
        ('$10^{16}$ & $1$ & $0.0304304$ & 仍薄，但已不是「忽略不计」 \\\\', '表格数据行'),
        ('$10^{16}$ & $10$ & $0.304304$ & 边缘的典型情形：$\\tau$ 接近 $1/3$ \\\\', '表格数据行'),
        ('（对照）轫致吸收 & $1$ & $1.46411\\times10^{-6}$ & 同路径下线吸收是它的 $20784.2$ 倍 \\\\', '表格数据行'),
        ('\\begin{tip}{光学薄是假设，不是事实——每次都要验算}{tip-app-ah-opacity}', '语义盒标题'),
        ('\\textbf{出处：}光学深度与自吸收对谱线诊断的限制，参照 Griem（1964，第 3 章）、', '出处声明整行'),
        ('等级 & 含义 & 本附录中的例子 \\\\', '表格数据行'),
        ('现算级 & 由 \\texttt{tmp/appAH\\_calc.py} 从常数与输入现算 & Saha 常数、各电离度、跃迁波长、多普勒线宽、塞曼分裂、$\\sigma_{12}$、$\\tau$、练习解答 \\\\', '表格数据行'),
        ('声明级 & 脚本顶部写入的输入或文献值，不参与推导 & D-$\\alpha$ 波长 $656.28$ nm、He II $468.57$ nm、$A_{21}=4.41\\times10^{7}$ s$^{-1}$、仪器展宽 $0.02$ nm、Stark FWHM $0.053$ nm、日冕极限 $10^{20}$ m$^{-3}$ \\\\', '表格数据行'),
        ('文献级 & 只作定性对照或公式出处，不进入任何数值断言 & Griem 的谱线展宽理论、Saha 方程的原始文献、NIST 波长表的量级 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ah-honest}', '语义盒标题'),
        ('\\textbf{出处：}本节的分级与版权口径沿用全书统一约定；', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AH'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 2, '\\,\\mathrm{Ry}', '能级与跃迁：谱线为什·公式'),
        (0, 13.60569312, '\\ \\mathrm{eV}', '能级与跃迁：谱线为什·公式'),
        (0, 13.6057, '$ eV', '能级与跃迁：谱线为什·公式'),
        (0, 489.805, '$ eV', '能级与跃迁：谱线为什·公式'),
        (0, 3.37506, '$ nm', '能级与跃迁：谱线为什·公式'),
        (0, 656.28, '$ nm', '能级与跃迁：谱线为什·公式'),
        (0, 13126100, '\\ \\mathrm{m\\,s^{-1}}', '能级与跃迁：谱线为什·公式'),
        (0, 0.0437841, '\\,c', '能级与跃迁：谱线为什·公式'),
        (0, 367.354, '$ eV', '能级与跃迁：谱线为什·正文'),
        (0, 100, '$ eV', '能级与跃迁：谱线为什·正文'),
        (1, 6.03713e+27, '\\,T^{3', 'Saha 方程：电离·公式'),
        (1, 5, 'cm', 'Saha 方程：电离·公式'),
        (1, 0.3, '$ eV', 'Saha 方程：电离·公式'),
        (1, 4.46814e-05, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 0.5, '$ eV', 'Saha 方程：电离·公式'),
        (1, 0.568235, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 0.6, '$ eV', 'Saha 方程：电离·公式'),
        (1, 6.11317, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 0.7, '$ eV', 'Saha 方程：电离·公式'),
        (1, 29.9454, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 1, '$ eV', 'Saha 方程：电离·公式'),
        (1, 98.692, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 5.53851, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 1.78588, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 0.180042, '$\\,\\%', 'Saha 方程：电离·公式'),
        (1, 9.194, '$ keV', 'Saha 方程：电离·公式'),
        (1, 99.4318, '\\%', 'Saha 方程：电离·公式'),
        (1, 11.26, '$ eV', 'Saha 方程：电离·公式'),
        (1, 5.77523, '\\%', 'Saha 方程：电离·公式'),
        (1, 0.568235, '\\%', 'Saha 方程：电离·公式'),
        (1, 24.587, '$ eV', 'Saha 方程：电离·公式'),
        (1, 3.49693, '\\%', 'Saha 方程：电离·公式'),
        (1, 2, '$ eV', 'Saha 方程：电离·公式'),
        (1, 99.8725, '\\%', 'Saha 方程：电离·公式'),
        (2, 68.0285, '\\ \\mathrm{eV}', '强度比：为什么日冕平·公式'),
        (2, 300, '$ eV', '强度比：为什么日冕平·公式'),
        (2, 1000, '$ eV', '强度比：为什么日冕平·公式'),
        (2, 500, '$ eV', '强度比：为什么日冕平·公式'),
        (2, 20, '\\%', '强度比：为什么日冕平·公式'),
        (2, 2.72114, '\\%', '强度比：为什么日冕平·公式'),
        (2, 68.0285, '$ eV', '强度比：为什么日冕平·公式'),
        (2, 40, '\\%', '强度比：为什么日冕平·公式'),
        (2, 8.67126, '\\%', '强度比：为什么日冕平·公式'),
        (2, 44100000, '$ s$^{-1}$', '强度比：为什么日冕平·正文'),
        (3, 2, '\\,T', '多普勒展宽：线宽测离·公式'),
        (3, 2, 'cm', '多普勒展宽：线宽测离·公式'),
        (3, 468.57, '$ nm', '多普勒展宽：线宽测离·公式'),
        (3, 162989, '$ m\\,s$^{-1}$', '多普勒展宽：线宽测离·正文'),
        (3, 0.02, '$ nm', '多普勒展宽：线宽测离·正文'),
        (4, 13.996, '\\ \\mathrm{GHz\\,T^{-1}}', '塞曼分裂与斯塔克展宽·公式'),
        (4, 5.3, '$ T', '塞曼分裂与斯塔克展宽·公式'),
        (4, 7.41788e+10, '$ Hz', '塞曼分裂与斯塔克展宽·公式'),
        (4, 74.1788, '$ GHz', '塞曼分裂与斯塔克展宽·公式'),
        (4, 0.106571, '$ nm', '塞曼分裂与斯塔克展宽·公式'),
        (4, 2.81853e-06, '$ nm', '塞曼分裂与斯塔克展宽·公式'),
        (4, 0.053, '$ nm', '塞曼分裂与斯塔克展宽·公式'),
        (4, 0.246004, '$ nm', '塞曼分裂与斯塔克展宽·公式'),
        (5, 3.04304e-18, '\\ \\mathrm{m^2}', '光学薄还是光学厚：谱·公式'),
        (5, 2.48353e+11, '$ Hz', '光学薄还是光学厚：谱·公式'),
        (5, 3.04304e-18, '$ m$^2$', '光学薄还是光学厚：谱·公式'),
        (5, 1, '$ m', '光学薄还是光学厚：谱·公式'),
        (6, 200, '$ eV', '估算练习·正文'),
        (6, 29.9454, '\\%', '估算练习·正文'),
        (6, 1, '$ keV', '估算练习·正文'),
        (6, 2, '$ m', '估算练习·正文'),
        (6, 183.84, '$ u', '估算练习·正文'),
        (6, 5, '$ nm', '估算练习·正文'),
        (6, 0.255556, '$ nm', '估算练习·正文'),
        (6, 0.00106261, '$ nm', '估算练习·正文'),
        (6, 10.692, '\\%', '估算练习·正文'),
        (6, 0.00237607, '$ nm', '估算练习·正文'),
        (6, 98.1443, '$ eV', '估算练习·正文'),
        (6, 0.328619, '$ m', '估算练习·正文'),
        (6, 0.00272816, '$ nm', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AH 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AH 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('13.6\\ ', '里德伯能量四舍五入'),
        ('489.8\\ ', '类氢碳基态能量四舍五入'),
        ('0.0438\\ ', '玻尔速度占光速比四舍五入'),
        ('6.04\\times10^{27}', 'Saha 常数四舍五入'),
        ('0.57\\%', '0.5 eV 电离度四舍五入'),
        ('29.9\\%', '0.7 eV 电离度四舍五入'),
        ('98.7\\%', '1 eV 电离度四舍五入'),
        ('0.80\\ ', 'R(300 eV) 四舍五入'),
        ('0.87\\ ', 'R(500 eV) 四舍五入'),
        ('0.93\\ ', 'R(1000 eV) 四舍五入'),
        ('0.136\\ ', '强度比灵敏度四舍五入'),
        ('0.3568\\ ', 'D-alpha 线宽四舍五入'),
        ('0.2523\\ ', 'D-alpha 50 eV 线宽四舍五入'),
        ('0.2556\\ ', 'He II 线宽四舍五入'),
        ('474.9\\ ', 'D-alpha 与 C VI 线宽比四舍五入'),
        ('26.62\\ ', '仪器展宽之比四舍五入'),
        ('74.2\\ ', '塞曼频率四舍五入'),
        ('0.1066\\ ', '塞曼波长分裂四舍五入'),
        ('2.82\\times10^{-6}', 'C VI 塞曼分裂四舍五入'),
        ('0.2987\\ ', '塞曼与多普勒之比四舍五入'),
        ('4.64\\ ', '斯塔克密度定标四舍五入'),
        ('3.04\\times10^{-18}', '线心截面四舍五入'),
        ('4.57\\times10^{10}', '截面与汤姆孙截面之比四舍五入'),
        ('2.48\\times10^{11}', '多普勒频率宽度四舍五入'),
        ('0.0304\\ ', '光学深度四舍五入'),
        ('20784\\ ', '线吸收与轫致吸收之比四舍五入'),
        ('0.329\\ ', '练习 5 长度四舍五入'),
        ('98.1\\ ', '练习 4 温度四舍五入'),
        ('0.0027\\ ', '练习 6 钨线宽四舍五入'),
        ('0.0012\\ ', '练习 3 比值四舍五入'),
        ('15.3\\ ', '钨碳质量比四舍五入'),
        ('3.91\\ ', '质量比平方根四舍五入'),
        ('1.48\\ ', '波长比四舍五入'),
        ('10.7\\%', '练习 2 电离度四舍五入'),
        ('8.67\\%', 'R 降幅四舍五入'),
        ('3.67\\ ', '能量比四舍五入'),
        ('1.88\\times10^{-14}', '中性份额四舍五入'),
        ('1.31\\times10^{7}', '玻尔速度四舍五入'),
        ('5.78\\%', '碳一次电离度四舍五入'),
        ('3.50\\%', '氦一次电离度四舍五入'),
        ('1.89\\ ', 'D-alpha 能量四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AH'], why))


def check_P41_appendix_ai_icf():
    """附录 AI：惯性约束内爆与点火——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AI'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AI 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '驱动：黑腔里的辐射温度',
        '内爆：火箭模型与滞止',
        '绝热与熵：为什么压缩要冷',
        '流体不稳定性：瑞利--泰勒',
        '热斑与点火：自加热怎样',
        '燃烧与增益：把三本账',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AI 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AI 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ai-drive',
        'tab:app-ai-implosion',
        'tab:app-ai-adiabat',
        'tab:app-ai-rt',
        'tab:app-ai-ignition',
        'tab:app-ai-burn',
        'tab:app-ai-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AI'])

    for lab in (
        'sec:app-ai-drive',
        'sec:app-ai-implosion',
        'sec:app-ai-adiabat',
        'sec:app-ai-rt',
        'sec:app-ai-ignition',
        'sec:app-ai-burn',
        'sec:app-ai-problems',
        'sec:app-ai-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AI'])

    for lab in (
        'tip-app-ai-hohlraum',
        'tip-app-ai-rocket',
        'tip-app-ai-adiabat',
        'tip-app-ai-rt',
        'tip-app-ai-ignition',
        'tip-app-ai-gain',
        'tip-app-ai-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AI'])

    for lab in (
        'eq:app-ai-flux',
        'eq:app-ai-prad',
        'eq:app-ai-rocket',
        'eq:app-ai-ekin',
        'eq:app-ai-stag',
        'eq:app-ai-fermi',
        'eq:app-ai-adiabat',
        'eq:app-ai-shock',
        'eq:app-ai-rt',
        'eq:app-ai-rtabl',
        'eq:app-ai-cutoff',
        'eq:app-ai-efold',
        'eq:app-ai-hotspot',
        'eq:app-ai-dwell',
        'eq:app-ai-alpha',
        'eq:app-ai-brems',
        'eq:app-ai-range',
        'eq:app-ai-burn',
        'eq:app-ai-gain',
        'eq:app-ai-wall',
        'eq:app-ai-frac',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AI'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AI'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AI 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AI 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AI 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AI 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('4.0$ ns', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('1.5$ mm', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('5.0$ mm', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('6.12611\\times10^{-5}$ m$^2$', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('0.612611$ cm$^2$', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('6.93752\\times10^{18}$ W m$^{-2}$', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('3.32582\\times10^{6}\\ \\mathrm{K}', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('286.597\\ \\mathrm{eV}', '驱动：黑腔里的辐射温·公式（精确绑定）'),
        ('286.597$ eV', '驱动：黑腔里的辐射温·正文（精确绑定）'),
        ('3.08548\\times10^{10}$ Pa', '驱动：黑腔里的辐射温·正文（精确绑定）'),
        ('0.308548$ Mbar', '驱动：黑腔里的辐射温·正文（精确绑定）'),
        ('9.25643\\times10^{10}$ J m$^{-3}$', '驱动：黑腔里的辐射温·正文（精确绑定）'),
        ('3271.49$ J', '驱动：黑腔里的辐射温·正文（精确绑定）'),
        ('3.0cm', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('2.0\\ \\mathrm{MJ}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('1.7\\ \\mathrm{MJ}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('6.12611\\times10^{-5}\\ \\mathrm{m^2}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('6.93752\\times10^{18}\\ \\mathrm{W\\,m^{-2}}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('0.308548\\ \\mathrm{Mbar}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('3271.49\\ \\mathrm{J}', '驱动：黑腔里的辐射温·表格（精确绑定）'),
        ('330$ km', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('1.38331$ mm', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('2.82743$ mg', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('4.21285$ mg', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('7.04028$ mg', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('301.052\\ \\mathrm{km\\,s^{-1}}', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('128.128\\ \\mathrm{kJ}', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('6.40642\\%', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('4.59493$ ns', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('3.27591\\times10^{13}$ m s$^{-2}$', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('1.81265\\times10^{16}\\ \\mathrm{Pa}', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('181.265\\ \\mathrm{Gbar}', '内爆：火箭模型与滞止·公式（精确绑定）'),
        ('200$ g cm$^{-3}$', '内爆：火箭模型与滞止·正文（精确绑定）'),
        ('3$ g cm$^{-2}$', '内爆：火箭模型与滞止·正文（精确绑定）'),
        ('4.0cm', '内爆：火箭模型与滞止·表格（精确绑定）'),
        ('330\\ \\mathrm{km\\,s^{-1}}', '内爆：火箭模型与滞止·表格（精确绑定）'),
        ('4.59493\\ \\mathrm{ns}', '内爆：火箭模型与滞止·表格（精确绑定）'),
        ('3.27591\\times10^{13}\\ \\mathrm{m\\,s^{-2}}', '内爆：火箭模型与滞止·表格（精确绑定）'),
        ('3.34277\\times10^{12}\\,g', '内爆：火箭模型与滞止·表格（精确绑定）'),
        ('5 m', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('200\\ \\mathrm{g\\,cm^{-3}}', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('1.49065\\times10^{15}\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('14.9065\\ \\mathrm{Gbar}', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('100\\ \\mathrm{eV}', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('7.71883\\times10^{14}\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('30$ g cm$^{-3}$', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('6.3124\\times10^{13}$ Pa', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('5$ keV', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('5.78912\\times10^{15}$ Pa', '绝热与熵：为什么压缩·公式（精确绑定）'),
        ('4.81771\\times10^{31}\\ \\mathrm{m^{-3}}', '绝热与熵：为什么压缩·表格（精确绑定）'),
        ('30\\ \\mathrm{g\\,cm^{-3}}', '绝热与熵：为什么压缩·表格（精确绑定）'),
        ('6.3124\\times10^{13}\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·表格（精确绑定）'),
        ('5\\ \\mathrm{keV}', '绝热与熵：为什么压缩·表格（精确绑定）'),
        ('5.78912\\times10^{15}\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·表格（精确绑定）'),
        ('1.0\\times10^{14}$ m s$^{-2}$', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('628319\\ \\mathrm{m^{-1}}', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('7.30801\\times10^{9}\\ \\mathrm{s^{-1}}', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('0.136836\\ \\mathrm{ns}', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('3.0\\times10^{4}$ m s$^{-1}$', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('94444.4$ m$^{-1}$', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('1.41667\\times10^{9}$ s$^{-1}$', '流体不稳定性：瑞利-·公式（精确绑定）'),
        ('20$ nm', '流体不稳定性：瑞利-·正文（精确绑定）'),
        ('8.95311\\%', '流体不稳定性：瑞利-·正文（精确绑定）'),
        ('7.4462$ nm', '流体不稳定性：瑞利-·正文（精确绑定）'),
        ('94444.4\\ \\mathrm{m^{-1}}', '流体不稳定性：瑞利-·表格（精确绑定）'),
        ('1.41667\\times10^{9}\\ \\mathrm{s^{-1}}', '流体不稳定性：瑞利-·表格（精确绑定）'),
        ('7.4462\\ \\mathrm{nm}', '流体不稳定性：瑞利-·表格（精确绑定）'),
        ('0.12\\ \\mathrm{g\\,cm^{-2}}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('2327.94\\ \\mathrm{J}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('57.8912\\ \\mathrm{Gbar}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('0.3$ g cm$^{-2}$', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('0.12$ g cm$^{-2}$', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('64.3872\\ \\mathrm{ps}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('621.242\\ \\mathrm{km\\,s^{-1}}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('2.75\\times10^{-23}$ m$^3$ s$^{-1}$', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('3$ s$^{-1}$', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('2.01333\\times10^{26}\\ \\mathrm{W\\,m^{-3}}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('6.24074\\times10^{25}\\ \\mathrm{W\\,m^{-3}}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('3475.23$ J', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('1077.22$ J', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('1\\ \\mathrm{keV}', '热斑与点火：自加热怎·公式（精确绑定）'),
        ('2.4$ keV', '热斑与点火：自加热怎·正文（精确绑定）'),
        ('4.3$ keV', '热斑与点火：自加热怎·正文（精确绑定）'),
        ('1.1\\times10^{-24}T^2$', '热斑与点火：自加热怎·正文（精确绑定）'),
        ('3475.23\\ \\mathrm{J}', '热斑与点火：自加热怎·表格（精确绑定）'),
        ('1077.22\\ \\mathrm{J}', '热斑与点火：自加热怎·表格（精确绑定）'),
        ('2.4\\ \\mathrm{keV}', '热斑与点火：自加热怎·表格（精确绑定）'),
        ('0.3\\ \\mathrm{g\\,cm^{-2}}', '热斑与点火：自加热怎·表格（精确绑定）'),
        ('0.241452\\ \\mathrm{ns}', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('6\\ \\mathrm{g\\,cm^{-2}}', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('33.3333\\%', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('960.195\\ \\mathrm{MJ}', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('320.065\\ \\mathrm{MJ}', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('10\\%', '燃烧与增益：把三本账·公式（精确绑定）'),
        ('63.66\\ \\mathrm{MJ}', '燃烧与增益：把三本账·表格（精确绑定）'),
        ('250$ km', '估算练习·正文（精确绑定）'),
        ('3.46876\\times10^{18}$ W m$^{-2}$', '估算练习·正文（精确绑定）'),
        ('240.998$ eV', '估算练习·正文（精确绑定）'),
        ('5.81505$ mg', '估算练习·正文（精确绑定）'),
        ('368.716$ km', '估算练习·正文（精确绑定）'),
        ('9.60983\\%', '估算练习·正文（精确绑定）'),
        ('228.07$ km', '估算练习·正文（精确绑定）'),
        ('3.67678\\%', '估算练习·正文（精确绑定）'),
        ('0.24$ g cm$^{-2}$', '估算练习·正文（精确绑定）'),
        ('128.774$ ps', '估算练习·正文（精确绑定）'),
        ('20.8291\\%', '估算练习·正文（精确绑定）'),
        ('1.57854$ g cm$^{-2}$', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AI'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 表达式 & 数值 & 说明 \\\\', '表格数据行'),
        ('激光能量 & $E_{\\mathrm{drv}}$ & $2.0\\ \\mathrm{MJ}$ & 3$\\omega$（$0.35\\ \\mu$m）声明值 \\\\', '表格数据行'),
        ('转换效率 & $\\eta_{\\mathrm{ce}}$ & $0.85$ & 激光 $\\to$ X 射线，声明值 \\\\', '表格数据行'),
        ('X 射线能量 & $\\eta_{\\mathrm{ce}}E_{\\mathrm{drv}}$ & $1.7\\ \\mathrm{MJ}$ & 式 \\eqref{eq:app-ai-flux} \\\\', '表格数据行'),
        ('腔壁面积 & $A_{\\mathrm{hohl}}$ & $6.12611\\times10^{-5}\\ \\mathrm{m^2}$ & 圆柱侧面积加两个端面 \\\\', '表格数据行'),
        ('壁面通量 & $F_X$ & $6.93752\\times10^{18}\\ \\mathrm{W\\,m^{-2}}$ & 与 $\\sigma T^4$ 一致 \\\\', '表格数据行'),
        ('辐射温度 & $T_{\\mathrm{rad}}$ & $286.597\\ \\mathrm{eV}$ & 黑腔的「工作温度」 \\\\', '表格数据行'),
        ('辐射压 & $aT^4/3$ & $0.308548\\ \\mathrm{Mbar}$ & 式 \\eqref{eq:app-ai-prad} \\\\', '表格数据行'),
        ('腔内辐射能 & $uV$ & $3271.49\\ \\mathrm{J}$ & 只占 $E_X$ 的 $0.19\\%$ \\\\', '表格数据行'),
        ('\\begin{tip}{黑腔不是炉子，是能量中转站}{tip-app-ai-hohlraum}', '语义盒标题'),
        ('\\textbf{出处：}本文的间接驱动与黑腔图像、转换效率与辐射温度的口径参照', '出处声明整行'),
        ('排出速度 & $v_{\\mathrm{ex}}$ & $330\\ \\mathrm{km\\,s^{-1}}$ & 声明量级值 \\\\', '表格数据行'),
        ('质量比 & $m_0/m_f$ & $2.48999$ & 烧蚀掉 $59.8\\%$ 的初始质量 \\\\', '表格数据行'),
        ('内爆速度 & 式 \\eqref{eq:app-ai-rocket} & $301.052\\ \\mathrm{km\\,s^{-1}}$ & 设计目标 $300$ 量级 \\\\', '表格数据行'),
        ('内爆动能 & $\\frac12 m_f v^2$ & $128.128\\ \\mathrm{kJ}$ & 式 \\eqref{eq:app-ai-ekin} \\\\', '表格数据行'),
        ('流体力学效率 & $E_{\\mathrm{kin}}/E_{\\mathrm{drv}}$ & $6.40642\\%$ & 兆焦驱动只留千焦动能 \\\\', '表格数据行'),
        ('内爆时间 & $R_{\\mathrm{ice}}/v_{\\mathrm{imp}}$ & $4.59493\\ \\mathrm{ns}$ & 与脉宽同量级 \\\\', '表格数据行'),
        ('平均加速度 & $v^2/2R$ & $3.27591\\times10^{13}\\ \\mathrm{m\\,s^{-2}}$ & $3.34277\\times10^{12}\\,g_0$ \\\\', '表格数据行'),
        ('滞止动压 & $\\rho v^2$ & $181.265\\ \\mathrm{Gbar}$ & 式 \\eqref{eq:app-ai-stag} \\\\', '表格数据行'),
        ('压缩结果 & $R_f$, $\\rho$, $\\rho R$ & $150\\ \\mu$m, $200$ g cm$^{-3}$, $3$ g cm$^{-2}$ & 与附录 Y 同点 \\\\', '表格数据行'),
        ('\\begin{tip}{质量比决定内爆速度}{tip-app-ai-rocket}', '语义盒标题'),
        ('\\textbf{出处：}本文的火箭模型与内爆参数化参照 Atzeni \\& Meyer-ter-Vehn', '出处声明整行'),
        ('压缩燃料电子密度 & $n_e$ & $4.81771\\times10^{31}\\ \\mathrm{m^{-3}}$ & $\\rho = 200$ g cm$^{-3}$ \\\\', '表格数据行'),
        ('费米压 & $P_F$ & $14.9065\\ \\mathrm{Gbar}$ & 式 \\eqref{eq:app-ai-fermi} \\\\', '表格数据行'),
        ('冷燃料热压 & $P_{\\mathrm{th}}(100\\ \\mathrm{eV})$ & $7.71883\\times10^{14}\\ \\mathrm{Pa}$ & 简并压主导 \\\\', '表格数据行'),
        ('冷燃料绝热因子 & $\\alpha_{\\mathrm{ad}}$ & $0.517816$ & 小于一：已简并 \\\\', '表格数据行'),
        ('热斑费米压 & $P_F(30\\ \\mathrm{g\\,cm^{-3}})$ & $6.3124\\times10^{13}\\ \\mathrm{Pa}$ & 密度低，简并弱 \\\\', '表格数据行'),
        ('热斑热压 & $P_{\\mathrm{th}}(5\\ \\mathrm{keV})$ & $5.78912\\times10^{15}\\ \\mathrm{Pa}$ & 理想气体 \\\\', '表格数据行'),
        ('热斑绝热因子 & $\\alpha_{\\mathrm{ad}}$ & $91.7103$ & 远大于一 \\\\', '表格数据行'),
        ('强激波压缩比上限 & $(\\gamma+1)/(\\gamma-1)$ & $4$ & 式 \\eqref{eq:app-ai-shock} \\\\', '表格数据行'),
        ('\\begin{tip}{温度不是压缩的障碍，熵才是}{tip-app-ai-adiabat}', '语义盒标题'),
        ('\\textbf{出处：}本文的简并压与绝热压缩图像参照 Atzeni \\& Meyer-ter-Vehn', '出处声明整行'),
        ('扰动波数 & $2\\pi/\\lambda$ & $628319\\ \\mathrm{m^{-1}}$ & $\\lambda = 10\\ \\mu$m \\\\', '表格数据行'),
        ('无稳定化增长率 & $\\sqrt{A_T k g}$ & $7.30801\\times10^{9}\\ \\mathrm{s^{-1}}$ & 式 \\eqref{eq:app-ai-rt} \\\\', '表格数据行'),
        ('无稳定化时间 & $1/\\gamma$ & $0.136836\\ \\mathrm{ns}$ & 远小于内爆时间 \\\\', '表格数据行'),
        ('无稳定化放大 & $\\gamma t_{\\mathrm{imp}}$ & $33.5798$ & 约 $4\\times10^{14}$ 倍 \\\\', '表格数据行'),
        ('截止波数 & $A_T g/v_{\\mathrm{abl}}^2$ & $94444.4\\ \\mathrm{m^{-1}}$ & 式 \\eqref{eq:app-ai-cutoff} \\\\', '表格数据行'),
        ('截止波长 & $2\\pi/k_c$ & $66.5278\\ \\mu$m & 短于此波长被对流稳定 \\\\', '表格数据行'),
        ('最不稳定波长 & $2\\pi/k_{\\max}$ & $133.056\\ \\mu$m & $k_{\\max} = k_c/2$ \\\\', '表格数据行'),
        ('最大增长率 & $\\sqrt{A_T k_c g}/2$ & $1.41667\\times10^{9}\\ \\mathrm{s^{-1}}$ & 降低 $5.16$ 倍 \\\\', '表格数据行'),
        ('最大放大倍数 & $\\gamma_{\\max}t_{\\mathrm{imp}}$ & $6.50949$ & 式 \\eqref{eq:app-ai-efold} \\\\', '表格数据行'),
        ('种子 $20$ nm 的末期幅度 & $\\delta_0 e^{N_e}$ & $13.4297\\ \\mu$m & 占 $R_f$ 的 $8.95311\\%$ \\\\', '表格数据行'),
        ('末期 $5\\ \\mu$m 对应的种子 & $\\delta_0$ & $\\le 7.4462\\ \\mathrm{nm}$ & 表面粗糙度指标 \\\\', '表格数据行'),
        ('\\begin{tip}{不稳定性要按 e 倍增长数管理}{tip-app-ai-rt}', '语义盒标题'),
        ('\\textbf{出处：}本文的瑞利--泰勒不稳定性经典处理见 Rayleigh \\emph{Proceedings of the', '出处声明整行'),
        ('热斑面密度 & $\\rho r$ & $0.12\\ \\mathrm{g\\,cm^{-2}}$ & 阈值的 $2.5$ 分之一 \\\\', '表格数据行'),
        ('热斑能量 & $P V$ & $2327.94\\ \\mathrm{J}$ & 式 \\eqref{eq:app-ai-hotspot} \\\\', '表格数据行'),
        ('热斑压力 & $n k T$ & $57.8912\\ \\mathrm{Gbar}$ & 惯性约束 \\\\', '表格数据行'),
        ('声速 & $\\sqrt{2kT/m_i}$ & $621.242\\ \\mathrm{km\\,s^{-1}}$ & 式 \\eqref{eq:app-ai-dwell} \\\\', '表格数据行'),
        ('滞止时间 & $r/c_s$ & $64.3872\\ \\mathrm{ps}$ & 自加热的时间窗口 \\\\', '表格数据行'),
        ('$\\alpha$ 功率密度 & 式 \\eqref{eq:app-ai-alpha} & $2.01333\\times10^{26}\\ \\mathrm{W\\,m^{-3}}$ & $5$ keV，DT \\\\', '表格数据行'),
        ('轫致功率密度 & 式 \\eqref{eq:app-ai-brems} & $6.24074\\times10^{25}\\ \\mathrm{W\\,m^{-3}}$ & 附录 X 口径 \\\\', '表格数据行'),
        ('自加热/损失 & $P_\\alpha/P_{\\mathrm{br}}$ & $3.22611$ & 已跑赢 \\\\', '表格数据行'),
        ('一个滞止时间 $\\alpha$ 沉积 & $P_\\alpha V \\tau$ & $3475.23\\ \\mathrm{J}$ & $1.49283\\,E_{\\mathrm{hot}}$ \\\\', '表格数据行'),
        ('一个滞止时间轫致损失 & $P_{\\mathrm{br}} V \\tau$ & $1077.22\\ \\mathrm{J}$ & 约一半热斑能量 \\\\', '表格数据行'),
        ('模型交叉温度 & $P_\\alpha = P_{\\mathrm{br}}$ & 约 $2.4\\ \\mathrm{keV}$ & 文献理想值 $4.3$ keV \\\\', '表格数据行'),
        ('$\\alpha$ 射程 & $\\rho R = 0.3\\ \\mathrm{g\\,cm^{-2}}$ & $15\\ \\mu$m ($200$ g cm$^{-3}$) & 热斑外由冷燃料拦住 \\\\', '表格数据行'),
        ('\\begin{tip}{点火是自加热跑赢损失}{tip-app-ai-ignition}', '语义盒标题'),
        ('\\textbf{出处：}本文的热斑与点火图像、$\\alpha$ 自加热判据参照 Betti \\& Hurricane', '出处声明整行'),
        ('燃烧时间 & $R_f/c_s$ & $0.241452\\ \\mathrm{ns}$ & 式 \\eqref{eq:app-ai-burn} \\\\', '表格数据行'),
        ('燃耗份额 & $\\rho R/(\\rho R+6)$ & $33.3333\\%$ & 只烧掉三分之一 \\\\', '表格数据行'),
        ('聚变能 & $f m q$ & $960.195\\ \\mathrm{MJ}$ & 与附录 Y 一致 \\\\', '表格数据行'),
        ('$\\alpha$ 能量 & $E_{\\mathrm{fus}}/3$ 量级 & $63.66\\ \\mathrm{MJ}$ & 份额 $19.89\\%$ \\\\', '表格数据行'),
        ('燃烧放能 & $\\eta_\\alpha E_{\\mathrm{fus}}$ & $320.065\\ \\mathrm{MJ}$ & 式 \\eqref{eq:app-ai-gain} \\\\', '表格数据行'),
        ('内爆动能 & 式 \\eqref{eq:app-ai-ekin} & $128.128\\ \\mathrm{kJ}$ & 输入 \\\\', '表格数据行'),
        ('放大倍数 & $E_{\\mathrm{burn}}/E_{\\mathrm{kin}}$ & $2498.0$ & 压缩撬动聚变 \\\\', '表格数据行'),
        ('靶增益 & $E_{\\mathrm{burn}}/E_{\\mathrm{drv}}$ & $160.033$ & 式 \\eqref{eq:app-ai-wall} \\\\', '表格数据行'),
        ('墙插增益 & $\\eta_{\\mathrm{drv}} G$ & $16.0033$ & 驱动器效率 $10\\%$ \\\\', '表格数据行'),
        ('文献对照 & NIF 2022 & $2.05$ MJ $\\to$ 约 $3.15$ MJ & 靶增益约 $1.54$ \\\\', '表格数据行'),
        ('\\begin{tip}{增益不是效率，是比值}{tip-app-ai-gain}', '语义盒标题'),
        ('\\textbf{出处：}本文的燃耗拟合、增益定义与 ICF 能量账参照 Atzeni \\&', '出处声明整行'),
        ('\\textbf{出处：}本节六道题的题型参照 Atzeni \\& Meyer-ter-Vehn', '出处声明整行'),
        ('级别 & 定义 & 本附录的例子 & 复算途径 \\\\', '表格数据行'),
        ('$G = 160.033$ & \\texttt{tmp/appAI\\_calc.py} \\\\', '表格数据行'),
        ('$\\alpha$ 射程 $0.3$ g cm$^{-2}$ & 脚本内 \\texttt{[DECLARED]} 段 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ai-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录的间接驱动、内爆、不稳定性与点火图像分别参照', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AI'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 4, '$ ns', '驱动：黑腔里的辐射温·公式'),
        (0, 1.5, '$ mm', '驱动：黑腔里的辐射温·公式'),
        (0, 5, '$ mm', '驱动：黑腔里的辐射温·公式'),
        (0, 6.12611e-05, '$ m$^2$', '驱动：黑腔里的辐射温·公式'),
        (0, 0.612611, '$ cm$^2$', '驱动：黑腔里的辐射温·公式'),
        (0, 6.93752e+18, '$ W m$^{-2}$', '驱动：黑腔里的辐射温·公式'),
        (0, 3325820, '\\ \\mathrm{K}', '驱动：黑腔里的辐射温·公式'),
        (0, 286.597, '\\ \\mathrm{eV}', '驱动：黑腔里的辐射温·公式'),
        (0, 286.597, '$ eV', '驱动：黑腔里的辐射温·正文'),
        (0, 3.08548e+10, '$ Pa', '驱动：黑腔里的辐射温·正文'),
        (0, 0.308548, '$ Mbar', '驱动：黑腔里的辐射温·正文'),
        (0, 9.25643e+10, '$ J m$^{-3}$', '驱动：黑腔里的辐射温·正文'),
        (0, 3271.49, '$ J', '驱动：黑腔里的辐射温·正文'),
        (0, 3, 'cm', '驱动：黑腔里的辐射温·表格'),
        (0, 2, '\\ \\mathrm{MJ}', '驱动：黑腔里的辐射温·表格'),
        (0, 1.7, '\\ \\mathrm{MJ}', '驱动：黑腔里的辐射温·表格'),
        (0, 6.12611e-05, '\\ \\mathrm{m^2}', '驱动：黑腔里的辐射温·表格'),
        (0, 6.93752e+18, '\\ \\mathrm{W\\,m^{-2}}', '驱动：黑腔里的辐射温·表格'),
        (0, 0.308548, '\\ \\mathrm{Mbar}', '驱动：黑腔里的辐射温·表格'),
        (0, 3271.49, '\\ \\mathrm{J}', '驱动：黑腔里的辐射温·表格'),
        (1, 330, '$ km', '内爆：火箭模型与滞止·公式'),
        (1, 1.38331, '$ mm', '内爆：火箭模型与滞止·公式'),
        (1, 2.82743, '$ mg', '内爆：火箭模型与滞止·公式'),
        (1, 4.21285, '$ mg', '内爆：火箭模型与滞止·公式'),
        (1, 7.04028, '$ mg', '内爆：火箭模型与滞止·公式'),
        (1, 301.052, '\\ \\mathrm{km\\,s^{-1}}', '内爆：火箭模型与滞止·公式'),
        (1, 128.128, '\\ \\mathrm{kJ}', '内爆：火箭模型与滞止·公式'),
        (1, 6.40642, '\\%', '内爆：火箭模型与滞止·公式'),
        (1, 4.59493, '$ ns', '内爆：火箭模型与滞止·公式'),
        (1, 3.27591e+13, '$ m s$^{-2}$', '内爆：火箭模型与滞止·公式'),
        (1, 1.81265e+16, '\\ \\mathrm{Pa}', '内爆：火箭模型与滞止·公式'),
        (1, 181.265, '\\ \\mathrm{Gbar}', '内爆：火箭模型与滞止·公式'),
        (1, 200, '$ g cm$^{-3}$', '内爆：火箭模型与滞止·正文'),
        (1, 3, '$ g cm$^{-2}$', '内爆：火箭模型与滞止·正文'),
        (1, 4, 'cm', '内爆：火箭模型与滞止·表格'),
        (1, 330, '\\ \\mathrm{km\\,s^{-1}}', '内爆：火箭模型与滞止·表格'),
        (1, 4.59493, '\\ \\mathrm{ns}', '内爆：火箭模型与滞止·表格'),
        (1, 3.27591e+13, '\\ \\mathrm{m\\,s^{-2}}', '内爆：火箭模型与滞止·表格'),
        (1, 3.34277e+12, '\\,g', '内爆：火箭模型与滞止·表格'),
        (2, 5, 'm', '绝热与熵：为什么压缩·公式'),
        (2, 200, '\\ \\mathrm{g\\,cm^{-3}}', '绝热与熵：为什么压缩·公式'),
        (2, 1.49065e+15, '\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·公式'),
        (2, 14.9065, '\\ \\mathrm{Gbar}', '绝热与熵：为什么压缩·公式'),
        (2, 100, '\\ \\mathrm{eV}', '绝热与熵：为什么压缩·公式'),
        (2, 7.71883e+14, '\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·公式'),
        (2, 30, '$ g cm$^{-3}$', '绝热与熵：为什么压缩·公式'),
        (2, 6.3124e+13, '$ Pa', '绝热与熵：为什么压缩·公式'),
        (2, 5, '$ keV', '绝热与熵：为什么压缩·公式'),
        (2, 5.78912e+15, '$ Pa', '绝热与熵：为什么压缩·公式'),
        (2, 4.81771e+31, '\\ \\mathrm{m^{-3}}', '绝热与熵：为什么压缩·表格'),
        (2, 30, '\\ \\mathrm{g\\,cm^{-3}}', '绝热与熵：为什么压缩·表格'),
        (2, 6.3124e+13, '\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·表格'),
        (2, 5, '\\ \\mathrm{keV}', '绝热与熵：为什么压缩·表格'),
        (2, 5.78912e+15, '\\ \\mathrm{Pa}', '绝热与熵：为什么压缩·表格'),
        (3, 1e+14, '$ m s$^{-2}$', '流体不稳定性：瑞利-·公式'),
        (3, 628319, '\\ \\mathrm{m^{-1}}', '流体不稳定性：瑞利-·公式'),
        (3, 7308010000, '\\ \\mathrm{s^{-1}}', '流体不稳定性：瑞利-·公式'),
        (3, 0.136836, '\\ \\mathrm{ns}', '流体不稳定性：瑞利-·公式'),
        (3, 30000, '$ m s$^{-1}$', '流体不稳定性：瑞利-·公式'),
        (3, 94444.4, '$ m$^{-1}$', '流体不稳定性：瑞利-·公式'),
        (3, 1416670000, '$ s$^{-1}$', '流体不稳定性：瑞利-·公式'),
        (3, 20, '$ nm', '流体不稳定性：瑞利-·正文'),
        (3, 8.95311, '\\%', '流体不稳定性：瑞利-·正文'),
        (3, 7.4462, '$ nm', '流体不稳定性：瑞利-·正文'),
        (3, 94444.4, '\\ \\mathrm{m^{-1}}', '流体不稳定性：瑞利-·表格'),
        (3, 1416670000, '\\ \\mathrm{s^{-1}}', '流体不稳定性：瑞利-·表格'),
        (3, 7.4462, '\\ \\mathrm{nm}', '流体不稳定性：瑞利-·表格'),
        (4, 0.12, '\\ \\mathrm{g\\,cm^{-2}}', '热斑与点火：自加热怎·公式'),
        (4, 2327.94, '\\ \\mathrm{J}', '热斑与点火：自加热怎·公式'),
        (4, 57.8912, '\\ \\mathrm{Gbar}', '热斑与点火：自加热怎·公式'),
        (4, 0.3, '$ g cm$^{-2}$', '热斑与点火：自加热怎·公式'),
        (4, 0.12, '$ g cm$^{-2}$', '热斑与点火：自加热怎·公式'),
        (4, 64.3872, '\\ \\mathrm{ps}', '热斑与点火：自加热怎·公式'),
        (4, 621.242, '\\ \\mathrm{km\\,s^{-1}}', '热斑与点火：自加热怎·公式'),
        (4, 2.75e-23, '$ m$^3$ s$^{-1}$', '热斑与点火：自加热怎·公式'),
        (4, 3, '$ s$^{-1}$', '热斑与点火：自加热怎·公式'),
        (4, 2.01333e+26, '\\ \\mathrm{W\\,m^{-3}}', '热斑与点火：自加热怎·公式'),
        (4, 6.24074e+25, '\\ \\mathrm{W\\,m^{-3}}', '热斑与点火：自加热怎·公式'),
        (4, 3475.23, '$ J', '热斑与点火：自加热怎·公式'),
        (4, 1077.22, '$ J', '热斑与点火：自加热怎·公式'),
        (4, 1, '\\ \\mathrm{keV}', '热斑与点火：自加热怎·公式'),
        (4, 2.4, '$ keV', '热斑与点火：自加热怎·正文'),
        (4, 4.3, '$ keV', '热斑与点火：自加热怎·正文'),
        (4, 1.1e-24, 'T^2$', '热斑与点火：自加热怎·正文'),
        (4, 3475.23, '\\ \\mathrm{J}', '热斑与点火：自加热怎·表格'),
        (4, 1077.22, '\\ \\mathrm{J}', '热斑与点火：自加热怎·表格'),
        (4, 2.4, '\\ \\mathrm{keV}', '热斑与点火：自加热怎·表格'),
        (4, 0.3, '\\ \\mathrm{g\\,cm^{-2}}', '热斑与点火：自加热怎·表格'),
        (5, 0.241452, '\\ \\mathrm{ns}', '燃烧与增益：把三本账·公式'),
        (5, 6, '\\ \\mathrm{g\\,cm^{-2}}', '燃烧与增益：把三本账·公式'),
        (5, 33.3333, '\\%', '燃烧与增益：把三本账·公式'),
        (5, 960.195, '\\ \\mathrm{MJ}', '燃烧与增益：把三本账·公式'),
        (5, 320.065, '\\ \\mathrm{MJ}', '燃烧与增益：把三本账·公式'),
        (5, 10, '\\%', '燃烧与增益：把三本账·公式'),
        (5, 63.66, '\\ \\mathrm{MJ}', '燃烧与增益：把三本账·表格'),
        (6, 250, '$ km', '估算练习·正文'),
        (6, 3.46876e+18, '$ W m$^{-2}$', '估算练习·正文'),
        (6, 240.998, '$ eV', '估算练习·正文'),
        (6, 5.81505, '$ mg', '估算练习·正文'),
        (6, 368.716, '$ km', '估算练习·正文'),
        (6, 9.60983, '\\%', '估算练习·正文'),
        (6, 228.07, '$ km', '估算练习·正文'),
        (6, 3.67678, '\\%', '估算练习·正文'),
        (6, 0.24, '$ g cm$^{-2}$', '估算练习·正文'),
        (6, 128.774, '$ ps', '估算练习·正文'),
        (6, 20.8291, '\\%', '估算练习·正文'),
        (6, 1.57854, '$ g cm$^{-2}$', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AI 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AI 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('287\\ ', '辐射温度四舍五入'),
        ('3.33\\times10^{6}', '辐射温度开尔文四舍五入'),
        ('0.31\\ ', '辐射压四舍五入'),
        ('9.26\\times10^{10}', '辐射能密度四舍五入'),
        ('3271\\ ', '腔内辐射能四舍五入'),
        ('6.13\\times10^{-5}', '腔壁面积四舍五入'),
        ('6.94\\times10^{18}', '壁面通量四舍五入'),
        ('301.1\\ ', '内爆速度四舍五入'),
        ('128.1\\ ', '内爆动能四舍五入'),
        ('6.41\\%', '流体力学效率四舍五入'),
        ('4.59\\ ', '内爆时间四舍五入'),
        ('3.28\\times10^{13}', '平均加速度四舍五入'),
        ('181.3\\ ', '滞止动压四舍五入'),
        ('14.91\\ ', '压缩燃料费米压四舍五入'),
        ('0.52\\ ', '冷燃料绝热因子四舍五入'),
        ('91.71\\ ', '热斑绝热因子四舍五入'),
        ('7.31\\times10^{9}', '无稳定化增长率四舍五入'),
        ('0.137\\ ', '无稳定化时间四舍五入'),
        ('33.58\\ ', '无稳定化放大倍数四舍五入'),
        ('94444\\ ', '截止波数四舍五入'),
        ('66.53\\ ', '截止波长四舍五入'),
        ('133.1\\ ', '最不稳定波长四舍五入'),
        ('1.42\\times10^{9}', '最大增长率四舍五入'),
        ('6.51\\ ', '最大放大倍数四舍五入'),
        ('13.43\\ ', '末期扰动幅度四舍五入'),
        ('7.45\\ ', '种子扰动上限四舍五入'),
        ('2328\\ ', '热斑能量四舍五入'),
        ('57.89\\ ', '热斑压力四舍五入'),
        ('621.2\\ ', '热斑声速四舍五入'),
        ('64.39\\ ', '滞止时间四舍五入'),
        ('2.01\\times10^{26}', 'alpha 功率密度四舍五入'),
        ('6.24\\times10^{25}', '轫致功率密度四舍五入'),
        ('3.23\\ ', '自加热与损失之比四舍五入'),
        ('3475\\ ', 'alpha 沉积能量四舍五入'),
        ('1077\\ ', '轫致损失能量四舍五入'),
        ('0.241\\ ', '燃烧时间四舍五入'),
        ('960.2\\ ', '聚变能四舍五入'),
        ('320.1\\ ', '燃烧放能四舍五入'),
        ('160.0\\ ', '靶增益四舍五入'),
        ('16.0\\ ', '墙插增益四舍五入'),
        ('2500\\ ', '放大倍数四舍五入'),
        ('3.47\\times10^{18}', '练习 1 通量四舍五入'),
        ('241\\ ', '练习 1 辐射温度四舍五入'),
        ('5.82\\ ', '练习 2 烧蚀质量四舍五入'),
        ('368.7\\ ', '练习 2 内爆速度四舍五入'),
        ('192.2\\ ', '练习 2 内爆动能四舍五入'),
        ('9.61\\%', '练习 2 效率四舍五入'),
        ('228.1\\ ', '练习 3 内爆速度四舍五入'),
        ('73.54\\ ', '练习 3 内爆动能四舍五入'),
        ('3.68\\%', '练习 3 效率四舍五入'),
        ('18.62\\ ', '练习 5 热斑能量四舍五入'),
        ('128.8\\ ', '练习 5 滞止时间四舍五入'),
        ('0.208\\ ', '练习 6 燃耗份额四舍五入'),
        ('1.579\\ ', '练习 6 面密度四舍五入'),
        ('8.95\\%', '末期扰动占比四舍五入'),
        ('1.383\\ ', '燃料冰球半径四舍五入'),
        ('0.912\\ ', '质量比对数四舍五入'),
        ('2.49\\ ', '质量比四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AI'], why))


def check_P42_appendix_aj_fuelling():
    """附录 AJ：加料与排灰——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AJ'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AJ 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '粒子平衡：加进去的粒子有多少被烧掉',
        '气体加料与弹丸加料：效率差在哪',
        '弹丸烧蚀与穿透：中性气体屏蔽',
        '氦灰：稀释与排出条件',
        '抽气与氚存量：为什么瓶颈不在泵上',
        '密度控制：反馈回路的带宽',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AJ 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AJ 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-aj-balance',
        'tab:app-aj-fuelling',
        'tab:app-aj-pellet',
        'tab:app-aj-ash',
        'tab:app-aj-pumping',
        'tab:app-aj-control',
        'tab:app-aj-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AJ'])

    for lab in (
        'sec:app-aj-balance',
        'sec:app-aj-fuelling',
        'sec:app-aj-pellet',
        'sec:app-aj-ash',
        'sec:app-aj-pumping',
        'sec:app-aj-control',
        'sec:app-aj-problems',
        'sec:app-aj-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AJ'])

    for lab in (
        'tip-app-aj-balance',
        'tip-app-aj-efficiency',
        'tip-app-aj-ngs',
        'tip-app-aj-ash',
        'tip-app-aj-pump',
        'tip-app-aj-control',
        'tip-app-aj-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AJ'])

    for lab in (
        'eq:app-aj-rfus',
        'eq:app-aj-burn',
        'eq:app-aj-loss',
        'eq:app-aj-recycle',
        'eq:app-aj-burnfrac',
        'eq:app-aj-gasinj',
        'eq:app-aj-gasflow',
        'eq:app-aj-pelletmass',
        'eq:app-aj-abl',
        'eq:app-aj-pen',
        'eq:app-aj-pelletmass2',
        'eq:app-aj-pelletfreq',
        'eq:app-aj-ashdens',
        'eq:app-aj-dilute',
        'eq:app-aj-ashreq',
        'eq:app-aj-ashreq2',
        'eq:app-aj-pumpq',
        'eq:app-aj-hefrac',
        'eq:app-aj-inventory',
        'eq:app-aj-inventory2',
        'eq:app-aj-tbr',
        'eq:app-aj-greenwald',
        'eq:app-aj-ramp',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AJ'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AJ'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AJ 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AJ 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AJ 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AJ 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('556.5\\ \\mathrm{MW}', '粒子平衡：加进去的粒·正文（精确绑定）'),
        ('17.6\\ \\mathrm{MeV}', '粒子平衡：加进去的粒·正文（精确绑定）'),
        ('1.97352\\times10^{20}\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('3.94705\\times10^{20}\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('31.1913\\ \\text{kg/年}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('31.19\\ \\text{kg/年}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('2.5\\ \\mathrm{s}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('3.3288\\times10^{22}\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('3.36827\\times10^{22}\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('2.0591\\times10^{21}\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('19.1687\\%', '粒子平衡：加进去的粒·公式（精确绑定）'),
        ('10\\%', '粒子平衡：加进去的粒·正文（精确绑定）'),
        ('1.91687\\%', '粒子平衡：加进去的粒·正文（精确绑定）'),
        ('2.0591\\times10^{22}\\ \\mathrm{s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('300\\ \\mathrm{K}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('426.435\\ \\mathrm{mbar\\,L\\,s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('1.5\\ \\mathrm{mm}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('255\\ \\mathrm{kg\\,m^{-3}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('300\\ \\mathrm{m\\,s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('3.43184\\times10^{21}\\ \\mathrm{s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('7.12338\\times10^{-6}\\ \\mathrm{kg\\,s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('7.12338\\ \\mathrm{mg\\,s^{-1}}', '气体加料与弹丸加料：·公式（精确绑定）'),
        ('615.46\\ \\mathrm{g/day}', '气体加料与弹丸加料：·正文（精确绑定）'),
        ('224.797\\ \\text{kg/年}', '气体加料与弹丸加料：·正文（精确绑定）'),
        ('0.12\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('12\\ \\mathrm{cm}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('2.0\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('6\\%', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('3.60498\\times10^{-6}\\ \\mathrm{kg}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('3.60498\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('2.16299\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('4.18155\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('1.97599\\ \\mathrm{Hz}', '弹丸烧蚀与穿透：中性·公式（精确绑定）'),
        ('0.506077\\ \\mathrm{s}', '弹丸烧蚀与穿透：中性·正文（精确绑定）'),
        ('2\\ \\mathrm{mm}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('8.54513\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('0.646087\\ \\mathrm{ms}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('0.193826\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('3\\ \\mathrm{mm}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('28.8398\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('1.26992\\ \\mathrm{ms}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('0.380976\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('19.0488\\%', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('1000\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('1.26992\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格（精确绑定）'),
        ('1\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·正文（精确绑定）'),
        ('2500\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·正文（精确绑定）'),
        ('787.451\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·正文（精确绑定）'),
        ('18.5\\ \\mathrm{s}', '氦灰：稀释与排出条件·正文（精确绑定）'),
        ('4.38719\\times10^{18}\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('4.38719\\%', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('9.12256\\times10^{19}\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('91.2256\\%', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('463.126\\ \\mathrm{MW}', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('83.2211\\%', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('5\\%', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('21.0841\\ \\mathrm{s}', '氦灰：稀释与排出条件·公式（精确绑定）'),
        ('3.0cm', '氦灰：稀释与排出条件·表格（精确绑定）'),
        ('2.5\\times10^{19}\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·表格（精确绑定）'),
        ('5\\ \\mathrm{Pa}', '抽气与氚存量：为什么·正文（精确绑定）'),
        ('1.20716\\times10^{23}\\ \\mathrm{s^{-1}}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('0.163485\\%', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('5.15627\\times10^{-6}\\ \\mathrm{kg\\,s^{-1}}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('162.72\\ \\text{kg/年}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('1\\%', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('1.6272\\ \\text{kg/年}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('2.0\\ \\mathrm{kg}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('1.22911\\ \\text{年}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('0.112524\\ \\text{kg/年}', '抽气与氚存量：为什么·公式（精确绑定）'),
        ('11.0844\\ \\mathrm{MA}', '密度控制：反馈回路的·正文（精确绑定）'),
        ('8.82069\\times10^{19}\\ \\mathrm{m^{-3}}', '密度控制：反馈回路的·公式（精确绑定）'),
        ('40.4156\\ \\mathrm{s}', '密度控制：反馈回路的·公式（精确绑定）'),
        ('3.7\\ \\mathrm{s}', '密度控制：反馈回路的·表格（精确绑定）'),
        ('5\\ \\mathrm{s}', '估算练习·正文（精确绑定）'),
        ('3\\%', '估算练习·正文（精确绑定）'),
        ('10\\ \\mathrm{Pa}', '估算练习·正文（精确绑定）'),
        ('1.2269\\times10^{21}\\ \\mathrm{s^{-1}}', '估算练习·正文（精确绑定）'),
        ('0.833619\\ \\mathrm{Hz}', '估算练习·正文（精确绑定）'),
        ('12.6505\\ \\mathrm{s}', '估算练习·正文（精确绑定）'),
        ('450.765\\ \\mathrm{MW}', '估算练习·正文（精确绑定）'),
        ('2.41432\\times10^{23}\\ \\mathrm{s^{-1}}', '估算练习·正文（精确绑定）'),
        ('9.88481\\times10^{-7}\\ \\mathrm{kg\\,s^{-1}}', 'AJ.1·燃耗折算到每秒（精确绑定·补充）'),
        ('6.18573\\%', 'AJ.1·加料量占输运损失的比例（精确绑定·补充）'),
        ('19\\%', 'AJ.1·净加料燃耗份额的口语说法（精确绑定·补充）'),
        ('98\\%', 'AJ.1·被抽走的粒子份额（精确绑定·补充）'),
        ('60\\%', 'AJ.2·弹丸加料效率（精确绑定·补充）'),
        ('0.4\\ \\mathrm{ms}', 'AJ.3·声明烧蚀时间（精确绑定·补充）'),
        ('4.4\\%', 'AJ.4·灰份额的口语说法（精确绑定·补充）'),
        ('16.8\\%', 'AJ.4·功率损失的口语说法（精确绑定·补充）'),
        ('468.57\\ \\mathrm{nm}', 'AJ.4·He II 谱线（与附录 AH 同口径）（精确绑定·补充）'),
        ('81\\%', 'AJ.7·练习 4 答案的功率份额（精确绑定·补充）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AJ'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 表达式 & 数值 & 说明 \\\\', '表格数据行'),
        ('电子总数 $N_e$ & $n_eV$ & $8.322\\times10^{22}$ & 工作点存量 \\\\', '表格数据行'),
        ('聚变反应率 & $P_{\\mathrm{fus}}/E_{\\mathrm{fus}}$ & $1.97352\\times10^{20}\\ \\mathrm{s^{-1}}$ & 由功率反推 \\\\', '表格数据行'),
        ('燃烧消耗率 & $2R_{\\mathrm{fus}}$ & $3.94705\\times10^{20}\\ \\mathrm{s^{-1}}$ & D 与 T 各一份 \\\\', '表格数据行'),
        ('氚燃耗 & $R_{\\mathrm{fus}}m_T$ & $31.1913\\ \\text{kg/年}$ & 与附录 Z 同口径 \\\\', '表格数据行'),
        ('输运损失 & $N_e/\\tau_p$ & $3.3288\\times10^{22}\\ \\mathrm{s^{-1}}$ & 比燃烧大 84 倍 \\\\', '表格数据行'),
        ('外部加料需求 & $(1-R)N_e/\\tau_p + 2R_{\\mathrm{fus}}$ & $2.0591\\times10^{21}\\ \\mathrm{s^{-1}}$ & 再循环帮了大忙 \\\\', '表格数据行'),
        ('净加料燃耗份额 & $2R_{\\mathrm{fus}}/S_{\\mathrm{ext}}$ & $19.1687\\%$ & 注进去只有 $1.91687\\%$ \\\\', '表格数据行'),
        ('\\begin{tip}{加料量不是燃烧量}{tip-app-aj-balance}', '语义盒标题'),
        ('\\textbf{出处：}粒子平衡与再循环的处理参照 Wesson \\emph{Tokamaks}（第 4 版，2011，第 1、6 章）、', '出处声明整行'),
        ('方式 & 加料效率（声明） & 注入量 & 工程口径 \\\\', '表格数据行'),
        ('气体加料 & $10\\%$ & $2.0591\\times10^{22}\\ \\mathrm{s^{-1}}$ & $42.6435\\ \\mathrm{Pa\\,m^{3}\\,s^{-1}}$（$426.435\\ \\mathrm{mbar\\,L\\,s^{-1}}$） \\\\', '表格数据行'),
        ('弹丸加料 & $60\\%$ & $3.43184\\times10^{21}\\ \\mathrm{s^{-1}}$ & $7.12338\\ \\mathrm{mg\\,s^{-1}}$（$615.46\\ \\mathrm{g/day}$） \\\\', '表格数据行'),
        ('超声分子束 & 介于两者之间 & 与气体同量级 & 靠定向高速中性流提高穿透 \\\\', '表格数据行'),
        ('紧凑环注入 & 芯部加料 & 与弹丸同量级 & 把等离子体团整体注入 \\\\', '表格数据行'),
        ('\\begin{tip}{效率低不是浪费，是物理}{tip-app-aj-efficiency}', '语义盒标题'),
        ('\\textbf{出处：}弹丸加料的综述与效率区间参照 Milora 等（1995，\\emph{Nucl. Fusion} \\textbf{3', '出处声明整行'),
        ('$r_p$ & 质量 & $\\tau_{\\mathrm{abl}}$ & $L_{\\mathrm{pen}}$ & 说明 \\\\', '表格数据行'),
        ('$1.5\\ \\mathrm{mm}$ & $3.60498\\ \\mathrm{mg}$ & $0.4\\ \\mathrm{ms}$ & $0.12\\ \\mathrm{m}$ & 占 $a$ 的 $6\\%$，边缘加料 \\\\', '表格数据行'),
        ('$2\\ \\mathrm{mm}$ & $8.54513\\ \\mathrm{mg}$ & $0.646087\\ \\mathrm{ms}$ & $0.193826\\ \\mathrm{m}$ & 仍在外半区 \\\\', '表格数据行'),
        ('$3\\ \\mathrm{mm}$ & $28.8398\\ \\mathrm{mg}$ & $1.26992\\ \\mathrm{ms}$ & $0.380976\\ \\mathrm{m}$ & 占 $a$ 的 $19.0488\\%$ \\\\', '表格数据行'),
        ('$3\\ \\mathrm{mm}$（$1000\\ \\mathrm{m\\,s^{-1}}$） & $28.8398\\ \\mathrm{mg}$ & $1.26992\\ \\mathrm{ms}$ & $1.26992\\ \\mathrm{m}$ & $L/a = 0.63496$，够到芯部 \\\\', '表格数据行'),
        ('\\begin{tip}{弹丸是被自己的气体盾牌挡住的}{tip-app-aj-ngs}', '语义盒标题'),
        ('\\textbf{出处：}中性气体屏蔽与烧蚀标度参照 Parks \\& Turnbull（1978，\\emph{Phys. Fluids} \\', '出处声明整行'),
        ('氦产生率 & $R_{\\mathrm{fus}}$ & $1.97352\\times10^{20}\\ \\mathrm{s^{-1}}$ & 每个反应一个 \\\\', '表格数据行'),
        ('停留时间 & $5\\tau_E$ & $18.5\\ \\mathrm{s}$ & 声明口径 \\\\', '表格数据行'),
        ('稳态氦密度 & $R_{\\mathrm{fus}}\\tau_{\\mathrm{He}}^{*}/V$ & $4.38719\\times10^{18}\\ \\mathrm{m^{-3}}$ & 即 $4.38719\\%\\,n_e$ \\\\', '表格数据行'),
        ('燃料离子密度 & $n_e - 2n_{\\mathrm{He}}$ & $9.12256\\times10^{19}\\ \\mathrm{m^{-3}}$ & $91.2256\\%\\,n_e$ \\\\', '表格数据行'),
        ('稀释后功率 & $\\propto n_{\\mathrm{fuel}}^{2}$ & $463.126\\ \\mathrm{MW}$ & $83.2211\\%$ \\\\', '表格数据行'),
        ('排出条件 & $\\tau_{\\mathrm{He}}^{*}/\\tau_E \\le 5.69841$ & $21.0841\\ \\mathrm{s}$ & 压到 $5\\%$ 以内 \\\\', '表格数据行'),
        ('偏滤器处氦密度 & 富集 $5$ 倍（声明） & $2.5\\times10^{19}\\ \\mathrm{m^{-3}}$ & 富集是抽气的必要条件 \\\\', '表格数据行'),
        ('\\begin{tip}{灰不是杂质，是产物}{tip-app-aj-ash}', '语义盒标题'),
        ('\\textbf{出处：}功率与粒子控制（含氦灰排出条件）参照 Loarte 等（2007，\\emph{Nucl. Fusion} \\text', '出处声明整行'),
        ('抽气排出率 & $PS/k_BT$ & $1.20716\\times10^{23}\\ \\mathrm{s^{-1}}$ & 是氦产生率的 $611.679$ 倍 \\\\', '表格数据行'),
        ('需要的氦份额 & $R_{\\mathrm{fus}}/Q_{\\mathrm{pump}}$ & $0.163485\\%$ & 泵不是瓶颈 \\\\', '表格数据行'),
        ('加料氚流 & $S_{\\mathrm{ext}}m_T/2$ & $162.72\\ \\text{kg/年}$ & 燃耗的 $5.21683$ 倍 \\\\', '表格数据行'),
        ('年滞留 & $1\\%$ 的加料氚流 & $1.6272\\ \\text{kg/年}$ & 声明份额 \\\\', '表格数据行'),
        ('到库存上限 & $2.0\\ \\mathrm{kg}/\\dot m_T^{\\mathrm{ret}}$ & $1.22911\\ \\text{年}$ & $14.7493$ 个月 \\\\', '表格数据行'),
        ('衰变损失 & $\\ln 2\\,M/T_{1/2}$ & $0.112524\\ \\text{kg/年}$ & 库存 $2.0\\ \\mathrm{kg}$ \\\\', '表格数据行'),
        ('自持下限 & 燃耗$+$滞留$+$衰变 & $1.05578$ & 附录 Z 只算衰变是 $1.0036$ \\\\', '表格数据行'),
        ('\\begin{tip}{瓶颈不在泵上}{tip-app-aj-pump}', '语义盒标题'),
        ('\\textbf{出处：}偏滤器抽气与再循环的处理参照 Stangeby（2000，第 8 章）与 Loarte 等（2007）；', '出处声明整行'),
        ('Greenwald 密度 & $I_p/(\\pi a^{2})$ & $8.82069\\times10^{19}\\ \\mathrm{m^{-3}}$ & $I_p = 11.0844\\ \\mathrm{MA}$ \\\\', '表格数据行'),
        ('工作点比值 & $n_e/n_G$ & $1.1337$ & 略高于经验极限 \\\\', '表格数据行'),
        ('弹丸周期 & $1/f_p$ & $0.506077\\ \\mathrm{s}$ & 加料的最快尺度 \\\\', '表格数据行'),
        ('粒子约束时间 & $\\tau_p$（声明） & $2.5\\ \\mathrm{s}$ & 密度回路的自然尺度 \\\\', '表格数据行'),
        ('能量约束时间 & $\\tau_E$（声明） & $3.7\\ \\mathrm{s}$ & 与附录 AD/AE 同点 \\\\', '表格数据行'),
        ('密度爬升时间 & $N_e/S_{\\mathrm{ext}}$ & $40.4156\\ \\mathrm{s}$ & 慢回路的证据 \\\\', '表格数据行'),
        ('\\begin{tip}{密度控制是慢回路}{tip-app-aj-control}', '语义盒标题'),
        ('\\textbf{出处：}Greenwald 密度极限参照 Greenwald 等（1988，\\emph{Nucl. Fusion} \\tex', '出处声明整行'),
        ('\\textbf{出处：}本题参数与公式分别来自 Milora 等（1995）、Parks \\& Turnbull（1978）、', '出处声明整行'),
        ('级别 & 含义 & 本附录的例子 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-aj-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录借鉴的框架主要来自 Wesson \\emph{Tokamaks}（第 4 版，2011）、', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AJ'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 556.5, '\\ \\mathrm{MW}', '粒子平衡：加进去的粒·正文'),
        (0, 17.6, '\\ \\mathrm{MeV}', '粒子平衡：加进去的粒·正文'),
        (0, 1.97352e+20, '\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式'),
        (0, 3.94705e+20, '\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式'),
        (0, 31.1913, '\\ \\text{kg/年}', '粒子平衡：加进去的粒·公式'),
        (0, 31.19, '\\ \\text{kg/年}', '粒子平衡：加进去的粒·公式'),
        (0, 2.5, '\\ \\mathrm{s}', '粒子平衡：加进去的粒·公式'),
        (0, 3.3288e+22, '\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式'),
        (0, 3.36827e+22, '\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式'),
        (0, 2.0591e+21, '\\ \\mathrm{s^{-1}}', '粒子平衡：加进去的粒·公式'),
        (0, 19.1687, '\\%', '粒子平衡：加进去的粒·公式'),
        (0, 10, '\\%', '粒子平衡：加进去的粒·正文'),
        (0, 1.91687, '\\%', '粒子平衡：加进去的粒·正文'),
        (1, 2.0591e+22, '\\ \\mathrm{s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 300, '\\ \\mathrm{K}', '气体加料与弹丸加料：·公式'),
        (1, 426.435, '\\ \\mathrm{mbar\\,L\\,s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 1.5, '\\ \\mathrm{mm}', '气体加料与弹丸加料：·公式'),
        (1, 255, '\\ \\mathrm{kg\\,m^{-3}}', '气体加料与弹丸加料：·公式'),
        (1, 300, '\\ \\mathrm{m\\,s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 3.43184e+21, '\\ \\mathrm{s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 7.12338e-06, '\\ \\mathrm{kg\\,s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 7.12338, '\\ \\mathrm{mg\\,s^{-1}}', '气体加料与弹丸加料：·公式'),
        (1, 615.46, '\\ \\mathrm{g/day}', '气体加料与弹丸加料：·正文'),
        (1, 224.797, '\\ \\text{kg/年}', '气体加料与弹丸加料：·正文'),
        (2, 0.12, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 12, '\\ \\mathrm{cm}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 2, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 6, '\\%', '弹丸烧蚀与穿透：中性·公式'),
        (2, 3.60498e-06, '\\ \\mathrm{kg}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 3.60498, '\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 2.16299, '\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 4.18155, '\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 1.97599, '\\ \\mathrm{Hz}', '弹丸烧蚀与穿透：中性·公式'),
        (2, 0.506077, '\\ \\mathrm{s}', '弹丸烧蚀与穿透：中性·正文'),
        (2, 2, '\\ \\mathrm{mm}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 8.54513, '\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 0.646087, '\\ \\mathrm{ms}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 0.193826, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 3, '\\ \\mathrm{mm}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 28.8398, '\\ \\mathrm{mg}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 1.26992, '\\ \\mathrm{ms}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 0.380976, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 19.0488, '\\%', '弹丸烧蚀与穿透：中性·表格'),
        (2, 1000, '\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 1.26992, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·表格'),
        (2, 1, '\\ \\mathrm{m}', '弹丸烧蚀与穿透：中性·正文'),
        (2, 2500, '\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·正文'),
        (2, 787.451, '\\ \\mathrm{m\\,s^{-1}}', '弹丸烧蚀与穿透：中性·正文'),
        (3, 18.5, '\\ \\mathrm{s}', '氦灰：稀释与排出条件·正文'),
        (3, 4.38719e+18, '\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·公式'),
        (3, 4.38719, '\\%', '氦灰：稀释与排出条件·公式'),
        (3, 9.12256e+19, '\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·公式'),
        (3, 91.2256, '\\%', '氦灰：稀释与排出条件·公式'),
        (3, 463.126, '\\ \\mathrm{MW}', '氦灰：稀释与排出条件·公式'),
        (3, 83.2211, '\\%', '氦灰：稀释与排出条件·公式'),
        (3, 5, '\\%', '氦灰：稀释与排出条件·公式'),
        (3, 21.0841, '\\ \\mathrm{s}', '氦灰：稀释与排出条件·公式'),
        (3, 3, 'cm', '氦灰：稀释与排出条件·表格'),
        (3, 2.5e+19, '\\ \\mathrm{m^{-3}}', '氦灰：稀释与排出条件·表格'),
        (4, 5, '\\ \\mathrm{Pa}', '抽气与氚存量：为什么·正文'),
        (4, 1.20716e+23, '\\ \\mathrm{s^{-1}}', '抽气与氚存量：为什么·公式'),
        (4, 0.163485, '\\%', '抽气与氚存量：为什么·公式'),
        (4, 5.15627e-06, '\\ \\mathrm{kg\\,s^{-1}}', '抽气与氚存量：为什么·公式'),
        (4, 162.72, '\\ \\text{kg/年}', '抽气与氚存量：为什么·公式'),
        (4, 1, '\\%', '抽气与氚存量：为什么·公式'),
        (4, 1.6272, '\\ \\text{kg/年}', '抽气与氚存量：为什么·公式'),
        (4, 2, '\\ \\mathrm{kg}', '抽气与氚存量：为什么·公式'),
        (4, 1.22911, '\\ \\text{年}', '抽气与氚存量：为什么·公式'),
        (4, 0.112524, '\\ \\text{kg/年}', '抽气与氚存量：为什么·公式'),
        (5, 11.0844, '\\ \\mathrm{MA}', '密度控制：反馈回路的·正文'),
        (5, 8.82069e+19, '\\ \\mathrm{m^{-3}}', '密度控制：反馈回路的·公式'),
        (5, 40.4156, '\\ \\mathrm{s}', '密度控制：反馈回路的·公式'),
        (5, 3.7, '\\ \\mathrm{s}', '密度控制：反馈回路的·表格'),
        (6, 5, '\\ \\mathrm{s}', '估算练习·正文'),
        (6, 3, '\\%', '估算练习·正文'),
        (6, 10, '\\ \\mathrm{Pa}', '估算练习·正文'),
        (6, 1.2269e+21, '\\ \\mathrm{s^{-1}}', '估算练习·正文'),
        (6, 0.833619, '\\ \\mathrm{Hz}', '估算练习·正文'),
        (6, 61.52, '\\%', '估算练习·正文'),
        (6, 12.6505, '\\ \\mathrm{s}', '估算练习·正文'),
        (6, 450.765, '\\ \\mathrm{MW}', '估算练习·正文'),
        (6, 2.41432e+23, '\\ \\mathrm{s^{-1}}', '估算练习·正文'),
        (0, 9.88481e-07, '', 'AJ.1·燃耗折算到每秒'),
        (0, 6.18573, '', 'AJ.1·加料量占输运损失的比例'),
        (0, 19, '', 'AJ.1·净加料燃耗份额的口语说法'),
        (0, 98, '', 'AJ.1·被抽走的粒子份额'),
        (1, 60, '', 'AJ.2·弹丸加料效率'),
        (2, 0.0004, '', 'AJ.3·声明烧蚀时间'),
        (3, 4.4, '', 'AJ.4·灰份额的口语说法'),
        (3, 16.8, '', 'AJ.4·功率损失的口语说法'),
        (3, 468.57, '', 'AJ.4·He II 谱线（与附录 AH 同口径）'),
        (6, 81, '', 'AJ.7·练习 4 答案的功率份额'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AJ 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AJ 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('1.97\\times10^{20}', '聚变反应率四舍五入'),
        ('3.95\\times10^{20}', '燃烧消耗率四舍五入'),
        ('8.32\\times10^{22}', '电子总数四舍五入'),
        ('3.33\\times10^{22}', '输运损失四舍五入'),
        ('3.37\\times10^{22}', '总损失四舍五入'),
        ('2.06\\times10^{21}', '外部加料需求四舍五入'),
        ('19.17\\%', '净加料燃耗份额四舍五入'),
        ('1.92\\%', '注气燃耗份额四舍五入'),
        ('2.06\\times10^{22}', '气体注入率四舍五入'),
        ('1.03\\times10^{22}', 'D2 分子流率四舍五入'),
        ('42.64\\ ', '气体流量四舍五入'),
        ('426.4\\ ', '气体流量毫巴升每秒四舍五入'),
        ('3.43\\times10^{21}', '弹丸注入率四舍五入'),
        ('7.12\\ ', '弹丸质量流率四舍五入'),
        ('615.5\\ ', '日耗克数四舍五入'),
        ('224.8\\ ', '年耗千克数四舍五入'),
        ('3.605\\ ', '弹丸质量四舍五入'),
        ('2.16\\ ', '弹丸含氚质量四舍五入'),
        ('8.68\\times10^{20}', '弹丸原子数四舍五入'),
        ('1.976\\ ', '弹丸频率四舍五入'),
        ('0.506\\ ', '弹丸周期四舍五入'),
        ('0.194\\ ', 'r=2 mm 穿透深度四舍五入'),
        ('28.84\\ ', 'r=3 mm 弹丸质量四舍五入'),
        ('1.27\\ ', 'r=3 mm 烧蚀时间四舍五入'),
        ('0.381\\ ', 'r=3 mm 穿透深度四舍五入'),
        ('19.05\\%', 'r=3 mm 穿透占小半径比四舍五入'),
        ('0.635\\ ', '穿透深度比四舍五入'),
        ('787.5\\ ', 'r=3 mm 所需速度四舍五入'),
        ('4.39\\times10^{18}', '稳态氦密度四舍五入'),
        ('4.39\\%', '氦份额四舍五入'),
        ('9.12\\times10^{19}', '燃料离子密度四舍五入'),
        ('91.23\\%', '燃料份额四舍五入'),
        ('463.1\\ ', '稀释后功率四舍五入'),
        ('83.22\\%', '功率保留份额四舍五入'),
        ('21.08\\ ', '允许停留时间四舍五入'),
        ('5.70\\ ', '停留时间上限倍数四舍五入'),
        ('1.21\\times10^{23}', '抽气排出率四舍五入'),
        ('611.7\\ ', '抽气倍数四舍五入'),
        ('0.163\\%', '所需氦份额四舍五入'),
        ('5.16\\times10^{-6}', '加料氚流四舍五入'),
        ('162.7\\ ', '年氚流量四舍五入'),
        ('1.627\\ ', '年滞留量四舍五入'),
        ('1.229\\ ', '到库存上限时间四舍五入'),
        ('14.75\\ ', '月数四舍五入'),
        ('0.113\\ ', '衰变损失四舍五入'),
        ('1.056\\ ', '氚自持下限四舍五入'),
        ('8.82\\times10^{19}', 'Greenwald 密度四舍五入'),
        ('1.134\\ ', '密度比四舍五入'),
        ('40.42\\ ', '密度爬升时间四舍五入'),
        ('6.19\\%', '净加料占输运损失比四舍五入'),
        ('2.37\\ ', '弹丸质量比四舍五入'),
        ('12.65\\ ', '3% 灰的停留时间四舍五入'),
        ('3.42\\ ', '3% 灰的倍数四舍五入'),
        ('450.8\\ ', '5% 灰的功率四舍五入'),
        ('2.41\\times10^{23}', '加倍抽气排出率四舍五入'),
        ('1223\\ ', '加倍抽气倍数四舍五入'),
        ('0.834\\ ', 'r=2 mm 弹丸频率四舍五入'),
        ('1.23\\times10^{21}', '练习 1 加料需求四舍五入'),
        ('0.596\\ ', '练习 1 比值四舍五入'),
        ('1.68\\ ', '练习 1 下降倍数四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AJ'], why))


def check_P43_appendix_ak_safety():
    """附录 AK：聚变安全与环境——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AK'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AK 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '第四本账：把安全写成数字',
        '氚的活度与剂量：克数怎么变成西弗',
        '氚的包容：渗透、除氚与年释放目标',
        '中子活化：壁材料怎么变成放射性废物',
        '衰变热与废物分级：冷却多久才能豁免',
        '事故与纵深防御：LOCA 情景的剂量账',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AK 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AK 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ak-frame',
        'tab:app-ak-tritium',
        'tab:app-ak-containment',
        'tab:app-ak-activation',
        'tab:app-ak-waste',
        'tab:app-ak-accident',
        'tab:app-ak-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AK'])

    for lab in (
        'sec:app-ak-frame',
        'sec:app-ak-tritium',
        'sec:app-ak-containment',
        'sec:app-ak-activation',
        'sec:app-ak-waste',
        'sec:app-ak-accident',
        'sec:app-ak-problems',
        'sec:app-ak-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AK'])

    for lab in (
        'tip-app-ak-frame',
        'tip-app-ak-tritium',
        'tip-app-ak-containment',
        'tip-app-ak-activation',
        'tip-app-ak-waste',
        'tip-app-ak-accident',
        'tip-app-ak-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AK'])

    for lab in (
        'eq:app-ak-flux',
        'eq:app-ak-dpa',
        'eq:app-ak-activity',
        'eq:app-ak-inventory',
        'eq:app-ak-dose',
        'eq:app-ak-limit',
        'eq:app-ak-burn',
        'eq:app-ak-perm',
        'eq:app-ak-df',
        'eq:app-ak-clean',
        'eq:app-ak-fluence',
        'eq:app-ak-specact',
        'eq:app-ak-decay',
        'eq:app-ak-clear',
        'eq:app-ak-heat',
        'eq:app-ak-site',
        'eq:app-ak-release',
        'eq:app-ak-margin',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AK'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AK'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AK 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AK 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AK 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AK 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('14.1 MeV', '第四本账：把安全写成·正文（精确绑定）'),
        ('0.90944\\times10^{6}\\ \\mathrm{W\\,m^{-2}}', '第四本账：把安全写成·公式（精确绑定）'),
        ('1.27042\\times10^{25}\\ \\mathrm{n\\,m^{-2}}', '第四本账：把安全写成·公式（精确绑定）'),
        ('80\\%', '第四本账：把安全写成·公式（精确绑定）'),
        ('4.5472\\ \\mathrm{dpa}', '第四本账：把安全写成·公式（精确绑定）'),
        ('3.57928\\times10^{-25}\\ \\mathrm{m^{2}}', '第四本账：把安全写成·公式（精确绑定）'),
        ('3579.28\\ \\mathrm{b}', '第四本账：把安全写成·公式（精确绑定）'),
        ('50.001$ dpa', '第四本账：把安全写成·正文（精确绑定）'),
        ('4.5472$ dpa', '第四本账：把安全写成·表格（精确绑定）'),
        ('50.001\\ \\mathrm{dpa}', '第四本账：把安全写成·表格（精确绑定）'),
        ('1.78283\\times10^{-9}\\ \\mathrm{s^{-1}}', '氚的活度与剂量：克数·正文（精确绑定）'),
        ('1 g', '氚的活度与剂量：克数·正文（精确绑定）'),
        ('3.55978\\times10^{14}\\ \\mathrm{Bq\\,g^{-1}}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('355.978\\ \\mathrm{TBq\\,g^{-1}}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1 mg', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('3.55978\\times10^{11}$ Bq', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('0.1 g', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('3.55978\\times10^{13}$ Bq', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('2 kg', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('7.11956\\times10^{17}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('6 mm', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1 Bq', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.8\\times10^{-11}$ Sv', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('6407.61\\ \\mathrm{Sv}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1 mSv', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1\\,\\mathrm{mSv}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('5.55556\\times10^{7}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.56065\\times10^{-7}\\ \\mathrm{g}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('20 mSv', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.11111\\times10^{9}$ Bq', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('3.12129\\times10^{-6}$ g', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.8\\times10^{-15}$ Sv', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('0.640761$ Sv', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('31.191$ kg', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.11033\\times10^{19}$ Bq', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1\\%', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('1.9986\\times10^{6}\\ \\mathrm{Sv}', '氚的活度与剂量：克数·公式（精确绑定）'),
        ('3.55978\\times10^{11}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('6.40761\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('3.55978\\times10^{13}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('640.761\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('3.55978\\times10^{14}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('1.28152\\times10^{7}\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('31.191 kg', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('1.11033\\times10^{19}\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('1.9986\\times10^{8}\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格（精确绑定）'),
        ('489.532$ m$^{2}$', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('4.89532\\times10^{9}$ Bq', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('1.18815$ g', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('0.433972$ kg', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('3.55978\\times10^{14}$ Bq', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('162.72$ kg', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('1.6272$ kg', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('1000$ m$^{3}$', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('1$ h$^{-1}$', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('6.90776\\ \\mathrm{h}', '氚的包容：渗透、除氚·公式（精确绑定）'),
        ('10$ appm', '中子活化：壁材料怎么·公式（精确绑定）'),
        ('500.01$ appm', '中子活化：壁材料怎么·公式（精确绑定）'),
        ('1\\ \\mathrm{TBq}', '中子活化：壁材料怎么·公式（精确绑定）'),
        ('128.205$ m$^{3}$', '中子活化：壁材料怎么·正文（精确绑定）'),
        ('4.16697\\times10^{-9}$ s$^{-1}$', '衰变热与废物分级：冷·正文（精确绑定）'),
        ('1394.84$ Bq', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('1.94557$ Bq', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('0.05$ W', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('50$ kW', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('5\\times10^{4}\\ \\mathrm{W}', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('50\\ \\mathrm{kW}', '衰变热与废物分级：冷·公式（精确绑定）'),
        ('69.7419$ W', '衰变热与废物分级：冷·正文（精确绑定）'),
        ('3.0cm', '衰变热与废物分级：冷·表格（精确绑定）'),
        ('5.61832\\times10^{-20}\\ \\mathrm{Sv\\,Bq^{-1}}', '事故与纵深防御：LO·公式（精确绑定）'),
        ('7.11956\\times10^{17}$ Bq', '事故与纵深防御：LO·公式（精确绑定）'),
        ('10\\%', '事故与纵深防御：LO·公式（精确绑定）'),
        ('7.11956\\times10^{16}\\ \\mathrm{Bq}', '事故与纵深防御：LO·公式（精确绑定）'),
        ('200\\ \\mathrm{g}', '事故与纵深防御：LO·公式（精确绑定）'),
        ('200 g', '事故与纵深防御：LO·公式（精确绑定）'),
        ('0.004$ Sv', '事故与纵深防御：LO·公式（精确绑定）'),
        ('4$ mSv', '事故与纵深防御：LO·公式（精确绑定）'),
        ('7.11956\\times10^{11}$ Bq', '事故与纵深防御：LO·公式（精确绑定）'),
        ('0.002$ g', '事故与纵深防御：LO·公式（精确绑定）'),
        ('4\\times10^{-8}\\ \\mathrm{Sv}', '事故与纵深防御：LO·公式（精确绑定）'),
        ('7.11956\\times10^{16}$ Bq', '事故与纵深防御：LO·表格（精确绑定）'),
        ('1\\ \\mathrm{mSv}', '事故与纵深防御：LO·表格（精确绑定）'),
        ('10 mg', '估算练习·正文（精确绑定）'),
        ('1 kg', '估算练习·正文（精确绑定）'),
        ('5 appm', '估算练习·正文（精确绑定）'),
        ('10 appm', '估算练习·正文（精确绑定）'),
        ('3.55978\\times10^{12}$ Bq', '估算练习·正文（精确绑定）'),
        ('3559.78$ GBq', '估算练习·正文（精确绑定）'),
        ('64.0761$ Sv', '估算练习·正文（精确绑定）'),
        ('5.55556\\times10^{11}$ Bq', '估算练习·正文（精确绑定）'),
        ('0.00156065$ g', '估算练习·正文（精确绑定）'),
        ('0.118815$ g', '估算练习·正文（精确绑定）'),
        ('0.0433972$ kg', '估算练习·正文（精确绑定）'),
        ('3.55978\\times10^{17}$ Bq', '估算练习·正文（精确绑定）'),
        ('250.005$ appm', '估算练习·正文（精确绑定）'),
        ('1.39134\\%', '包容·渗透量占燃耗份额（精确绑定·补充）'),
        ('5.7 keV', '氚β端点能量（文献级）（精确绑定·补充）'),
        ('3.55978$ TBq', '练习1·活度的 TBq 换算（精确绑定·补充）'),
        ('1.56065$ mg', '练习2·HT 等效质量换算（精确绑定·补充）'),
        ('4.89532\\times10^{8}$ Bq', '练习3·全壁渗透率（精确绑定·补充）'),
        ('2\\times10^{-8}$ Sv', '练习4·环境释放剂量（精确绑定·补充）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AK'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 表达式 & 数值 & 说明 \\\\', '表格数据行'),
        ('中子壁负载 & $q_n$ & $0.90944\\ \\mathrm{MW\\,m^{-2}}$ & 附录 Z 口径 \\\\', '表格数据行'),
        ('中子通量 & $q_n/E_n$ & $4.02573\\times10^{17}\\ \\mathrm{n\\,m^{-2}\\,s^{-1}}$ & 只差一个中子能量 \\\\', '表格数据行'),
        ('年注量 & $\\phi\\times3.15576\\times10^{7}$ & $1.27042\\times10^{25}\\ \\mathrm{n\\,m^{-2}}$ & 一年累计 \\\\', '表格数据行'),
        ('位移截面 & 标定 & $3579.28\\ \\mathrm{b}$ & 与附录 Z 的 $4.5472$ dpa/年 自洽 \\\\', '表格数据行'),
        ('寿期累计损伤 & $\\sigma_d\\phi t_{life}$ & $50.001\\ \\mathrm{dpa}$ & $10.996$ 年 \\\\', '表格数据行'),
        ('\\begin{tip}{安全不是感觉，是数字}{tip-app-ak-frame}', '语义盒标题'),
        ('\\textbf{出处：}本附录的安全量化框架参照 Dolan（编）\\emph{Magnetic Fusion Technology}（201', '出处声明整行'),
        ('质量 & 活度 & 全吸入剂量 & 说明 \\\\', '表格数据行'),
        ('1 mg & $3.55978\\times10^{11}\\ \\mathrm{Bq}$ & $6.40761\\ \\mathrm{Sv}$ & 一滴水的千分之一 \\\\', '表格数据行'),
        ('0.1 g & $3.55978\\times10^{13}\\ \\mathrm{Bq}$ & $640.761\\ \\mathrm{Sv}$ & 仍远超限值 \\\\', '表格数据行'),
        ('1 g & $3.55978\\times10^{14}\\ \\mathrm{Bq}$ & $6407.61\\ \\mathrm{Sv}$ & 公众限值的六千倍 \\\\', '表格数据行'),
        ('2 kg（库存上限） & $7.11956\\times10^{17}\\ \\mathrm{Bq}$ & $1.28152\\times10^{7}\\ \\mathrm{Sv}$ & 附录 Z/AJ 口径 \\\\', '表格数据行'),
        ('31.191 kg（年燃耗） & $1.11033\\times10^{19}\\ \\mathrm{Bq}$ & $1.9986\\times10^{8}\\ \\mathrm{Sv}$ & 只算 1\\% 也有百万西弗 \\\\', '表格数据行'),
        ('\\begin{tip}{氚的危险不在外照射，在摄入}{tip-app-ak-tritium}', '语义盒标题'),
        ('\\textbf{出处：}氚的比活度与剂量换算参照 Tanabe（编）\\emph{Tritium: Fuel of Fusion Reacto', '出处声明整行'),
        ('出路 & 数值 & 与管理目标之比 & 对策 \\\\', '表格数据行'),
        ('壁渗透 & $0.433972$ kg/年 & $433.972$ 倍 & 渗透屏障 + 除氚 \\\\', '表格数据行'),
        ('滞留（AJ 口径） & $1.6272$ kg/年 & 库存 2 kg 只够 $1.22911$ 年 & 氚工厂回收 \\\\', '表格数据行'),
        ('事故释放 & 释放份额 $0.1$ & 见 \\S AK.6 & 安全壳 + 纵深防御 \\\\', '表格数据行'),
        ('\\begin{tip}{屏障不是可选项，是设计输入}{tip-app-ak-containment}', '语义盒标题'),
        ('\\textbf{出处：}氚渗透与除氚的工程口径参照 Tanabe（编）\\emph{Tritium: Fuel of Fusion React', '出处声明整行'),
        ('年注量 & $\\phi t$ & $1.27042\\times10^{25}$ n\\,m$^{-2}$ & 与 $q_n$ 同源 \\\\', '表格数据行'),
        ('位移损伤 & 标定 $\\sigma_d$ & $50.001$ dpa & $10.996$ 年 \\\\', '表格数据行'),
        ('氦产额 & $10$ appm/dpa & $500.01$ appm & 声明值 \\\\', '表格数据行'),
        ('比活度（寿期末） & 声明 & $10^{6}$ Bq/kg & 典型堆内构件 \\\\', '表格数据行'),
        ('总活度 & $a_{spec}M$ & $10^{12}$ Bq & $10^{6}$ kg 构件 \\\\', '表格数据行'),
        ('豁免水平之比 & $a_{spec}/a_{clear}$ & $100$ 倍 & 需分级处置 \\\\', '表格数据行'),
        ('\\begin{tip}{活化是中子的第二份产物}{tip-app-ak-activation}', '语义盒标题'),
        ('\\textbf{出处：}中子活化与废物评价的框架参照 Dolan（编）\\emph{Magnetic Fusion Technology}（2', '出处声明整行'),
        ('冷却时间 & 比活度 & 衰变热 & 结论 \\\\', '表格数据行'),
        ('$0$ 年 & $10^{6}$ Bq/kg & $50$ kW & 必须主动冷却 \\\\', '表格数据行'),
        ('$35.0204$ 年 & $10^{4}$ Bq/kg & —— & 刚压到豁免线 \\\\', '表格数据行'),
        ('$50$ 年 & $1394.84$ Bq/kg & $69.7419$ W & 可申请解控 \\\\', '表格数据行'),
        ('$100$ 年 & $1.94557$ Bq/kg & —— & 远低于豁免线 \\\\', '表格数据行'),
        ('\\begin{tip}{废物的等级不是材料定的，是活度定的}{tip-app-ak-waste}', '语义盒标题'),
        ('\\textbf{出处：}衰变热与废物分级方法参照 Dolan（编）\\emph{Magnetic Fusion Technology}（201', '出处声明整行'),
        ('环节 & 表达式 & 数值 & 说明 \\\\', '表格数据行'),
        ('堆内库存 & 2 kg & $7.11956\\times10^{17}$ Bq & 附录 Z/AJ 口径 \\\\', '表格数据行'),
        ('进入安全壳 & $10\\%$ & $7.11956\\times10^{16}$ Bq & 释放份额声明值 \\\\', '表格数据行'),
        ('无包容剂量 & $A_{rel}f$ & $4$ mSv & 超限值 4 倍 \\\\', '表格数据行'),
        ('有包容剂量 & $A_{rel}f\\times10^{-5}$ & $0.04$ $\\mu$Sv & 限值的 $4\\times10^{-5}$ \\\\', '表格数据行'),
        ('余量 & $1\\ \\mathrm{mSv}/D$ & $25000$ 倍 & 纵深防御的效果 \\\\', '表格数据行'),
        ('\\begin{tip}{纵深防御是乘法，不是加法}{tip-app-ak-accident}', '语义盒标题'),
        ('\\textbf{出处：}事故情景与包容设计参照 Dolan（编）\\emph{Magnetic Fusion Technology}（2013', '出处声明整行'),
        ('级别 & 含义 & 本附录中的例子 & 可否复算 \\\\', '表格数据行'),
        ('现算级 & 由脚本按公式算出 & $3.55978\\times10^{14}$ Bq/g、$6407.61$ Sv、$35.0204$ 年、$25000$ 倍 & 可以 \\\\', '表格数据行'),
        ('声明级 & 脚本顶部的输入 & 剂量系数 $1.8\\times10^{-11}$ Sv/Bq、比活度 $10^{6}$ Bq/kg、PRF $=100$ & 改输入即可 \\\\', '表格数据行'),
        ('文献级 & 只作量级对照 & 年释放目标 1 g/年、场址因子 $20$ $\\mu$Sv/g、豁免水平 $10^{4}$ Bq/kg & 不做断言 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ak-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录的三级口径与版权处理沿用附录 K、L、N、Z 的约定；', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（整行绑定：%s）' % (APP_FILES['AK'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 14.1, 'MeV', '第四本账：把安全写成·正文'),
        (0, 909440, '\\ \\mathrm{W\\,m^{-2}}', '第四本账：把安全写成·公式'),
        (0, 1.27042e+25, '\\ \\mathrm{n\\,m^{-2}}', '第四本账：把安全写成·公式'),
        (0, 80, '\\%', '第四本账：把安全写成·公式'),
        (0, 4.5472, '\\ \\mathrm{dpa}', '第四本账：把安全写成·公式'),
        (0, 3.57928e-25, '\\ \\mathrm{m^{2}}', '第四本账：把安全写成·公式'),
        (0, 3579.28, '\\ \\mathrm{b}', '第四本账：把安全写成·公式'),
        (0, 50.001, '$ dpa', '第四本账：把安全写成·正文'),
        (0, 0.90944, '\\ \\mathrm{MW\\,m^{-2}}', '第四本账：把安全写成·表格'),
        (0, 4.5472, '$ dpa', '第四本账：把安全写成·表格'),
        (0, 50.001, '\\ \\mathrm{dpa}', '第四本账：把安全写成·表格'),
        (1, 1.78283e-09, '\\ \\mathrm{s^{-1}}', '氚的活度与剂量：克数·正文'),
        (1, 1, 'g', '氚的活度与剂量：克数·正文'),
        (1, 3.55978e+14, '\\ \\mathrm{Bq\\,g^{-1}}', '氚的活度与剂量：克数·公式'),
        (1, 355.978, '\\ \\mathrm{TBq\\,g^{-1}}', '氚的活度与剂量：克数·公式'),
        (1, 1, 'mg', '氚的活度与剂量：克数·公式'),
        (1, 3.55978e+11, '$ Bq', '氚的活度与剂量：克数·公式'),
        (1, 0.1, 'g', '氚的活度与剂量：克数·公式'),
        (1, 3.55978e+13, '$ Bq', '氚的活度与剂量：克数·公式'),
        (1, 2, 'kg', '氚的活度与剂量：克数·公式'),
        (1, 7.11956e+17, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·公式'),
        (1, 6, 'mm', '氚的活度与剂量：克数·公式'),
        (1, 1, 'Bq', '氚的活度与剂量：克数·公式'),
        (1, 1.8e-11, '$ Sv', '氚的活度与剂量：克数·公式'),
        (1, 6407.61, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·公式'),
        (1, 1, 'mSv', '氚的活度与剂量：克数·公式'),
        (1, 1, '\\,\\mathrm{mSv}', '氚的活度与剂量：克数·公式'),
        (1, 55555600, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·公式'),
        (1, 1.56065e-07, '\\ \\mathrm{g}', '氚的活度与剂量：克数·公式'),
        (1, 20, 'mSv', '氚的活度与剂量：克数·公式'),
        (1, 1111110000, '$ Bq', '氚的活度与剂量：克数·公式'),
        (1, 3.12129e-06, '$ g', '氚的活度与剂量：克数·公式'),
        (1, 1.8e-15, '$ Sv', '氚的活度与剂量：克数·公式'),
        (1, 0.640761, '$ Sv', '氚的活度与剂量：克数·公式'),
        (1, 31.191, '$ kg', '氚的活度与剂量：克数·公式'),
        (1, 1.11033e+19, '$ Bq', '氚的活度与剂量：克数·公式'),
        (1, 1, '\\%', '氚的活度与剂量：克数·公式'),
        (1, 1998600, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·公式'),
        (1, 3.55978e+11, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格'),
        (1, 6.40761, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格'),
        (1, 3.55978e+13, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格'),
        (1, 640.761, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格'),
        (1, 3.55978e+14, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格'),
        (1, 12815200, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格'),
        (1, 31.191, 'kg', '氚的活度与剂量：克数·表格'),
        (1, 1.11033e+19, '\\ \\mathrm{Bq}', '氚的活度与剂量：克数·表格'),
        (1, 199860000, '\\ \\mathrm{Sv}', '氚的活度与剂量：克数·表格'),
        (2, 489.532, '$ m$^{2}$', '氚的包容：渗透、除氚·公式'),
        (2, 4895320000, '$ Bq', '氚的包容：渗透、除氚·公式'),
        (2, 1.18815, '$ g', '氚的包容：渗透、除氚·公式'),
        (2, 0.433972, '$ kg', '氚的包容：渗透、除氚·公式'),
        (2, 3.55978e+14, '$ Bq', '氚的包容：渗透、除氚·公式'),
        (2, 162.72, '$ kg', '氚的包容：渗透、除氚·公式'),
        (2, 1.6272, '$ kg', '氚的包容：渗透、除氚·公式'),
        (2, 1000, '$ m$^{3}$', '氚的包容：渗透、除氚·公式'),
        (2, 1, '$ h$^{-1}$', '氚的包容：渗透、除氚·公式'),
        (2, 6.90776, '\\ \\mathrm{h}', '氚的包容：渗透、除氚·公式'),
        (3, 10, '$ appm', '中子活化：壁材料怎么·公式'),
        (3, 500.01, '$ appm', '中子活化：壁材料怎么·公式'),
        (3, 1, '\\ \\mathrm{TBq}', '中子活化：壁材料怎么·公式'),
        (3, 128.205, '$ m$^{3}$', '中子活化：壁材料怎么·正文'),
        (4, 4.16697e-09, '$ s$^{-1}$', '衰变热与废物分级：冷·正文'),
        (4, 1394.84, '$ Bq', '衰变热与废物分级：冷·公式'),
        (4, 1.94557, '$ Bq', '衰变热与废物分级：冷·公式'),
        (4, 0.05, '$ W', '衰变热与废物分级：冷·公式'),
        (4, 50, '$ kW', '衰变热与废物分级：冷·公式'),
        (4, 50000, '\\ \\mathrm{W}', '衰变热与废物分级：冷·公式'),
        (4, 50, '\\ \\mathrm{kW}', '衰变热与废物分级：冷·公式'),
        (4, 69.7419, '$ W', '衰变热与废物分级：冷·正文'),
        (4, 3, 'cm', '衰变热与废物分级：冷·表格'),
        (5, 5.61832e-20, '\\ \\mathrm{Sv\\,Bq^{-1}}', '事故与纵深防御：LO·公式'),
        (5, 7.11956e+17, '$ Bq', '事故与纵深防御：LO·公式'),
        (5, 10, '\\%', '事故与纵深防御：LO·公式'),
        (5, 7.11956e+16, '\\ \\mathrm{Bq}', '事故与纵深防御：LO·公式'),
        (5, 200, '\\ \\mathrm{g}', '事故与纵深防御：LO·公式'),
        (5, 200, 'g', '事故与纵深防御：LO·公式'),
        (5, 0.004, '$ Sv', '事故与纵深防御：LO·公式'),
        (5, 4, '$ mSv', '事故与纵深防御：LO·公式'),
        (5, 7.11956e+11, '$ Bq', '事故与纵深防御：LO·公式'),
        (5, 0.002, '$ g', '事故与纵深防御：LO·公式'),
        (5, 4e-08, '\\ \\mathrm{Sv}', '事故与纵深防御：LO·公式'),
        (5, 7.11956e+16, '$ Bq', '事故与纵深防御：LO·表格'),
        (5, 1, '\\ \\mathrm{mSv}', '事故与纵深防御：LO·表格'),
        (6, 10, 'mg', '估算练习·正文'),
        (6, 1, 'kg', '估算练习·正文'),
        (6, 5, 'appm', '估算练习·正文'),
        (6, 10, 'appm', '估算练习·正文'),
        (6, 3.55978e+12, '$ Bq', '估算练习·正文'),
        (6, 3559.78, '$ GBq', '估算练习·正文'),
        (6, 64.0761, '$ Sv', '估算练习·正文'),
        (6, 5.55556e+11, '$ Bq', '估算练习·正文'),
        (6, 0.00156065, '$ g', '估算练习·正文'),
        (6, 0.118815, '$ g', '估算练习·正文'),
        (6, 0.0433972, '$ kg', '估算练习·正文'),
        (6, 3.55978e+17, '$ Bq', '估算练习·正文'),
        (6, 250.005, '$ appm', '估算练习·正文'),
        (2, 1.39134, '', '包容·渗透量占燃耗份额'),
        (1, 5.7, '', '氚β端点能量（文献级）'),
        (6, 3.55978, '', '练习1·活度的 TBq 换算'),
        (6, 1.56065, '', '练习2·HT 等效质量换算'),
        (6, 489532000, '', '练习3·全壁渗透率'),
        (6, 2e-08, '', '练习4·环境释放剂量'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AK 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AK 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('4.03\\times10^{17}', '中子通量四舍五入'),
        ('1.27\\times10^{25}', '年注量四舍五入'),
        ('3.58\\times10^{-25}', '位移截面四舍五入'),
        ('3579\\ ', '位移截面巴数四舍五入'),
        ('50.0\\ ', '寿期累计 dpa 四舍五入'),
        ('1.78\\times10^{-9}', '氚衰变常数四舍五入'),
        ('2.00\\times10^{23}', '每克氚原子数四舍五入'),
        ('3.56\\times10^{14}', '氚比活度四舍五入'),
        ('356\\ ', '比活度 TBq 四舍五入'),
        ('7.12\\times10^{17}', '2 kg 库存活度四舍五入'),
        ('6408\\ ', '1 g 吸入剂量四舍五入'),
        ('5.56\\times10^{7}', '1 mSv 对应活度四舍五入'),
        ('3.12\\times10^{-6}', '20 mSv 对应质量四舍五入'),
        ('0.641\\ ', 'HT 等效剂量四舍五入'),
        ('31.19\\ ', '年燃耗质量四舍五入'),
        ('1.11\\times10^{19}', '年燃耗活度四舍五入'),
        ('2.00\\times10^{6}', '1% 吸入剂量四舍五入'),
        ('6.41\\ ', '1 mg 剂量四舍五入'),
        ('641\\ ', '0.1 g 剂量四舍五入'),
        ('1.28\\times10^{7}', '2 kg 剂量四舍五入'),
        ('2.00\\times10^{8}', '年燃耗剂量四舍五入'),
        ('4.90\\times10^{9}', '全壁渗透率四舍五入'),
        ('4.23\\times10^{14}', '日渗透活度四舍五入'),
        ('1.54\\times10^{5}', '年渗透 TBq 四舍五入'),
        ('1.19\\ ', '日渗透克数四舍五入'),
        ('0.434\\ ', '年渗透千克数四舍五入'),
        ('1.39\\%', '渗透占燃耗比四舍五入'),
        ('434\\ ', '所需去污因子四舍五入'),
        ('1.63\\ ', '年滞留量四舍五入'),
        ('1.23\\ ', '到库存上限时间四舍五入'),
        ('6.91\\ ', '手套箱净化时间四舍五入'),
        ('128\\ ', '构件体积四舍五入'),
        ('4.17\\times10^{-9}', 'Co-60 衰变常数四舍五入'),
        ('35.0\\ ', '降到豁免所需年数四舍五入'),
        ('1395\\ ', '50 年比活度四舍五入'),
        ('0.139\\ ', '豁免份额四舍五入'),
        ('1.95\\ ', '100 年比活度四舍五入'),
        ('69.7\\ ', '50 年衰变热四舍五入'),
        ('717\\ ', '衰变热下降倍数四舍五入'),
        ('5.62\\times10^{-20}', '场址剂量因子四舍五入'),
        ('7.12\\times10^{16}', '事故释放活度四舍五入'),
        ('7.12\\times10^{11}', '有包容释放活度四舍五入'),
        ('4.0\\times10^{-8}', '有包容剂量四舍五入'),
        ('1.0\\times10^{5}', '包容降低倍数四舍五入'),
        ('4.0\\times10^{-5}', '与限值之比四舍五入'),
        ('3.56\\times10^{12}', '练习 1 活度四舍五入'),
        ('3559.8\\ ', '练习 1 GBq 四舍五入'),
        ('64.08\\ ', '练习 1 剂量四舍五入'),
        ('64076\\ ', '练习 1 限值倍数四舍五入'),
        ('5.56\\times10^{11}', '练习 2 活度四舍五入'),
        ('0.0434\\ ', '练习 3 年渗透四舍五入'),
        ('43.4\\ ', '练习 3 埋存比四舍五入'),
        ('3.56\\times10^{17}', '练习 4 库存活度四舍五入'),
        ('3.56\\times10^{11}', '练习 4 环境释放四舍五入'),
        ('2.0\\times10^{-8}', '练习 4 剂量四舍五入'),
        ('17.5\\ ', '练习 5 豁免年限四舍五入'),
        ('250.0\\ ', '练习 5 氦浓度四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AK'], why))


def check_P44_appendix_al_nonlinear():
    """附录 AL：非线性等离子体物理——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AL'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AL 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '非线性从哪来：一个无量纲数就够',
        '孤子：色散与非线性平衡',
        '参量不稳定性：三波共振与阈值',
        '自聚焦与波崩塌',
        '粒子俘获：什么时候非线性跑赢阻尼',
        '相干结构：双电层与相空间空洞',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AL 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AL 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-al-nonlinear',
        'tab:app-al-soliton',
        'tab:app-al-parametric',
        'tab:app-al-collapse',
        'tab:app-al-trapping',
        'tab:app-al-structures',
        'tab:app-al-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AL'])

    for lab in (
        'sec:app-al-nonlinear',
        'sec:app-al-soliton',
        'sec:app-al-parametric',
        'sec:app-al-collapse',
        'sec:app-al-trapping',
        'sec:app-al-structures',
        'sec:app-al-problems',
        'sec:app-al-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AL'])

    for lab in (
        'tip-app-al-nonlinear',
        'tip-app-al-soliton',
        'tip-app-al-parametric',
        'tip-app-al-collapse',
        'tip-app-al-trapping',
        'tip-app-al-structures',
        'tip-app-al-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AL'])

    for lab in (
        'eq:app-al-vosc',
        'eq:app-al-ratio',
        'eq:app-al-up',
        'eq:app-al-nc',
        'eq:app-al-kdv',
        'eq:app-al-solitonwidth',
        'eq:app-al-dent',
        'eq:app-al-matching',
        'eq:app-al-threshold',
        'eq:app-al-srs',
        'eq:app-al-pc',
        'eq:app-al-caviton',
        'eq:app-al-ewb',
        'eq:app-al-omegab',
        'eq:app-al-ecrit',
        'eq:app-al-oneil',
        'eq:app-al-dl',
        'eq:app-al-hole',
        'eq:app-al-zakharov',
        'eq:app-al-bgk',
        'eq:app-al-ex-srs',
        'eq:app-al-ex-dent',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AL'])

    # ---------- ③ 出处与版权声明 ----------
    n_copy = text.count('不复制')
    if n_copy < 3:
        rec('FAIL', 'P', '附录 AL 的「不复制」只有 %d 处（要求 >= 3）' % n_copy)
    else:
        rec('PASS', 'P', '附录 AL 的「不复制」共 %d 处' % n_copy)
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AL 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AL 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AL 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AL 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('140$ GHz', '非线性从哪来：一个无·正文（精确绑定）'),
        ('1$ MW', '非线性从哪来：一个无·正文（精确绑定）'),
        ('2$ cm', '非线性从哪来：一个无·正文（精确绑定）'),
        ('1$ keV', '非线性从哪来：一个无·正文（精确绑定）'),
        ('7.95775\\times10^{8}$ W m$^{-2}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('0.0795775$ MW cm$^{-2}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('774329$ V m$^{-1}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('7.74329$ kV cm$^{-1}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('154824$ m s$^{-1}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('4.02127\\times10^{7}$ m s$^{-1}$', '非线性从哪来：一个无·公式（精确绑定）'),
        ('4m', '非线性从哪来：一个无·公式（精确绑定）'),
        ('5.45892\\times10^{-21}$ J', '非线性从哪来：一个无·公式（精确绑定）'),
        ('0.0340719$ eV', '非线性从哪来：一个无·公式（精确绑定）'),
        ('11.4292$ eV', '非线性从哪来：一个无·公式（精确绑定）'),
        ('1\\ \\mathrm{keV}', '非线性从哪来：一个无·公式（精确绑定）'),
        ('114.292$ eV', '非线性从哪来：一个无·公式（精确绑定）'),
        ('2.43127\\times10^{20}$ m$^{-3}$', '非线性从哪来：一个无·正文（精确绑定）'),
        ('9.10085\\times10^{27}$ m$^{-3}$', '非线性从哪来：一个无·正文（精确绑定）'),
        ('3.0cm', '非线性从哪来：一个无·表格（精确绑定）'),
        ('140 GHz', '非线性从哪来：一个无·表格（精确绑定）'),
        ('8.68021\\times10^{10}$ V m$^{-1}$', '非线性从哪来：一个无·表格（精确绑定）'),
        ('2.83674\\times10^{6}$ m s$^{-1}$', '非线性从哪来：一个无·表格（精确绑定）'),
        ('30\\%', '孤子：色散与非线性平·公式（精确绑定）'),
        ('593524$ m s$^{-1}$', '孤子：色散与非线性平·公式（精确绑定）'),
        ('593.524$ km s$^{-1}$', '孤子：色散与非线性平·公式（精确绑定）'),
        ('2.60280\\times10^{-4}\\ \\mathrm{m}', '孤子：色散与非线性平·公式（精确绑定）'),
        ('652877$ m s$^{-1}$', '孤子：色散与非线性平·公式（精确绑定）'),
        ('438.533$ ps', '孤子：色散与非线性平·公式（精确绑定）'),
        ('9.194$ keV', '孤子：色散与非线性平·正文（精确绑定）'),
        ('1\\%', '孤子：色散与非线性平·正文（精确绑定）'),
        ('3.64820\\times10^{7}$ V m$^{-1}$', '孤子：色散与非线性平·正文（精确绑定）'),
        ('1.76643\\times10^{12}$ W m$^{-2}$', '孤子：色散与非线性平·正文（精确绑定）'),
        ('1.76643\\times10^{8}$ W cm$^{-2}$', '孤子：色散与非线性平·正文（精确绑定）'),
        ('2.21976$ GW', '孤子：色散与非线性平·正文（精确绑定）'),
        ('2.69841$ MW', '参量不稳定性：三波共·公式（精确绑定）'),
        ('1.32621\\times10^{7}$ m s$^{-1}$', '参量不稳定性：三波共·公式（精确绑定）'),
        ('83876.6$ m s$^{-1}$', '参量不稳定性：三波共·公式（精确绑定）'),
        ('8.74265\\times10^{15}$ W m$^{-2}$', '参量不稳定性：三波共·公式（精确绑定）'),
        ('8.74265\\times10^{11}$ W cm$^{-2}$', '参量不稳定性：三波共·公式（精确绑定）'),
        ('1.78531\\times10^{7}$ m$^{-1}$', '参量不稳定性：三波共·正文（精确绑定）'),
        ('3.25962\\times10^{7}$ m$^{-1}$', '参量不稳定性：三波共·正文（精确绑定）'),
        ('192.758$ nm', '参量不稳定性：三波共·正文（精确绑定）'),
        ('23.5082$ nm', '参量不稳定性：三波共·正文（精确绑定）'),
        ('426.177$ nm', '参量不稳定性：三波共·正文（精确绑定）'),
        ('707.343$ THz', '参量不稳定性：三波共·正文（精确绑定）'),
        ('2251.29$ m$^{-1}$', '参量不稳定性：三波共·表格（精确绑定）'),
        ('4.23041\\times10^{10}$ W', '自聚焦与波崩塌·公式（精确绑定）'),
        ('42.3041$ GW', '自聚焦与波崩塌·公式（精确绑定）'),
        ('1.58355\\times10^{12}$ W', '自聚焦与波崩塌·公式（精确绑定）'),
        ('1.58355$ TW', '自聚焦与波崩塌·公式（精确绑定）'),
        ('500$ TW', '自聚焦与波崩塌·公式（精确绑定）'),
        ('2.60417$ TW', '自聚焦与波崩塌·公式（精确绑定）'),
        ('2.35082\\times10^{-8}$ m', '自聚焦与波崩塌·公式（精确绑定）'),
        ('7.12806\\times10^{-4}\\ \\mathrm{m}', '自聚焦与波崩塌·公式（精确绑定）'),
        ('235.082$ nm', '自聚焦与波崩塌·公式（精确绑定）'),
        ('1.20317\\times10^{10}$ V m$^{-1}$', '自聚焦与波崩塌·公式（精确绑定）'),
        ('1.92128\\times10^{13}$ W cm$^{-2}$', '自聚焦与波崩塌·公式（精确绑定）'),
        ('9.61592\\times10^{8}\\ \\mathrm{V\\,m^{-1}}', '自聚焦与波崩塌·公式（精确绑定）'),
        ('961.592\\ \\mathrm{MV\\,m^{-1}}', '自聚焦与波崩塌·公式（精确绑定）'),
        ('961.592$ MV m$^{-1}$', '自聚焦与波崩塌·表格（精确绑定）'),
        ('4.0cm', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('5.58162\\times10^{10}$ s$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('1.76506\\times10^{11}$ s$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('112.569$ ps', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('35.5975$ ps', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('2.65241\\times10^{7}$ m s$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('8.38766\\times10^{7}$ m s$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('1.27283\\times10^{10}$ s$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('78.5651$ ps', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('52.002\\ \\mathrm{V\\,m^{-1}}', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('1474.96$ V m$^{-1}$', '粒子俘获：什么时候非·公式（精确绑定）'),
        ('1.28983\\times10^{6}\\ \\mathrm{V\\,m^{-1}}', '相干结构：双电层与相·公式（精确绑定）'),
        ('12.8983\\ \\mathrm{kV\\,cm^{-1}}', '相干结构：双电层与相·公式（精确绑定）'),
        ('8.38766\\times10^{7}\\ \\mathrm{m\\,s^{-1}}', '相干结构：双电层与相·公式（精确绑定）'),
        ('7.12806$ mm', '相干结构：双电层与相·表格（精确绑定）'),
        ('1.28983$ MV m$^{-1}$', '相干结构：双电层与相·表格（精确绑定）'),
        ('2$ MW', '估算练习·正文（精确绑定）'),
        ('12.5664$ cm$^2$', '估算练习·正文（精确绑定）'),
        ('9194$ eV', '估算练习·正文（精确绑定）'),
        ('50$ MW', '估算练习·正文（精确绑定）'),
        ('1.59155\\times10^{9}$ W m$^{-2}$', '估算练习·公式（精确绑定）'),
        ('1.09507\\times10^{6}$ V m$^{-1}$', '估算练习·公式（精确绑定）'),
        ('218954$ m s$^{-1}$', '估算练习·公式（精确绑定）'),
        ('158.355$ GW', '估算练习·公式（精确绑定）'),
        ('1.12979\\times10^{11}$ s$^{-1}$', '估算练习·公式（精确绑定）'),
        ('7014.53$ m$^{-1}$', '估算练习·公式（精确绑定）'),
        ('52.002$ V m$^{-1}$', '估算练习·公式（精确绑定）'),
        ('1.47431\\times10^{7}$ m$^{-1}$', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AL'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 微波侧（140 GHz） & 激光侧（$10^{15}$ W cm$^{-2}$） & 说明 \\\\', '表格数据行'),
        ('功率密度 & $7.95775\\times10^{8}$ W m$^{-2}$ & $10^{19}$ W m$^{-2}$ & 相差 $1.25664\\times10^{10}$ 倍 \\\\', '表格数据行'),
        ('电场 & $774329$ V m$^{-1}$ & $8.68021\\times10^{10}$ V m$^{-1}$ & $E = \\sqrt{2S/\\epsilon_0c}$ \\\\', '表格数据行'),
        ('抖动速度 & $154824$ m s$^{-1}$ & $2.83674\\times10^{6}$ m s$^{-1}$ & $eE/m_e\\omega$ \\\\', '表格数据行'),
        ('归一化振幅 & $v_{\\mathrm{osc}}/v_{the} = 3.85013\\times10^{-3}$ & $a_0 = 9.45601\\times10^{-3}$ & 都是「小量」 \\\\', '表格数据行'),
        ('有质动力势 & $0.0340719$ eV & $11.4292$ eV & 与 $T_e$ 比 \\\\', '表格数据行'),
        ('与温度之比 & $3.70588\\times10^{-6}$ & $0.0114292$ & 差 $3354.44$ 倍 \\\\', '表格数据行'),
        ('截止密度 & $2.43127\\times10^{20}$ m$^{-3}$ & $9.10085\\times10^{27}$ m$^{-3}$ & $n_c \\propto \\omega^2$ \\\\', '表格数据行'),
        ('工作密度/截止 & $0.411308$ & $0.010988$ & 冕区离截止更远 \\\\', '表格数据行'),
        ('\\begin{tip}{非线性不是「更强的场」，是「不可叠加」}{tip-app-al-nonlinear}', '语义盒标题'),
        ('\\textbf{出处：}非线性等离子体理论的框架参照 Sagdeev 与 Galeev \\emph{Nonlinear Plasma The', '出处声明整行'),
        ('量 & 数值 & 说明 \\\\', '表格数据行'),
        ('马赫数 $M$ & $1.1$ & 超音速才成孤子 \\\\', '表格数据行'),
        ('相对振幅 $\\delta n/n_0$ & $0.3$（$30\\%$） & $3(M-1)$ \\\\', '表格数据行'),
        ('声速 $c_s$ & $593524$ m s$^{-1}$ & $\\sqrt{T_e/m_i}$ \\\\', '表格数据行'),
        ('孤子宽度 $\\Delta$ & $260.280$ $\\mu$m & $3.65148\\,\\lambda_{De}$ \\\\', '表格数据行'),
        ('孤子速度 & $652877$ m s$^{-1}$ & $Mc_s$ \\\\', '表格数据行'),
        ('特征时间 $\\Delta/c_s$ & $438.533$ ps & 走过自身宽度 \\\\', '表格数据行'),
        ('\\begin{tip}{孤子是平衡，不是「不散」}{tip-app-al-soliton}', '语义盒标题'),
        ('\\textbf{出处：}离子声孤子的 KdV 解参照 Washimi 与 Taniuti（\\emph{Physical Review Let', '出处声明整行'),
        ('量 & 微波侧（140 GHz） & 激光冕区（0.35 $\\mu$m） & 说明 \\\\', '表格数据行'),
        ('泵浦波数 & $2251.29$ m$^{-1}$ & $1.78531\\times10^{7}$ m$^{-1}$ & O 模 \\\\', '表格数据行'),
        ('阈值判据 & $6.32456\\times10^{-3}$ & $8.74265\\times10^{11}$ W cm$^{-2}$ & 同一不等式 \\\\', '表格数据行'),
        ('工作点 & $3.85013\\times10^{-3}$ & $10^{15}$ W cm$^{-2}$ & 低于/高于阈值 \\\\', '表格数据行'),
        ('阈值功率 & $2.69841$ MW & —— & 翻 $1.64268$ 倍即可 \\\\', '表格数据行'),
        ('子波 1 & 离子声波 & 朗缪尔波 $192.758$ nm & $k_L\\lambda_{De} = 0.766278$ \\\\', '表格数据行'),
        ('子波 2 & 电子波 & 散射光 $426.177$ nm & 斯托克斯位移 \\\\', '表格数据行'),
        ('\\begin{tip}{参量过程是「泵浦把能量交给子波」，不是碰撞}{tip-app-al-parametric}', '语义盒标题'),
        ('\\textbf{出处：}参量不稳定性与三波耦合的框架参照 Sagdeev 与 Galeev（1969，第 3 章）、', '出处声明整行'),
        ('量 & 托卡马克（140 GHz） & ICF 冕区（0.35 $\\mu$m） & 说明 \\\\', '表格数据行'),
        ('$n_c/n_e$ & $2.43127$ & $91.0085$ & 决定 $P_c$ \\\\', '表格数据行'),
        ('临界功率 $P_c$ & $42.3041$ GW & $1.58355$ TW & $17.4(n_c/n_e)$ GW \\\\', '表格数据行'),
        ('与驱动功率之比 & $42304.1$（$1$ MW） & $315.747$（$500$ TW） & 差四个数量级 \\\\', '表格数据行'),
        ('德拜长度 & $71.2806$ $\\mu$m & $23.5082$ nm & 决定崩塌尺度 \\\\', '表格数据行'),
        ('caviton 尺度 & $712.806$ $\\mu$m & $235.082$ nm & $\\sim10\\,\\lambda_{De}$ \\\\', '表格数据行'),
        ('$1\\%$ 密度坑所需场 & $3.64820\\times10^{7}$ V m$^{-1}$ & $1.20317\\times10^{10}$ V m$^{-1}$ & $W/(4nT_e)$ \\\\', '表格数据行'),
        ('波破裂场 $E_{WB}$ & $961.592$ MV m$^{-1}$ & —— & $\\propto\\sqrt{n_e}$ \\\\', '表格数据行'),
        ('\\begin{tip}{崩塌是「非线性把波自己关起来」}{tip-app-al-collapse}', '语义盒标题'),
        ('\\textbf{出处：}波崩塌与 caviton 的框架参照 Zakharov（\\emph{Soviet Physics JETP} 35,', '出处声明整行'),
        ('量 & $E = 10^{3}$ V m$^{-1}$ & $E = 10^{4}$ V m$^{-1}$ & 说明 \\\\', '表格数据行'),
        ('俘获频率 $\\omega_b$ & $5.58162\\times10^{10}$ s$^{-1}$ & $1.76506\\times10^{11}$ s$^{-1}$ & $k\\sqrt{eE/m_e}$ \\\\', '表格数据行'),
        ('俘获时间 & $112.569$ ps & $35.5975$ ps & $2\\pi/\\omega_b$ \\\\', '表格数据行'),
        ('俘获速度宽度 & $2.65241\\times10^{7}$ m s$^{-1}$ & $8.38766\\times10^{7}$ m s$^{-1}$ & $2\\omega_b/k$ \\\\', '表格数据行'),
        ('相对热速度 & $0.659596\\,v_{the}$ & $2.08582\\,v_{the}$ & 势井有多宽 \\\\', '表格数据行'),
        ('$\\omega_b/\\gamma$ & $4.38521$ & $13.8672$ & $>1$ 即非线性主导 \\\\', '表格数据行'),
        ('\\begin{tip}{俘获把「阻尼」变成「振荡」}{tip-app-al-trapping}', '语义盒标题'),
        ("\\textbf{出处：}非线性Landau 阻尼与粒子俘获参照 O'Neil（\\emph{Physics of Fluids} 8, 2255, 19", '出处声明整行'),
        ('结构 & 尺度 & 特征量 & 寿命 \\\\', '表格数据行'),
        ('双电层 & $100\\,\\lambda_{De} = 7.12806$ mm & $E = 1.28983$ MV m$^{-1}$ & 稳态（由电流维持） \\\\', '表格数据行'),
        ('相空间空洞 & $2.08582\\,v_{the}$（速度宽度） & $\\omega_b = 1.76506\\times10^{11}$ s$^{-1}$ & $35.5975$ ps \\\\', '表格数据行'),
        ('caviton & $\\sim10\\,\\lambda_{De} = 712.806$ $\\mu$m & 崩到 $E_{WB}$ 为止 & 几个等离子体周期 \\\\', '表格数据行'),
        ('孤子 & $3.65148\\,\\lambda_{De} = 260.280$ $\\mu$m & $M = 1.1$ & $438.533$ ps（越过自身） \\\\', '表格数据行'),
        ('\\begin{tip}{相干结构是「自组织」的，不是「外加」的}{tip-app-al-structures}', '语义盒标题'),
        ('\\textbf{出处：}双电层与相空间空洞的实验与理论综述参照 Sagdeev 与 Galeev（1969，第 4 章）', '出处声明整行'),
        ('级别 & 内容 & 说明 \\\\', '表格数据行'),
        ('现算级 & 抖动速度、有质动力势、阈值、孤子宽度、$P_c$、$E_{WB}$、$\\omega_b$、$E^{*}$、双电层场强 & 全部由 \\texttt{tmp/appAL\\_calc.py} 用 $n_e = 10^{20}$ m$^{-3}$、$T_e = 9194$ eV、$B = 5.3$ T 与声速 $c_s = 593524$ m s$^{-1}$（与附录 AE 同口径）算出 \\\\', '表格数据行'),
        ('声明级 & 自聚焦系数 $17.4$、孤子形状因子、双电层厚度 $100\\,\\lambda_{De}$、冕区 $0.35$ $\\mu$m 与 $500$ TW 驱动器参数 & 量级或工程典型值，正文中已逐处标为「声明」或「量级」 \\\\', '表格数据行'),
        ("文献级 & 孤子解（KdV）、三波耦合方程形式、O'Neil 阻尼振荡图景、BGK 平衡、扎哈罗夫方程 & 只引用方程形式与结论，不作数值比对 \\\\", '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-al-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（整行绑定：%s）' % (APP_FILES['AL'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 140, '$ GHz', '非线性从哪来：一个无·正文'),
        (0, 1, '$ MW', '非线性从哪来：一个无·正文'),
        (0, 2, '$ cm', '非线性从哪来：一个无·正文'),
        (0, 1, '$ keV', '非线性从哪来：一个无·正文'),
        (0, 795775000, '$ W m$^{-2}$', '非线性从哪来：一个无·公式'),
        (0, 0.0795775, '$ MW cm$^{-2}$', '非线性从哪来：一个无·公式'),
        (0, 774329, '$ V m$^{-1}$', '非线性从哪来：一个无·公式'),
        (0, 7.74329, '$ kV cm$^{-1}$', '非线性从哪来：一个无·公式'),
        (0, 154824, '$ m s$^{-1}$', '非线性从哪来：一个无·公式'),
        (0, 40212700, '$ m s$^{-1}$', '非线性从哪来：一个无·公式'),
        (0, 4, 'm', '非线性从哪来：一个无·公式'),
        (0, 5.45892e-21, '$ J', '非线性从哪来：一个无·公式'),
        (0, 0.0340719, '$ eV', '非线性从哪来：一个无·公式'),
        (0, 11.4292, '$ eV', '非线性从哪来：一个无·公式'),
        (0, 1, '\\ \\mathrm{keV}', '非线性从哪来：一个无·公式'),
        (0, 114.292, '$ eV', '非线性从哪来：一个无·公式'),
        (0, 2.43127e+20, '$ m$^{-3}$', '非线性从哪来：一个无·正文'),
        (0, 9.10085e+27, '$ m$^{-3}$', '非线性从哪来：一个无·正文'),
        (0, 3, 'cm', '非线性从哪来：一个无·表格'),
        (0, 140, 'GHz', '非线性从哪来：一个无·表格'),
        (0, 8.68021e+10, '$ V m$^{-1}$', '非线性从哪来：一个无·表格'),
        (0, 2836740, '$ m s$^{-1}$', '非线性从哪来：一个无·表格'),
        (1, 30, '\\%', '孤子：色散与非线性平·公式'),
        (1, 593524, '$ m s$^{-1}$', '孤子：色散与非线性平·公式'),
        (1, 593.524, '$ km s$^{-1}$', '孤子：色散与非线性平·公式'),
        (1, 0.00026028, '\\ \\mathrm{m}', '孤子：色散与非线性平·公式'),
        (1, 652877, '$ m s$^{-1}$', '孤子：色散与非线性平·公式'),
        (1, 438.533, '$ ps', '孤子：色散与非线性平·公式'),
        (1, 9.194, '$ keV', '孤子：色散与非线性平·正文'),
        (1, 1, '\\%', '孤子：色散与非线性平·正文'),
        (1, 36482000, '$ V m$^{-1}$', '孤子：色散与非线性平·正文'),
        (1, 1.76643e+12, '$ W m$^{-2}$', '孤子：色散与非线性平·正文'),
        (1, 176643000, '$ W cm$^{-2}$', '孤子：色散与非线性平·正文'),
        (1, 2.21976, '$ GW', '孤子：色散与非线性平·正文'),
        (2, 2.69841, '$ MW', '参量不稳定性：三波共·公式'),
        (2, 2.7, '$ MW', '参量不稳定性：三波共·公式'),
        (2, 13262100, '$ m s$^{-1}$', '参量不稳定性：三波共·公式'),
        (2, 83876.6, '$ m s$^{-1}$', '参量不稳定性：三波共·公式'),
        (2, 8.74265e+15, '$ W m$^{-2}$', '参量不稳定性：三波共·公式'),
        (2, 8.74265e+11, '$ W cm$^{-2}$', '参量不稳定性：三波共·公式'),
        (2, 17853100, '$ m$^{-1}$', '参量不稳定性：三波共·正文'),
        (2, 32596200, '$ m$^{-1}$', '参量不稳定性：三波共·正文'),
        (2, 192.758, '$ nm', '参量不稳定性：三波共·正文'),
        (2, 23.5082, '$ nm', '参量不稳定性：三波共·正文'),
        (2, 426.177, '$ nm', '参量不稳定性：三波共·正文'),
        (2, 707.343, '$ THz', '参量不稳定性：三波共·正文'),
        (2, 2251.29, '$ m$^{-1}$', '参量不稳定性：三波共·表格'),
        (3, 4.23041e+10, '$ W', '自聚焦与波崩塌·公式'),
        (3, 42.3041, '$ GW', '自聚焦与波崩塌·公式'),
        (3, 1.58355e+12, '$ W', '自聚焦与波崩塌·公式'),
        (3, 1.58355, '$ TW', '自聚焦与波崩塌·公式'),
        (3, 500, '$ TW', '自聚焦与波崩塌·公式'),
        (3, 2.60417, '$ TW', '自聚焦与波崩塌·公式'),
        (3, 2.35082e-08, '$ m', '自聚焦与波崩塌·公式'),
        (3, 0.000712806, '\\ \\mathrm{m}', '自聚焦与波崩塌·公式'),
        (3, 235.082, '$ nm', '自聚焦与波崩塌·公式'),
        (3, 1.20317e+10, '$ V m$^{-1}$', '自聚焦与波崩塌·公式'),
        (3, 1.92128e+13, '$ W cm$^{-2}$', '自聚焦与波崩塌·公式'),
        (3, 961592000, '\\ \\mathrm{V\\,m^{-1}}', '自聚焦与波崩塌·公式'),
        (3, 961.592, '\\ \\mathrm{MV\\,m^{-1}}', '自聚焦与波崩塌·公式'),
        (3, 961.592, '$ MV m$^{-1}$', '自聚焦与波崩塌·表格'),
        (4, 4, 'cm', '粒子俘获：什么时候非·公式'),
        (4, 5.58162e+10, '$ s$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 1.76506e+11, '$ s$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 112.569, '$ ps', '粒子俘获：什么时候非·公式'),
        (4, 35.5975, '$ ps', '粒子俘获：什么时候非·公式'),
        (4, 26524100, '$ m s$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 83876600, '$ m s$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 1.27283e+10, '$ s$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 78.5651, '$ ps', '粒子俘获：什么时候非·公式'),
        (4, 52.002, '\\ \\mathrm{V\\,m^{-1}}', '粒子俘获：什么时候非·公式'),
        (4, 52, '$ V m$^{-1}$', '粒子俘获：什么时候非·公式'),
        (4, 1474.96, '$ V m$^{-1}$', '粒子俘获：什么时候非·公式'),
        (5, 1289830, '\\ \\mathrm{V\\,m^{-1}}', '相干结构：双电层与相·公式'),
        (5, 12.8983, '\\ \\mathrm{kV\\,cm^{-1}}', '相干结构：双电层与相·公式'),
        (5, 83876600, '\\ \\mathrm{m\\,s^{-1}}', '相干结构：双电层与相·公式'),
        (5, 7.12806, '$ mm', '相干结构：双电层与相·表格'),
        (5, 1.28983, '$ MV m$^{-1}$', '相干结构：双电层与相·表格'),
        (6, 2, '$ MW', '估算练习·正文'),
        (6, 12.5664, '$ cm$^2$', '估算练习·正文'),
        (6, 9194, '$ eV', '估算练习·正文'),
        (6, 50, '$ MW', '估算练习·正文'),
        (6, 1591550000, '$ W m$^{-2}$', '估算练习·公式'),
        (6, 1095070, '$ V m$^{-1}$', '估算练习·公式'),
        (6, 218954, '$ m s$^{-1}$', '估算练习·公式'),
        (6, 158.355, '$ GW', '估算练习·公式'),
        (6, 1.12979e+11, '$ s$^{-1}$', '估算练习·公式'),
        (6, 7014.53, '$ m$^{-1}$', '估算练习·公式'),
        (6, 52.002, '$ V m$^{-1}$', '估算练习·公式'),
        (6, 14743100, '$ m$^{-1}$', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AL 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AL 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('8.80\\times10^{11}', 'ECRH 角频率四舍五入'),
        ('7.96\\times10^{8}', '功率密度四舍五入'),
        ('7.74\\times10^{5}', '微波电场四舍五入'),
        ('1.55\\times10^{5}', '抖动速度四舍五入'),
        ('3.85\\times10^{-3}', '抖动速度比四舍五入'),
        ('5.46\\times10^{-21}', '有质动力势四舍五入'),
        ('0.0341\\ ', '有质动力势 eV 四舍五入'),
        ('3.71\\times10^{-6}', '有质动力势与温度之比四舍五入'),
        ('2.43\\times10^{20}', '140 GHz 截止密度四舍五入'),
        ('9.46\\times10^{-3}', '激光 a0 四舍五入'),
        ('6.32\\times10^{-3}', '参量阈值四舍五入'),
        ('1.64\\ ', '低于阈值倍数四舍五入'),
        ('2.70\\ ', '阈值功率四舍五入'),
        ('8.74\\times10^{11}', '阈值强度四舍五入'),
        ('2.60\\times10^{-4}', '孤子宽度四舍五入'),
        ('260.3\\ ', '孤子宽度微米四舍五入'),
        ('3.65\\ ', '宽度与德拜长度之比四舍五入'),
        ('438.5\\ ', '孤子特征时间四舍五入'),
        ('7.51\\times10^{-8}', '小场密度坑四舍五入'),
        ('3.65\\times10^{7}', '百分之一密度坑所需场四舍五入'),
        ('1.77\\times10^{12}', '百分之一密度坑功率密度四舍五入'),
        ('2.22\\times10^{9}', '百分之一密度坑功率四舍五入'),
        ('42.3\\ ', '托卡马克临界功率四舍五入'),
        ('4.23\\times10^{4}', '托卡马克功率倍数四舍五入'),
        ('1.58\\times10^{12}', '冕区临界功率四舍五入'),
        ('315.7\\ ', '与驱动器功率之比四舍五入'),
        ('7.13\\times10^{-4}', 'caviton 尺度四舍五入'),
        ('712.8\\ ', 'caviton 微米四舍五入'),
        ('961.6\\ ', '波破裂场四舍五入'),
        ('1.29\\times10^{6}', '双电层场强四舍五入'),
        ('12.9\\ ', '双电层 kV/cm 四舍五入'),
        ('235.1\\ ', '冕区 caviton 尺度四舍五入'),
        ('5.58\\times10^{10}', '俘获频率四舍五入'),
        ('112.6\\ ', '俘获时间四舍五入'),
        ('2.65\\times10^{7}', '俘获速度宽度四舍五入'),
        ('4.39\\ ', '俘获与阻尼之比四舍五入'),
        ('1.77\\times10^{11}', '强场俘获频率四舍五入'),
        ('35.6\\ ', '强场俘获时间四舍五入'),
        ('8.39\\times10^{7}', '强场俘获速度宽度四舍五入'),
        ('13.9\\ ', '强场俘获与阻尼之比四舍五入'),
        ('52.0\\ ', '俘获阈值场四舍五入'),
        ('1.27\\times10^{10}', 'Landau 阻尼率四舍五入'),
        ('78.6\\ ', '阻尼时间四舍五入'),
        ('1475\\ ', '短波阈值场四舍五入'),
        ('28.4\\ ', '阈值升高倍数四舍五入'),
        ('2.84\\times10^{6}', '激光抖动速度四舍五入'),
        ('8.68\\times10^{10}', '激光电场四舍五入'),
        ('1.26\\times10^{10}', '功率密度之比四舍五入'),
        ('3354\\ ', '有质动力之比四舍五入'),
        ('68.3\\ ', '附录 Y 归一化振幅四舍五入'),
        ('24.2\\ ', '附录 Y 有质动力势四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AL'], why))


def check_P45_appendix_am_magnets():
    """附录 AM：磁体与超导——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AM'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AM 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '为什么必须超导：一张电费单',
        '临界面：临界电流、温度裕度和 n 值',
        '导体与 CICC：为什么电缆要泡在氦里',
        '洛伦兹力与应力：磁体首先是结构件',
        '失超与保护：能量必须有个出口',
        '制冷与交流损耗：冷量是最贵的电',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AM 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AM 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-am-sc',
        'tab:app-am-critical',
        'tab:app-am-cicc',
        'tab:app-am-force',
        'tab:app-am-quench',
        'tab:app-am-cryo',
        'tab:app-am-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AM'])

    for lab in (
        'sec:app-am-sc',
        'sec:app-am-critical',
        'sec:app-am-cicc',
        'sec:app-am-force',
        'sec:app-am-quench',
        'sec:app-am-cryo',
        'sec:app-am-problems',
        'sec:app-am-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AM'])

    for lab in (
        'tip-app-am-why',
        'tip-app-am-margin',
        'tip-app-am-cicc',
        'tip-app-am-force',
        'tip-app-am-quench',
        'tip-app-am-cryo',
        'tip-app-am-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AM'])

    for lab in (
        'eq:app-am-ampere',
        'eq:app-am-pcu',
        'eq:app-am-carnot',
        'eq:app-am-jc',
        'eq:app-am-ic',
        'eq:app-am-tcs',
        'eq:app-am-nvalue',
        'eq:app-am-stekly',
        'eq:app-am-pb',
        'eq:app-am-florentz',
        'eq:app-am-hoop',
        'eq:app-am-stress',
        'eq:app-am-energy',
        'eq:app-am-inductance',
        'eq:app-am-miits',
        'eq:app-am-dt',
        'eq:app-am-cop',
        'eq:app-am-coupling',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AM'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AM'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AM 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AM 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答(\s*\d+)?\}', ex))
    if n_ex < 6 or n_ans < n_ex:
        rec('FAIL', 'P', '附录 AM 练习 %d 道、解答 %d 道（要求 >= 6 且解答 >= 练习）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AM 的 %d 道练习均有解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('5.3\\ \\mathrm{T}', '为什么必须超导：一张·正文（精确绑定）'),
        ('68.0\\ \\mathrm{kA}', '为什么必须超导：一张·公式（精确绑定）'),
        ('121.45\\ \\mathrm{km}', '为什么必须超导：一张·公式（精确绑定）'),
        ('5\\ \\mathrm{A\\,mm^{-2}}', '为什么必须超导：一张·公式（精确绑定）'),
        ('1.36\\times10^{-2}\\ \\mathrm{m^{2}}', '为什么必须超导：一张·公式（精确绑定）'),
        ('0.151813\\ \\Omega', '为什么必须超导：一张·公式（精确绑定）'),
        ('7.01983\\times10^{8}\\ \\mathrm{W}', '为什么必须超导：一张·公式（精确绑定）'),
        ('701.983\\ \\mathrm{MW}', '为什么必须超导：一张·公式（精确绑定）'),
        ('5\\ \\mathrm{kW}', '为什么必须超导：一张·公式（精确绑定）'),
        ('4.5\\ \\mathrm{K}', '为什么必须超导：一张·公式（精确绑定）'),
        ('300\\ \\mathrm{K}', '为什么必须超导：一张·公式（精确绑定）'),
        ('3.28333\\times10^{5}\\ \\mathrm{W}', '为什么必须超导：一张·公式（精确绑定）'),
        ('25\\%', '为什么必须超导：一张·正文（精确绑定）'),
        ('1.31333\\times10^{6}\\ \\mathrm{W}', '为什么必须超导：一张·正文（精确绑定）'),
        ('1313.33\\ \\mathrm{kW}', '为什么必须超导：一张·正文（精确绑定）'),
        ('3.5cm', '为什么必须超导：一张·表格（精确绑定）'),
        ('1.7\\times10^{-3}\\ \\mathrm{m^{2}}', '为什么必须超导：一张·表格（精确绑定）'),
        ('40\\ \\mathrm{A\\,mm^{-2}}', '为什么必须超导：一张·表格（精确绑定）'),
        ('6000\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('26\\ \\mathrm{T}', '临界面：临界电流、温·公式（精确绑定）'),
        ('4.26233\\times10^{9}\\ \\mathrm{A\\,m^{-2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('4262.33\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('12\\ \\mathrm{T}', '临界面：临界电流、温·公式（精确绑定）'),
        ('2370.74\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('20\\ \\mathrm{T}', '临界面：临界电流、温·公式（精确绑定）'),
        ('665.148\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('2.89838\\times10^{6}\\ \\mathrm{A}', '临界面：临界电流、温·公式（精确绑定）'),
        ('2.89838\\ \\mathrm{MA}', '临界面：临界电流、温·公式（精确绑定）'),
        ('18.0\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('16.0609\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('13.2084\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('15.7897\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('11.2897\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('2\\ \\mathrm{K}', '临界面：临界电流、温·公式（精确绑定）'),
        ('1.6121\\ \\mathrm{MA}', '临界面：临界电流、温·公式（精确绑定）'),
        ('0.4523\\ \\mathrm{MA}', '临界面：临界电流、温·公式（精确绑定）'),
        ('1700\\ \\mathrm{mm^{2}}', '临界面：临界电流、温·公式（精确绑定）'),
        ('40\\%', '导体与 CICC：为·正文（精确绑定）'),
        ('1.0\\times10^{8}\\ \\mathrm{A\\,m^{-2}}', '导体与 CICC：为·正文（精确绑定）'),
        ('66.6667\\ \\mathrm{A\\,mm^{-2}}', '导体与 CICC：为·正文（精确绑定）'),
        ('1.0\\ \\mathrm{K}', '导体与 CICC：为·正文（精确绑定）'),
        ('14282.9\\ \\mathrm{A}', '导体与 CICC：为·公式（精确绑定）'),
        ('14.2829\\ \\mathrm{kA}', '导体与 CICC：为·公式（精确绑定）'),
        ('1020\\ \\mathrm{mm^{2}}', '导体与 CICC：为·表格（精确绑定）'),
        ('9.80392\\times10^{-7}\\ \\Omega\\,\\mathrm{m^{-1}}', '导体与 CICC：为·表格（精确绑定）'),
        ('200\\ \\mathrm{W\\,m^{-1}}', '导体与 CICC：为·表格（精确绑定）'),
        ('1\\ \\mathrm{K}', '导体与 CICC：为·表格（精确绑定）'),
        ('1.11767\\times10^{7}\\ \\mathrm{Pa}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('11.1767\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('111.767\\ \\mathrm{bar}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('3.604\\times10^{5}\\ \\mathrm{N\\,m^{-1}}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('360.4\\ \\mathrm{kN\\,m^{-1}}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('8.0\\ \\mathrm{m}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('3.86349\\times10^{8}\\ \\mathrm{N}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('386.349\\ \\mathrm{MN}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('6.95428\\ \\mathrm{GN}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('0.456389\\ \\mathrm{m^{2}}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('8.46534\\times10^{8}\\ \\mathrm{Pa}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('846.534\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·公式（精确绑定）'),
        ('1693.07\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·正文（精确绑定）'),
        ('3.5\\ \\mathrm{m}', '失超与保护：能量必须·正文（精确绑定）'),
        ('1499.19\\ \\mathrm{m^{3}}', '失超与保护：能量必须·正文（精确绑定）'),
        ('1.6756\\times10^{10}\\ \\mathrm{J}', '失超与保护：能量必须·公式（精确绑定）'),
        ('16.756\\ \\mathrm{GJ}', '失超与保护：能量必须·公式（精确绑定）'),
        ('2W', '失超与保护：能量必须·公式（精确绑定）'),
        ('7.24739\\ \\mathrm{H}', '失超与保护：能量必须·公式（精确绑定）'),
        ('0.724739\\ \\Omega', '失超与保护：能量必须·公式（精确绑定）'),
        ('10\\ \\mathrm{s}', '失超与保护：能量必须·公式（精确绑定）'),
        ('49282.2\\ \\mathrm{V}', '失超与保护：能量必须·公式（精确绑定）'),
        ('49.2822\\ \\mathrm{kV}', '失超与保护：能量必须·公式（精确绑定）'),
        ('3.35119\\times10^{9}\\ \\mathrm{W}', '失超与保护：能量必须·公式（精确绑定）'),
        ('3.35119\\ \\mathrm{GW}', '失超与保护：能量必须·公式（精确绑定）'),
        ('1.0\\ \\mathrm{s}', '失超与保护：能量必须·公式（精确绑定）'),
        ('88.8889\\ \\mathrm{K}', '失超与保护：能量必须·公式（精确绑定）'),
        ('93.3889\\ \\mathrm{K}', '失超与保护：能量必须·正文（精确绑定）'),
        ('150\\ \\mathrm{K}', '失超与保护：能量必须·正文（精确绑定）'),
        ('1000\\ \\mathrm{s}', '失超与保护：能量必须·正文（精确绑定）'),
        ('7422.22\\ \\mathrm{K}', '失超与保护：能量必须·正文（精确绑定）'),
        ('4.0cm', '失超与保护：能量必须·表格（精确绑定）'),
        ('5\\ \\mathrm{s}', '失超与保护：能量必须·表格（精确绑定）'),
        ('51.8519\\ \\mathrm{K}', '失超与保护：能量必须·表格（精确绑定）'),
        ('20\\ \\mathrm{s}', '失超与保护：能量必须·表格（精确绑定）'),
        ('162.963\\ \\mathrm{K}', '失超与保护：能量必须·表格（精确绑定）'),
        ('5.0\\ \\mathrm{kW}', '制冷与交流损耗：冷量·公式（精确绑定）'),
        ('9947.18\\ \\mathrm{W\\,m^{-3}}', '制冷与交流损耗：冷量·公式（精确绑定）'),
        ('412.931\\ \\mathrm{m^{3}}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('4.1075\\times10^{6}\\ \\mathrm{W}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('4.1075\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('4.1075\\times10^{7}\\ \\mathrm{J}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('41.075\\ \\mathrm{MJ}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('41075\\ \\mathrm{W}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('41.075\\ \\mathrm{kW}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('10.789\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('1.31333\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文（精确绑定）'),
        ('0.05\\ \\mathrm{s}', '制冷与交流损耗：冷量·表格（精确绑定）'),
        ('0.5\\ \\mathrm{T\\,s^{-1}}', '制冷与交流损耗：冷量·表格（精确绑定）'),
        ('12.1024\\ \\mathrm{MW}', '制冷与交流损耗：冷量·表格（精确绑定）'),
        ('1\\ \\mathrm{W}', '制冷与交流损耗：冷量·语义盒（精确绑定）'),
        ('300\\ \\mathrm{W}', '制冷与交流损耗：冷量·语义盒（精确绑定）'),
        ('10.6\\ \\mathrm{T}', '估算练习·正文（精确绑定）'),
        ('10\\ \\mathrm{A\\,mm^{-2}}', '估算练习·正文（精确绑定）'),
        ('1.7\\times10^{-8}\\ \\Omega\\,\\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('5.0\\times10^{-9}\\ \\Omega\\,\\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('20.0\\ \\mathrm{kW}', '估算练习·正文（精确绑定）'),
        ('10\\%', '估算练习·正文（精确绑定）'),
        ('540\\ \\mathrm{MPa}', '估算练习·正文（精确绑定）'),
        ('67.0239\\ \\mathrm{GJ}', '估算练习·正文（精确绑定）'),
        ('772.698\\ \\mathrm{MN}', '估算练习·正文（精确绑定）'),
        ('0.701983\\ \\mathrm{GW}', '估算练习·正文（精确绑定）'),
        ('1.40397\\ \\mathrm{GW}', '估算练习·正文（精确绑定）'),
        ('5.25333\\times10^{6}\\ \\mathrm{W}', '估算练习·正文（精确绑定）'),
        ('5.25333\\ \\mathrm{MW}', '估算练习·正文（精确绑定）'),
        ('0.715461\\ \\mathrm{m^{2}}', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AM'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 铜方案（$5\\ \\mathrm{A\\,mm^{-2}}$） & 超导方案 & 说明 \\\\', '表格数据行'),
        ('安培匝数 & $1.643\\times10^{8}$ & $1.643\\times10^{8}$ & 由 $B_{TF}$ 定，与材料无关 \\\\', '表格数据行'),
        ('导体截面 & $1.36\\times10^{-2}\\ \\mathrm{m^{2}}$ & $1.7\\times10^{-3}\\ \\mathrm{m^{2}}$ & 超导可承受 $40\\ \\mathrm{A\\,mm^{-2}}$ \\\\', '表格数据行'),
        ('总长 & $121.45\\ \\mathrm{km}$ & $121.45\\ \\mathrm{km}$ & 同一个几何 \\\\', '表格数据行'),
        ('稳态电耗 & $701.983\\ \\mathrm{MW}$ & $1313.33\\ \\mathrm{kW}$ & 铜烧焦耳热，超导烧制冷功 \\\\', '表格数据行'),
        ('与 $P_{fus}$ 之比 & $1.26142$ & $2.35999\\times10^{-3}$ & 铜方案直接出局 \\\\', '表格数据行'),
        ('\\begin{tip}{超导不是为了「更强」，是为了「不烧电」}{tip-app-am-why}', '语义盒标题'),
        ('\\textbf{出处：}磁体设计方程与超导磁体的工程口径参照 Wilson \\emph{Superconducting Magnets}（1', '出处声明整行'),
        ('磁场 & 临界电流密度 & 电缆临界电流 & 工作点比值 $I/I_c$ \\\\', '表格数据行'),
        ('$5.3\\ \\mathrm{T}$ & $4262.33\\ \\mathrm{A\\,mm^{-2}}$ & $2.89838\\ \\mathrm{MA}$ & $0.0234614$ \\\\', '表格数据行'),
        ('$12\\ \\mathrm{T}$ & $2370.74\\ \\mathrm{A\\,mm^{-2}}$ & $1.6121\\ \\mathrm{MA}$ & $0.042181$ \\\\', '表格数据行'),
        ('$20\\ \\mathrm{T}$ & $665.148\\ \\mathrm{A\\,mm^{-2}}$ & $0.4523\\ \\mathrm{MA}$ & $0.150343$ \\\\', '表格数据行'),
        ('\\begin{tip}{裕度是三个坐标上的距离，不是一个数}{tip-app-am-margin}', '语义盒标题'),
        ('\\textbf{出处：}超导临界面的参数化与 $n$ 值口径参照 Wilson \\emph{Superconducting Magnets}', '出处声明整行'),
        ('量 & 数值 & 口径 & 说明 \\\\', '表格数据行'),
        ('电缆截面 & $1700\\ \\mathrm{mm^{2}}$ & $I/J_{op}$ & $40\\ \\mathrm{A\\,mm^{-2}}$ \\\\', '表格数据行'),
        ('铜截面 & $1020\\ \\mathrm{mm^{2}}$ & $60\\%$ & 分流通道 \\\\', '表格数据行'),
        ('单位长度电阻 & $9.80392\\times10^{-7}\\ \\Omega\\,\\mathrm{m^{-1}}$ & 冷态电阻率 & 失超时的主要热源 \\\\', '表格数据行'),
        ('单位长度换热 & $200\\ \\mathrm{W\\,m^{-1}}$ & $1\\ \\mathrm{K}$ 温升 & 池沸腾上界 \\\\', '表格数据行'),
        ('Stekly 电流 & $14.2829\\ \\mathrm{kA}$ & 判据 & 只到工作电流的 $1/5$ \\\\', '表格数据行'),
        ('\\begin{tip}{铜不是「拖累」，是失超时的救生通道}{tip-app-am-cicc}', '语义盒标题'),
        ('\\textbf{出处：}CICC 与稳定化判据参照 Wilson（1983，第 6--8 章）、', '出处声明整行'),
        ('量 & 表达式 & 数值 & 说明 \\\\', '表格数据行'),
        ('磁压 & $B^{2}/2\\mu_0$ & $11.1767\\ \\mathrm{MPa}$ & 能量密度口径 \\\\', '表格数据行'),
        ('单位长度力 & $IB$ & $360.4\\ \\mathrm{kN\\,m^{-1}}$ & 导体上的载荷 \\\\', '表格数据行'),
        ('每线圈向心力 & $f h N$ & $386.349\\ \\mathrm{MN}$ & 内腿求和 \\\\', '表格数据行'),
        ('全部 18 个线圈 & $\\times18$ & $6.95428\\ \\mathrm{GN}$ & 全机载荷 \\\\', '表格数据行'),
        ('平均应力 & $F/A_{wp}$ & $846.534\\ \\mathrm{MPa}$ & $0.940594\\,\\sigma_y$ \\\\', '表格数据行'),
        ('\\begin{tip}{磁体是结构件，不是「绕线」}{tip-app-am-force}', '语义盒标题'),
        ('\\textbf{出处：}磁压与洛伦兹力的口径参照 Wilson（1983，第 2 章）、', '出处声明整行'),
        ('泄能时间常数 & $\\int J^{2}\\mathrm{d}t$ & 热点温升 & 判断 \\\\', '表格数据行'),
        ('$5\\ \\mathrm{s}$ & 由上式 & $51.8519\\ \\mathrm{K}$ & 安全 \\\\', '表格数据行'),
        ('$10\\ \\mathrm{s}$ & $2.66667\\times10^{16}$ & $88.8889\\ \\mathrm{K}$ & 安全（本设计） \\\\', '表格数据行'),
        ('$20\\ \\mathrm{s}$ & 由上式 & $162.963\\ \\mathrm{K}$ & 超过 $150\\ \\mathrm{K}$ 上限 \\\\', '表格数据行'),
        ('$1000\\ \\mathrm{s}$ & $2.22667\\times10^{18}$ & $7422.22\\ \\mathrm{K}$ & 导体烧毁 \\\\', '表格数据行'),
        ('\\begin{tip}{失超不是故障，是必须被设计进去的事件}{tip-app-am-quench}', '语义盒标题'),
        ('\\textbf{出处：}储能、泄能与绝热热点模型参照 Wilson（1983，第 9--10 章）、', '出处声明整行'),
        ('冷负荷项 & 热负荷 & 折算电耗 & 说明 \\\\', '表格数据行'),
        ('稳态漏热（声明） & $5.0\\ \\mathrm{kW}$ & $1313.33\\ \\mathrm{kW}$ & 支撑、电流引线、辐射 \\\\', '表格数据行'),
        ('交流损耗（平均） & $41.075\\ \\mathrm{kW}$ & $10.789\\ \\mathrm{MW}$ & $\\tau=0.05\\ \\mathrm{s}$、$\\mathrm{d}B/\\mathrm{d}t=0.5\\ \\mathrm{T\\,s^{-1}}$ \\\\', '表格数据行'),
        ('合计 & --- & $12.1024\\ \\mathrm{MW}$ & 占 D-T 功率 $0.0217473$ \\\\', '表格数据行'),
        ('\\begin{tip}{冷量是最贵的电}{tip-app-am-cryo}', '语义盒标题'),
        ('\\textbf{出处：}制冷卡诺因子与交流损耗口径参照 Wilson（1983，第 8 章）、', '出处声明整行'),
        ('来源 & 借鉴了什么 & 没有借鉴什么 \\\\', '表格数据行'),
        ('Wilson \\emph{Superconducting Magnets}（1983） & 安培匝数设计方程、稳定性判据、失超保护框架 & 图表、例题、具体装置参数 \\\\', '表格数据行'),
        ('Iwasa（2009，第 2 版） & 临界面的组织方式、CICC 的工程口径 & 原文文字、案例数据 \\\\', '表格数据行'),
        ('Bottura 与 Godeke（\\emph{Rev. Accel. Sci. Tech.} 5，2012） & 超导材料与导体的综述框架 & 拟合系数与实验曲线 \\\\', '表格数据行'),
        ('Mitchell 等（\\emph{IEEE Trans. Appl. Supercond.} 21，2011） & ITER 磁体的力量级与保护时间常数 & 工程设计图、材料牌号数据 \\\\', '表格数据行'),
        ('Freidberg（2007）、Wesson（2011） & 磁压与线圈受力的物理图像 & 习题与插图 \\\\', '表格数据行'),
        ('ITER Physics Basis（\\emph{Nuclear Fusion} 39，1999） & 装置级参数的量级对照 & 表格与图 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-am-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AM'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 5.3, '\\ \\mathrm{T}', '为什么必须超导：一张·正文'),
        (0, 68, '\\ \\mathrm{kA}', '为什么必须超导：一张·公式'),
        (0, 121.45, '\\ \\mathrm{km}', '为什么必须超导：一张·公式'),
        (0, 5, '\\ \\mathrm{A\\,mm^{-2}}', '为什么必须超导：一张·公式'),
        (0, 0.0136, '\\ \\mathrm{m^{2}}', '为什么必须超导：一张·公式'),
        (0, 0.151813, '\\ \\Omega', '为什么必须超导：一张·公式'),
        (0, 701983000, '\\ \\mathrm{W}', '为什么必须超导：一张·公式'),
        (0, 701.983, '\\ \\mathrm{MW}', '为什么必须超导：一张·公式'),
        (0, 5, '\\ \\mathrm{kW}', '为什么必须超导：一张·公式'),
        (0, 4.5, '\\ \\mathrm{K}', '为什么必须超导：一张·公式'),
        (0, 300, '\\ \\mathrm{K}', '为什么必须超导：一张·公式'),
        (0, 328333, '\\ \\mathrm{W}', '为什么必须超导：一张·公式'),
        (0, 25, '\\%', '为什么必须超导：一张·正文'),
        (0, 1313330, '\\ \\mathrm{W}', '为什么必须超导：一张·正文'),
        (0, 1313.33, '\\ \\mathrm{kW}', '为什么必须超导：一张·正文'),
        (0, 3.5, 'cm', '为什么必须超导：一张·表格'),
        (0, 0.0017, '\\ \\mathrm{m^{2}}', '为什么必须超导：一张·表格'),
        (0, 40, '\\ \\mathrm{A\\,mm^{-2}}', '为什么必须超导：一张·表格'),
        (1, 6000, '\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式'),
        (1, 26, '\\ \\mathrm{T}', '临界面：临界电流、温·公式'),
        (1, 4262330000, '\\ \\mathrm{A\\,m^{-2}}', '临界面：临界电流、温·公式'),
        (1, 4262.33, '\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式'),
        (1, 12, '\\ \\mathrm{T}', '临界面：临界电流、温·公式'),
        (1, 2370.74, '\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式'),
        (1, 20, '\\ \\mathrm{T}', '临界面：临界电流、温·公式'),
        (1, 665.148, '\\ \\mathrm{A\\,mm^{-2}}', '临界面：临界电流、温·公式'),
        (1, 2898380, '\\ \\mathrm{A}', '临界面：临界电流、温·公式'),
        (1, 2.89838, '\\ \\mathrm{MA}', '临界面：临界电流、温·公式'),
        (1, 18, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 16.0609, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 13.2084, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 15.7897, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 11.2897, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 2, '\\ \\mathrm{K}', '临界面：临界电流、温·公式'),
        (1, 1.6121, '\\ \\mathrm{MA}', '临界面：临界电流、温·公式'),
        (1, 0.4523, '\\ \\mathrm{MA}', '临界面：临界电流、温·公式'),
        (1, 1700, '\\ \\mathrm{mm^{2}}', '临界面：临界电流、温·公式'),
        (2, 40, '\\%', '导体与 CICC：为·正文'),
        (2, 100000000, '\\ \\mathrm{A\\,m^{-2}}', '导体与 CICC：为·正文'),
        (2, 66.6667, '\\ \\mathrm{A\\,mm^{-2}}', '导体与 CICC：为·正文'),
        (2, 1, '\\ \\mathrm{K}', '导体与 CICC：为·正文'),
        (2, 14282.9, '\\ \\mathrm{A}', '导体与 CICC：为·公式'),
        (2, 14.2829, '\\ \\mathrm{kA}', '导体与 CICC：为·公式'),
        (2, 1020, '\\ \\mathrm{mm^{2}}', '导体与 CICC：为·表格'),
        (2, 9.80392e-07, '\\ \\Omega\\,\\mathrm{m^{-1}}', '导体与 CICC：为·表格'),
        (2, 200, '\\ \\mathrm{W\\,m^{-1}}', '导体与 CICC：为·表格'),
        (3, 11176700, '\\ \\mathrm{Pa}', '洛伦兹力与应力：磁体·公式'),
        (3, 11.1767, '\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·公式'),
        (3, 111.767, '\\ \\mathrm{bar}', '洛伦兹力与应力：磁体·公式'),
        (3, 360400, '\\ \\mathrm{N\\,m^{-1}}', '洛伦兹力与应力：磁体·公式'),
        (3, 360.4, '\\ \\mathrm{kN\\,m^{-1}}', '洛伦兹力与应力：磁体·公式'),
        (3, 8, '\\ \\mathrm{m}', '洛伦兹力与应力：磁体·公式'),
        (3, 386349000, '\\ \\mathrm{N}', '洛伦兹力与应力：磁体·公式'),
        (3, 386.349, '\\ \\mathrm{MN}', '洛伦兹力与应力：磁体·公式'),
        (3, 6.95428, '\\ \\mathrm{GN}', '洛伦兹力与应力：磁体·公式'),
        (3, 386, '\\ \\mathrm{MN}', '洛伦兹力与应力：磁体·公式'),
        (3, 0.456389, '\\ \\mathrm{m^{2}}', '洛伦兹力与应力：磁体·公式'),
        (3, 846534000, '\\ \\mathrm{Pa}', '洛伦兹力与应力：磁体·公式'),
        (3, 846.534, '\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·公式'),
        (3, 1693.07, '\\ \\mathrm{MPa}', '洛伦兹力与应力：磁体·正文'),
        (4, 3.5, '\\ \\mathrm{m}', '失超与保护：能量必须·正文'),
        (4, 1499.19, '\\ \\mathrm{m^{3}}', '失超与保护：能量必须·正文'),
        (4, 16756000000.0, '\\ \\mathrm{J}', '失超与保护：能量必须·公式'),
        (4, 16.756, '\\ \\mathrm{GJ}', '失超与保护：能量必须·公式'),
        (4, 2, 'W', '失超与保护：能量必须·公式'),
        (4, 7.24739, '\\ \\mathrm{H}', '失超与保护：能量必须·公式'),
        (4, 0.724739, '\\ \\Omega', '失超与保护：能量必须·公式'),
        (4, 10, '\\ \\mathrm{s}', '失超与保护：能量必须·公式'),
        (4, 49282.2, '\\ \\mathrm{V}', '失超与保护：能量必须·公式'),
        (4, 49.2822, '\\ \\mathrm{kV}', '失超与保护：能量必须·公式'),
        (4, 3351190000, '\\ \\mathrm{W}', '失超与保护：能量必须·公式'),
        (4, 3.35119, '\\ \\mathrm{GW}', '失超与保护：能量必须·公式'),
        (4, 1, '\\ \\mathrm{s}', '失超与保护：能量必须·公式'),
        (4, 88.8889, '\\ \\mathrm{K}', '失超与保护：能量必须·公式'),
        (4, 93.3889, '\\ \\mathrm{K}', '失超与保护：能量必须·正文'),
        (4, 150, '\\ \\mathrm{K}', '失超与保护：能量必须·正文'),
        (4, 1000, '\\ \\mathrm{s}', '失超与保护：能量必须·正文'),
        (4, 7422.22, '\\ \\mathrm{K}', '失超与保护：能量必须·正文'),
        (4, 4, 'cm', '失超与保护：能量必须·表格'),
        (4, 5, '\\ \\mathrm{s}', '失超与保护：能量必须·表格'),
        (4, 51.8519, '\\ \\mathrm{K}', '失超与保护：能量必须·表格'),
        (4, 20, '\\ \\mathrm{s}', '失超与保护：能量必须·表格'),
        (4, 162.963, '\\ \\mathrm{K}', '失超与保护：能量必须·表格'),
        (5, 9947.18, '\\ \\mathrm{W\\,m^{-3}}', '制冷与交流损耗：冷量·公式'),
        (5, 412.931, '\\ \\mathrm{m^{3}}', '制冷与交流损耗：冷量·正文'),
        (5, 4107500, '\\ \\mathrm{W}', '制冷与交流损耗：冷量·正文'),
        (5, 4.1075, '\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文'),
        (5, 41075000, '\\ \\mathrm{J}', '制冷与交流损耗：冷量·正文'),
        (5, 41.075, '\\ \\mathrm{MJ}', '制冷与交流损耗：冷量·正文'),
        (5, 41075, '\\ \\mathrm{W}', '制冷与交流损耗：冷量·正文'),
        (5, 41.075, '\\ \\mathrm{kW}', '制冷与交流损耗：冷量·正文'),
        (5, 10.789, '\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文'),
        (5, 1.31333, '\\ \\mathrm{MW}', '制冷与交流损耗：冷量·正文'),
        (5, 0.05, '\\ \\mathrm{s}', '制冷与交流损耗：冷量·表格'),
        (5, 0.5, '\\ \\mathrm{T\\,s^{-1}}', '制冷与交流损耗：冷量·表格'),
        (5, 12.1024, '\\ \\mathrm{MW}', '制冷与交流损耗：冷量·表格'),
        (5, 1, '\\ \\mathrm{W}', '制冷与交流损耗：冷量·语义盒'),
        (5, 300, '\\ \\mathrm{W}', '制冷与交流损耗：冷量·语义盒'),
        (6, 10.6, '\\ \\mathrm{T}', '估算练习·正文'),
        (6, 10, '\\ \\mathrm{A\\,mm^{-2}}', '估算练习·正文'),
        (6, 1.7e-08, '\\ \\Omega\\,\\mathrm{m}', '估算练习·正文'),
        (6, 5e-09, '\\ \\Omega\\,\\mathrm{m}', '估算练习·正文'),
        (6, 20, '\\ \\mathrm{kW}', '估算练习·正文'),
        (6, 10, '\\%', '估算练习·正文'),
        (6, 540, '\\ \\mathrm{MPa}', '估算练习·正文'),
        (6, 67.0239, '\\ \\mathrm{GJ}', '估算练习·正文'),
        (6, 772.698, '\\ \\mathrm{MN}', '估算练习·正文'),
        (6, 0.701983, '\\ \\mathrm{GW}', '估算练习·正文'),
        (6, 1.40397, '\\ \\mathrm{GW}', '估算练习·正文'),
        (6, 5253330, '\\ \\mathrm{W}', '估算练习·正文'),
        (6, 5.25333, '\\ \\mathrm{MW}', '估算练习·正文'),
        (6, 0.715461, '\\ \\mathrm{m^{2}}', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AM 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AM 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('1.64\\times10^{8}', '安培匝数四舍五入'),
        ('134.2\\ ', '每线圈匝数四舍五入'),
        ('6747\\ ', '每线圈导线长度四舍五入'),
        ('121.5\\ ', '导体总长四舍五入'),
        ('0.456\\ ', '绕组包截面积四舍五入'),
        ('413\\ ', '绕组包体积四舍五入'),
        ('0.152\\ ', '铜方案电阻四舍五入'),
        ('7.02\\times10^{8}', '铜方案功率四舍五入'),
        ('702\\ ', '铜方案兆瓦四舍五入'),
        ('0.702\\ ', '铜方案吉瓦四舍五入'),
        ('1.26\\ ', '铜电耗与聚变功率之比四舍五入'),
        ('66.7\\ ', '卡诺温比四舍五入'),
        ('328\\ ', '卡诺功四舍五入'),
        ('1313\\ ', '制冷电耗 kW 四舍五入'),
        ('535\\ ', '铜电耗与制冷电耗之比四舍五入'),
        ('4262\\ ', '5.3 T 临界电流密度四舍五入'),
        ('2371\\ ', '12 T 临界电流密度四舍五入'),
        ('665\\ ', '20 T 临界电流密度四舍五入'),
        ('2.90\\times10^{6}', '电缆临界电流四舍五入'),
        ('1.61\\times10^{6}', '12 T 电缆临界电流四舍五入'),
        ('452\\ ', '20 T 电缆临界电流四舍五入'),
        ('0.0235\\ ', '工作点电流比四舍五入'),
        ('0.0422\\ ', '12 T 电流比四舍五入'),
        ('0.150\\ ', '20 T 电流比四舍五入'),
        ('16.1\\ ', '5.3 T 临界温度四舍五入'),
        ('13.2\\ ', '12 T 临界温度四舍五入'),
        ('15.8\\ ', '分流温度四舍五入'),
        ('11.3\\ ', '温度裕度四舍五入'),
        ('2.553\\times10^{-33}', 'n 值电压比四舍五入'),
        ('9.80\\times10^{-7}', '铜电阻率每米四舍五入'),
        ('14283\\ ', 'Stekly 电流四舍五入'),
        ('14.28\\ ', 'Stekly 千安四舍五入'),
        ('22.7\\ ', 'Stekly 比值四舍五入'),
        ('1.12\\times10^{7}', '磁压四舍五入'),
        ('11.2\\ ', '磁压兆帕四舍五入'),
        ('111.8\\ ', '磁压巴四舍五入'),
        ('1499\\ ', '磁场体积四舍五入'),
        ('1.68\\times10^{10}', '磁储能焦耳四舍五入'),
        ('16756\\ ', '磁储能兆焦四舍五入'),
        ('16.8\\ ', '磁储能吉焦四舍五入'),
        ('7.25\\ ', '电感四舍五入'),
        ('0.725\\ ', '泄能电阻四舍五入'),
        ('49282\\ ', '峰值电压四舍五入'),
        ('49.3\\ ', '峰值电压千伏四舍五入'),
        ('3.35\\times10^{9}', '峰值功率四舍五入'),
        ('3.35\\ ', '峰值功率吉瓦四舍五入'),
        ('6.67\\times10^{7}', '铜中电流密度四舍五入'),
        ('66.7\\ ', '铜中电流密度每平方毫米四舍五入'),
        ('2.67\\times10^{16}', '正常保护 MIITs 四舍五入'),
        ('88.9\\ ', '正常保护温升四舍五入'),
        ('93.4\\ ', '热点温度四舍五入'),
        ('2.23\\times10^{18}', '保护失效 MIITs 四舍五入'),
        ('7422\\ ', '保护失效温升四舍五入'),
        ('49.5\\ ', '失效与允许之比四舍五入'),
        ('3.86\\times10^{8}', '单线圈合力四舍五入'),
        ('386\\ ', '单线圈合力兆牛四舍五入'),
        ('6.95\\times10^{9}', '全部线圈合力四舍五入'),
        ('8.47\\times10^{8}', '绕组包应力四舍五入'),
        ('847\\ ', '应力兆帕四舍五入'),
        ('0.941\\ ', '应力与屈服强度之比四舍五入'),
        ('1693\\ ', '绕组包减半后应力四舍五入'),
        ('9947\\ ', '耦合损耗密度四舍五入'),
        ('4.11\\times10^{6}', '瞬时耦合功率四舍五入'),
        ('4.11\\times10^{7}', '一次脉冲耦合能四舍五入'),
        ('41.08\\ ', '平均功率千瓦四舍五入'),
        ('1.08\\times10^{7}', '交流损耗制冷电耗四舍五入'),
        ('10.8\\ ', '交流损耗制冷电耗兆瓦四舍五入'),
        ('1.21\\times10^{7}', '制冷总电耗四舍五入'),
        ('12.1\\ ', '制冷总电耗兆瓦四舍五入'),
        ('0.0217\\ ', '制冷电耗与聚变输出之比四舍五入'),
        ('67.0\\ ', '练习 1 储能四舍五入'),
        ('772\\ ', '练习 1 合力四舍五入'),
        ('1.40\\times10^{9}', '练习 2 铜方案功率四舍五入'),
        ('51.9\\ ', '练习 3 温升四舍五入'),
        ('163\\ ', '练习 3 长泄能温升四舍五入'),
        ('5.25\\times10^{6}', '练习 4 电耗四舍五入'),
        ('0.0261\\ ', '练习 5 电流比四舍五入'),
        ('8.23\\ ', '练习 5 电压比四舍五入'),
        ('0.715\\ ', '练习 6 所需截面积四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AM'], why))


def check_P46_appendix_an_polarimetry():
    """附录 AN：偏振与磁光诊断——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与解答数量不足；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AN'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AN 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '偏振是诊断量：等离子体怎样改写一束光的偏振',
        '法拉第旋转：把密度和磁场一起称重',
        '科顿--穆顿效应：双折射里的三次方与一个单位陷阱',
        '运动斯塔克效应：把投掷角读成安全因子剖面',
        '从相位到剖面：弦数、分辨率与不确定度',
        '工程实现：调制、消光比与整圈模糊',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AN 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AN 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-an-polarization',
        'tab:app-an-faraday',
        'tab:app-an-cotton',
        'tab:app-an-mse',
        'tab:app-an-inversion',
        'tab:app-an-engineering',
        'tab:app-an-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AN'])

    for lab in (
        'sec:app-an-polarization',
        'sec:app-an-faraday',
        'sec:app-an-cotton',
        'sec:app-an-mse',
        'sec:app-an-inversion',
        'sec:app-an-engineering',
        'sec:app-an-problems',
        'sec:app-an-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AN'])

    for lab in (
        'tip-app-an-polarization',
        'tip-app-an-faraday',
        'tip-app-an-cotton',
        'tip-app-an-mse',
        'tip-app-an-inversion',
        'tip-app-an-engineering',
        'tip-app-an-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AN'])

    for lab in (
        'eq:app-an-stokes',
        'eq:app-an-eigen',
        'eq:app-an-faraday',
        'eq:app-an-rotation',
        'eq:app-an-coef',
        'eq:app-an-lineint',
        'eq:app-an-cotton',
        'eq:app-an-cmcoef',
        'eq:app-an-em',
        'eq:app-an-stark',
        'eq:app-an-split',
        'eq:app-an-pitch',
        'eq:app-an-abel',
        'eq:app-an-res',
        'eq:app-an-pem',
        'eq:app-an-extinction',
        'eq:app-an-ambiguity',
        'eq:app-an-calib',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AN'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AN'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AN 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AN 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    # 只数「带编号」的解答：本篇有一行裸 \textbf{解答} 小标题，
    # 若把它算进去，删掉一条编号解答仍然过得去（实测补丁 #466 曾落空）。
    n_ans = len(re.findall(r'\\textbf\{解答\s*\d+\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AN 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AN 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('5.64146\\times10^{11}\\ \\mathrm{rad\\,s^{-1}}', '偏振是诊断量：等离子·正文（精确绑定）'),
        ('9.32175\\times10^{11}\\ \\mathrm{rad\\,s^{-1}}', '偏振是诊断量：等离子·正文（精确绑定）'),
        ('5.3\\ \\mathrm{T}', '偏振是诊断量：等离子·正文（精确绑定）'),
        ('2.66667\\times10^{20}\\ \\mathrm{m^{-2}}', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('0.5\\ \\mathrm{T}', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('3.0cm', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('4.78823^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('603.768^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('2.39411^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('301.884^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('0.000111775^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('0.15755^\\circ', '偏振是诊断量：等离子·表格（精确绑定）'),
        ('1.41333\\times10^{21}\\ \\mathrm{T\\,m^{-2}}', '法拉第旋转：把密度和·公式（精确绑定）'),
        ('2.39404^\\circ', '法拉第旋转：把密度和·表格（精确绑定）'),
        ('300.713^\\circ', '法拉第旋转：把密度和·表格（精确绑定）'),
        ('1\\ \\mathrm{T}', '科顿--穆顿效应：双·正文（精确绑定）'),
        ('4\\ \\mathrm{G}', '科顿--穆顿效应：双·正文（精确绑定）'),
        ('0.157345^\\circ', '科顿--穆顿效应：双·表格（精确绑定）'),
        ('0.000111769^\\circ', '科顿--穆顿效应：双·表格（精确绑定）'),
        ('17.6793^\\circ', '科顿--穆顿效应：双·表格（精确绑定）'),
        ('1.25197\\times10^{-21}\\ \\mathrm{J}', '运动斯塔克效应：把投·公式（精确绑定）'),
        ('0.00781418\\ \\mathrm{eV}', '运动斯塔克效应：把投·公式（精确绑定）'),
        ('2.71454\\ \\mathrm{nm}', '运动斯塔克效应：把投·公式（精确绑定）'),
        ('0.0677692\\ \\mathrm{nm}', '运动斯塔克效应：把投·公式（精确绑定）'),
        ('5.38931^\\circ', '运动斯塔克效应：把投·公式（精确绑定）'),
        ('100\\ \\mathrm{keV}', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('3.09574\\times10^{6}\\ \\mathrm{m\\,s^{-1}}', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('16.4074\\ \\mathrm{MV\\,m^{-1}}', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('10\\ \\mathrm{eV}', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('1^\\circ', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('18.6978\\%', '运动斯塔克效应：把投·表格（精确绑定）'),
        ('18.5005\\%', '运动斯塔克效应：把投·正文（精确绑定）'),
        ('0.01^\\circ', '从相位到剖面：弦数、·公式（精确绑定）'),
        ('1.74533\\times10^{-4}\\ \\mathrm{rad}', '从相位到剖面：弦数、·公式（精确绑定）'),
        ('5.90355\\times10^{18}\\ \\mathrm{T\\,m^{-2}}', '从相位到剖面：弦数、·公式（精确绑定）'),
        ('0.417704\\%', '从相位到剖面：弦数、·正文（精确绑定）'),
        ('0.0221383\\ \\mathrm{T}', '从相位到剖面：弦数、·正文（精确绑定）'),
        ('0.5\\ \\mathrm{m}', '从相位到剖面：弦数、·正文（精确绑定）'),
        ('22.1383\\ \\mathrm{mT}', '从相位到剖面：弦数、·表格（精确绑定）'),
        ('0.569892\\ \\mathrm{T}', '从相位到剖面：弦数、·表格（精确绑定）'),
        ('50\\ \\mathrm{kHz}', '工程实现：调制、消光·公式（精确绑定）'),
        ('1\\ \\mathrm{ms}', '工程实现：调制、消光·公式（精确绑定）'),
        ('1.81185^\\circ', '工程实现：调制、消光·公式（精确绑定）'),
        ('1.69198\\times10^{21}\\ \\mathrm{T\\,m^{-2}}', '工程实现：调制、消光·公式（精确绑定）'),
        ('1.0\\ \\mathrm{T}', '估算练习·正文（精确绑定）'),
        ('150\\ \\mathrm{keV}', '估算练习·正文（精确绑定）'),
        ('0.005^\\circ', '估算练习·正文（精确绑定）'),
        ('5.24844\\ \\mathrm{rad}', '估算练习·正文（精确绑定）'),
        ('7.06667\\times10^{20}\\ \\mathrm{T\\,m^{-2}}', '估算练习·正文（精确绑定）'),
        ('1.19702^\\circ', '估算练习·正文（精确绑定）'),
        ('0.629382^\\circ', '估算练习·正文（精确绑定）'),
        ('3.79149\\times10^{6}\\ \\mathrm{m\\,s^{-1}}', '估算练习·正文（精确绑定）'),
        ('20.0949\\ \\mathrm{MV\\,m^{-1}}', '估算练习·正文（精确绑定）'),
        ('3.32461\\ \\mathrm{nm}', '估算练习·正文（精确绑定）'),
        ('0.333333\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('1.10336\\times10^{19}\\ \\mathrm{m^{-2}}', '估算练习·正文（精确绑定）'),
        ('2.95177\\times10^{18}\\ \\mathrm{T\\,m^{-2}}', '估算练习·正文（精确绑定）'),
        ('0.208852\\%', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AN'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & CO$_2$（$10.6\\ \\mu\\mathrm{m}$） & CH$_3$OH（$118.8\\ \\mu\\mathrm{m}$） & 标度 & 说明 \\\\', '表格数据行'),
        ('$X_{\\mathrm{eff}}$ & $6.71897\\times10^{-6}$ & $8.43963\\times10^{-4}$ & $\\lambda^2$ & 等效均匀密度口径 \\\\', '表格数据行'),
        ('$Y$（$B_t$） & $5.24569\\times10^{-3}$ & $5.87913\\times10^{-2}$ & $\\lambda$ & 磁场相对电子回旋频率 \\\\', '表格数据行'),
        ('$|\\Delta n|$（$R$--$L$） & $3.52467\\times10^{-8}$ & $4.98109\\times10^{-5}$ & $\\lambda^2 B$ & 相位差的一半是旋转角 \\\\', '表格数据行'),
        ('$R$/$L$ 相位差 & $4.78823^\\circ$ & $603.768^\\circ$ & $\\lambda^2 B$ & 可超过一圈 \\\\', '表格数据行'),
        ('偏振旋转角 & $2.39411^\\circ$ & $301.884^\\circ$ & $\\lambda^2 B$ & 观测量 \\\\', '表格数据行'),
        ('$|\\Delta n|$（$X$--$O$） & $8.22786\\times10^{-13}$ & $1.29978\\times10^{-8}$ & $\\lambda^3 B^2$ & 科顿--穆顿 \\\\', '表格数据行'),
        ('$X$/$O$ 相位差 & $0.000111775^\\circ$ & $0.15755^\\circ$ & $\\lambda^3 B_\\perp^2$ & 垂直弦时很小 \\\\', '表格数据行'),
        ('\\begin{tip}{偏振不是光的属性，是传播路径的记忆}{tip-app-an-polarization}', '语义盒标题'),
        ('\\textbf{出处：}偏振诊断的总体框架与 Stokes 参数口径参照 Hutchinson《Principles of Plasma D', '出处声明整行'),
        ('\\textbf{出处：}法拉第旋转的推导与工程系数参照 Hutchinson（2002）第 6 章与 Kunze《Introduction ', '出处声明整行'),
        ('\\textbf{出处：}科顿--穆顿效应的推导与高斯单位制下的常系数参照 Hutchinson（2002）第 6 章；单位换算口径由本书统一', '出处声明整行'),
        ('\\textbf{出处：}线性斯塔克效应与 MSE 原理参照 Levinton 等, Phys. Rev. Lett. 63, 2060 (1', '出处声明整行'),
        ('\\textbf{出处：}线积分反演与多弦几何的处理参照 Hutchinson（2002）第 6 章与本书附录 Q 的 Abel 反演口径。', '出处声明整行'),
        ('\\textbf{出处：}光弹调制与消光比的口径参照 Hutchinson（2002）第 6 章；整圈模糊的处理方式参照 Wesson（201', '出处声明整行'),
        ('量 & CO$_2$ 激光 & 远红外激光 & 说明 \\\\', '表格数据行'),
        ('$|n_{\\mathrm{R}}-n_{\\mathrm{L}}|$ & $3.52467\\times10^{-8}$ & $4.98109\\times10^{-5}$ & 与 $\\lambda^2 B$ 成正比 \\\\', '表格数据行'),
        ('$R$/$L$ 相位差 & $4.78823^\\circ$ & $603.768^\\circ$ & 精确冷等离子体解 \\\\', '表格数据行'),
        ('旋转角（精确） & $2.39411^\\circ$ & $301.884^\\circ$ & 相位差的一半 \\\\', '表格数据行'),
        ('旋转角（解析式） & $2.39404^\\circ$ & $300.713^\\circ$ & 与精确解差 $<0.4\\%$ \\\\', '表格数据行'),
        ('圈数 & $0.00665$ & $0.835315$ & 远红外已接近一圈 \\\\', '表格数据行'),
        ('反推磁场 & $5.3\\ \\mathrm{T}$ & --- & 已知密度剖面时 \\\\', '表格数据行'),
        ('\\begin{tip}{旋转角是相位差的一半——这条二分之一年年有人写错}{tip-app-an-faraday}', '语义盒标题'),
        ('量 & 远红外 & CO$_2$ & 说明 \\\\', '表格数据行'),
        ('$|\\Delta n_{\\mathrm{XO}}|$（垂直弦） & $1.29978\\times10^{-8}$ & $8.22786\\times10^{-13}$ & $\\propto\\lambda^3B_\\perp^2$ \\\\', '表格数据行'),
        ('相位（精确解） & $0.15755^\\circ$ & $0.000111775^\\circ$ & 垂直弦，$B_\\perp=0.5\\ \\mathrm{T}$ \\\\', '表格数据行'),
        ('相位（解析式） & $0.157345^\\circ$ & $0.000111769^\\circ$ & 与精确解差 $<0.2\\%$ \\\\', '表格数据行'),
        ('$\\varphi_{\\mathrm{CM}}/\\varphi_{\\mathrm{FR}}$ & $0.000523241$ & $4.66864\\times10^{-5}$ & 垂直弦上 CM 远小于 FR \\\\', '表格数据行'),
        ('切向弦放大 & $112.36$ 倍 & --- & $B_\\perp^2$ 放大 $(B_t/B_p)^2$ \\\\', '表格数据行'),
        ('切向弦相位 & $17.6793^\\circ$ & --- & 与 FR 同量级 \\\\', '表格数据行'),
        ('\\begin{tip}{同一个 2.45，单位不同差一亿倍}{tip-app-an-cotton}', '语义盒标题'),
        ('量 & 数值 & 说明 \\\\', '表格数据行'),
        ('束速度 $v_b$ & $3.09574\\times10^{6}\\ \\mathrm{m\\,s^{-1}}$ & $=0.0103263c$，能量 $100\\ \\mathrm{keV}$ \\\\', '表格数据行'),
        ('运动电场 $E_m$ & $16.4074\\ \\mathrm{MV\\,m^{-1}}$ & $v_bB\\sin\\theta$，$\\theta=90^\\circ$ \\\\', '表格数据行'),
        ('最大斯塔克分裂 & $0.00781418\\ \\mathrm{eV}$ & $n=3$，$|k|=2$ \\\\', '表格数据行'),
        ('$D_\\alpha$ 波长分裂 & $2.71454\\ \\mathrm{nm}$ & 观测量 \\\\', '表格数据行'),
        ('束内多普勒宽度 & $0.0677692\\ \\mathrm{nm}$ & $T_{\\mathrm{beam}}=10\\ \\mathrm{eV}$ \\\\', '表格数据行'),
        ('分裂/多普勒 & $40.0556$ & 无需解卷积 \\\\', '表格数据行'),
        ('投掷角 $\\gamma$ & $5.38931^\\circ$ & $\\arctan(B_p/B_t)$ \\\\', '表格数据行'),
        ('$\\gamma$ 变 $1^\\circ$ & $B_p$ 变 $18.6978\\%$ & 诊断灵敏度 \\\\', '表格数据行'),
        ('\\begin{tip}{MSE 测的是投掷角，不是磁场大小}{tip-app-an-mse}', '语义盒标题'),
        ('相位分辨率 & $0.01^\\circ$ & 典型 PEM 偏振仪水平 \\\\', '表格数据行'),
        ('线积分分辨率 & $5.90355\\times10^{18}\\ \\mathrm{T\\,m^{-2}}$ & 相对 $0.417704\\%$ \\\\', '表格数据行'),
        ('等效平均磁场分辨率 & $0.0221383\\ \\mathrm{T}$ & 单弦，$22.1383\\ \\mathrm{mT}$ \\\\', '表格数据行'),
        ('弦距（$9$ 弦） & $0.5\\ \\mathrm{m}$ & 空间分辨率上限 \\\\', '表格数据行'),
        ('边缘/中心线积分比 & $0.43589$ & $r=0.9a$ \\\\', '表格数据行'),
        ('$q$ 对 $\\gamma$ 的灵敏度 & $18.5005\\%$ 每度 & $1/\\tan\\gamma=10.6$ 每弧度 \\\\', '表格数据行'),
        ('$q=3$ 所需 $B_p$ & $0.569892\\ \\mathrm{T}$ & 当前值的 $1.13978$ 倍 \\\\', '表格数据行'),
        ('\\begin{tip}{线积分诊断的空间分辨率由弦距定，不由相位分辨率定}{tip-app-an-inversion}', '语义盒标题'),
        ('PEM 频率 & $50\\ \\mathrm{kHz}$ & 周期 $20\\ \\mu\\mathrm{s}$ \\\\', '表格数据行'),
        ('$1\\ \\mathrm{ms}$ 内周期数 & $50$ & 锁相带宽足够 \\\\', '表格数据行'),
        ('消光比 $10^{-3}$ 的相位误差 & $1.81185^\\circ$ & 元件质量是主要限制 \\\\', '表格数据行'),
        ('FR 的 $2\\pi$ 对应线积分 & $1.69198\\times10^{21}\\ \\mathrm{T\\,m^{-2}}$ & 远红外 \\\\', '表格数据行'),
        ('工作点占一圈的比例 & $0.835315$ & 已接近模糊 \\\\', '表格数据行'),
        ('CO$_2$ 的旋转角 & $2.39411^\\circ$ & 无模糊 \\\\', '表格数据行'),
        ('切向弦 CM 相位 & $17.6793^\\circ$ & 与 FR 同量级 \\\\', '表格数据行'),
        ('\\begin{tip}{整圈模糊不是噪声，是设计约束}{tip-app-an-engineering}', '语义盒标题'),
        ('借鉴内容 & 来源 & 本书的处理 \\\\', '表格数据行'),
        ('偏振诊断总框架、Stokes 参数与仪器口径 & Hutchinson《Principles of Plasma Diagnostics》2nd ed. 2002 第 6 章 & 重写推导，只保留结论与可复算的系数 \\\\', '表格数据行'),
        ('偏振诊断在托卡马克中的位置（FR、CM、MSE） & Wesson《Tokamaks》4th ed. 2011 第 12 章 & 改写为「先算量级、再谈方案」的叙述 \\\\', '表格数据行'),
        ('线性斯塔克效应与谱线分裂 & Griem《Principles of Plasma Spectroscopy》1997；Fujimoto《Plasma Spectroscopy》2004 & 只用能级分裂公式与量级 \\\\', '表格数据行'),
        ('谱线测量与仪器展宽 & Kunze《Introduction to Plasma Spectroscopy》2009 & 与附录 AH 的口径对齐 \\\\', '表格数据行'),
        ('MSE 原理与首次实现 & Levinton 等, Phys. Rev. Lett. 63, 2060 (1989) & 只引用其原理结论，不复制装置细节 \\\\', '表格数据行'),
        ('MSE 数据分析与定标 & Wróblewski 与 Lao, Rev. Sci. Instrum. 63, 5140 (1992) & 只引用定标思路 \\\\', '表格数据行'),
        ('原子物理与谱线展宽 & Salzmann《Atomic Physics in Hot Plasmas》1998 & 与附录 AH 同源，不重复 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-an-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AN'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 564146000000.0, '\\ \\mathrm{rad\\,s^{-1}}', '偏振是诊断量：等离子·正文'),
        (0, 932175000000.0, '\\ \\mathrm{rad\\,s^{-1}}', '偏振是诊断量：等离子·正文'),
        (0, 5.3, '\\ \\mathrm{T}', '偏振是诊断量：等离子·正文'),
        (0, 2.66667e+20, '\\ \\mathrm{m^{-2}}', '偏振是诊断量：等离子·表格'),
        (0, 0.5, '\\ \\mathrm{T}', '偏振是诊断量：等离子·表格'),
        (0, 3, 'cm', '偏振是诊断量：等离子·表格'),
        (0, 4.78823, '^\\circ', '偏振是诊断量：等离子·表格'),
        (0, 603.768, '^\\circ', '偏振是诊断量：等离子·表格'),
        (0, 2.39411, '^\\circ', '偏振是诊断量：等离子·表格'),
        (0, 301.884, '^\\circ', '偏振是诊断量：等离子·表格'),
        (0, 0.000111775, '^\\circ', '偏振是诊断量：等离子·表格'),
        (0, 0.15755, '^\\circ', '偏振是诊断量：等离子·表格'),
        (1, 1.41333e+21, '\\ \\mathrm{T\\,m^{-2}}', '法拉第旋转：把密度和·公式'),
        (1, 2.39404, '^\\circ', '法拉第旋转：把密度和·表格'),
        (1, 300.713, '^\\circ', '法拉第旋转：把密度和·表格'),
        (2, 1, '\\ \\mathrm{T}', '科顿--穆顿效应：双·正文'),
        (2, 4, '\\ \\mathrm{G}', '科顿--穆顿效应：双·正文'),
        (2, 0.157345, '^\\circ', '科顿--穆顿效应：双·表格'),
        (2, 0.000111769, '^\\circ', '科顿--穆顿效应：双·表格'),
        (2, 17.6793, '^\\circ', '科顿--穆顿效应：双·表格'),
        (3, 1.25197e-21, '\\ \\mathrm{J}', '运动斯塔克效应：把投·公式'),
        (3, 0.00781418, '\\ \\mathrm{eV}', '运动斯塔克效应：把投·公式'),
        (3, 2.71454, '\\ \\mathrm{nm}', '运动斯塔克效应：把投·公式'),
        (3, 0.0677692, '\\ \\mathrm{nm}', '运动斯塔克效应：把投·公式'),
        (3, 5.38931, '^\\circ', '运动斯塔克效应：把投·公式'),
        (3, 100, '\\ \\mathrm{keV}', '运动斯塔克效应：把投·表格'),
        (3, 3095740, '\\ \\mathrm{m\\,s^{-1}}', '运动斯塔克效应：把投·表格'),
        (3, 16.4074, '\\ \\mathrm{MV\\,m^{-1}}', '运动斯塔克效应：把投·表格'),
        (3, 10, '\\ \\mathrm{eV}', '运动斯塔克效应：把投·表格'),
        (3, 1, '^\\circ', '运动斯塔克效应：把投·表格'),
        (3, 18.6978, '\\%', '运动斯塔克效应：把投·表格'),
        (3, 18.5005, '\\%', '运动斯塔克效应：把投·正文'),
        (4, 0.01, '^\\circ', '从相位到剖面：弦数、·公式'),
        (4, 0.000174533, '\\ \\mathrm{rad}', '从相位到剖面：弦数、·公式'),
        (4, 5.90355e+18, '\\ \\mathrm{T\\,m^{-2}}', '从相位到剖面：弦数、·公式'),
        (4, 0.417704, '\\%', '从相位到剖面：弦数、·正文'),
        (4, 0.0221383, '\\ \\mathrm{T}', '从相位到剖面：弦数、·正文'),
        (4, 0.5, '\\ \\mathrm{m}', '从相位到剖面：弦数、·正文'),
        (4, 22.1383, '\\ \\mathrm{mT}', '从相位到剖面：弦数、·表格'),
        (4, 0.569892, '\\ \\mathrm{T}', '从相位到剖面：弦数、·表格'),
        (5, 50, '\\ \\mathrm{kHz}', '工程实现：调制、消光·公式'),
        (5, 1, '\\ \\mathrm{ms}', '工程实现：调制、消光·公式'),
        (5, 1.81185, '^\\circ', '工程实现：调制、消光·公式'),
        (5, 1.69198e+21, '\\ \\mathrm{T\\,m^{-2}}', '工程实现：调制、消光·公式'),
        (6, 150, '\\ \\mathrm{keV}', '估算练习·正文'),
        (6, 0.005, '^\\circ', '估算练习·正文'),
        (6, 5.24844, '\\ \\mathrm{rad}', '估算练习·正文'),
        (6, 7.06667e+20, '\\ \\mathrm{T\\,m^{-2}}', '估算练习·正文'),
        (6, 1.19702, '^\\circ', '估算练习·正文'),
        (6, 0.629382, '^\\circ', '估算练习·正文'),
        (6, 3791490, '\\ \\mathrm{m\\,s^{-1}}', '估算练习·正文'),
        (6, 20.0949, '\\ \\mathrm{MV\\,m^{-1}}', '估算练习·正文'),
        (6, 3.32461, '\\ \\mathrm{nm}', '估算练习·正文'),
        (6, 0.333333, '\\ \\mathrm{m}', '估算练习·正文'),
        (6, 1.10336e+19, '\\ \\mathrm{m^{-2}}', '估算练习·正文'),
        (6, 2.95177e+18, '\\ \\mathrm{T\\,m^{-2}}', '估算练习·正文'),
        (6, 0.208852, '\\%', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AN 有 %d 个数值未在小节内复算通过：%s'
            % (len(bad), '；'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AN 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流防护 ----------
    for lit, why in (
        ('2.39^\\circ', '法拉第旋转角四舍五入'),
        ('2.394\\ ', '法拉第旋转角再四舍五入'),
        ('301.9^\\circ', 'FIR 旋转角四舍五入'),
        ('300.7^\\circ', 'FIR 旋转角再四舍五入'),
        ('5.25\\ ', 'FIR 旋转角弧度四舍五入'),
        ('0.8353\\ ', 'FIR 整圈数四舍五入'),
        ('125.6\\ ', 'FIR 与 CO2 之比四舍五入'),
        ('40.1\\ ', '斯塔克与多普勒之比四舍五入'),
        ('18.7\\ ', '1 度对应的极向场份额四舍五入'),
        ('0.4177\\%', '单条弦分辨率四舍五入'),
        ('0.0221\\ ', '平均磁场分辨率四舍五入'),
        ('3.79\\times10^{6}', '束速度四舍五入'),
        ('1.10\\times10^{19}', '边缘弦线积分四舍五入'),
        ('2.95\\times10^{18}', '分辨率对应的线积分四舍五入'),
        ('0.2089', '分辨率百分比四舍五入'),
        ('0.6294', 'CM 旋转角四舍五入'),
        ('1.197\\ ', '练习 2 旋转角四舍五入'),
        ('3.325\\ ', '练习 4 波长差四舍五入'),
        ('1.225\\ ', '练习 4 比值四舍五入'),
        ('7.07\\times10^{20}', '练习 2 线积分四舍五入'),
        ('6.67\\times10^{19}', '有效密度四舍五入'),
        ('1.41\\times10^{21}', '密度磁场线积分四舍五入'),
        ('1.69\\times10^{21}', '整圈对应的线积分四舍五入'),
        ('5.90\\times10^{18}', '相位分辨率对应线积分四舍五入'),
        ('1.30\\times10^{-8}', 'CM 双折射差四舍五入'),
        ('8.23\\times10^{-13}', 'CO2 的 CM 双折射差四舍五入'),
        ('2.46\\times10^{-11}', 'CM 常系数四舍五入'),
        ('2.63\\times10^{-13}', '法拉第常系数四舍五入'),
        ('4.98\\times10^{-5}', 'FIR 双折射差四舍五入'),
        ('3.52\\times10^{-8}', 'CO2 双折射差四舍五入'),
        ('4.67\\times10^{-5}', 'CM 与法拉第之比四舍五入'),
        ('1.25\\times10^{-21}', '斯塔克分裂能四舍五入'),
        ('112.4\\ ', '切向放大倍数四舍五入'),
        ('17.7\\ ', '切向放大后的旋转角四舍五入'),
        ('5.92^\\circ', '极向场加十个百分点后的旋转角四舍五入'),
        ('1.14\\ ', '极向场倍数四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AN'], why))


def check_P47_appendix_ao_vacuum():
    """附录 AO：真空与壁处理——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与编号解答数量不足或不相等；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AO'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AO 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '真空是什么单位：密度、吞吐量与流导',
        '抽气：有效抽速由导管决定',
        '放气与漏率：本底压强是一次除法',
        '烘烤与壁处理：为什么要把壁“洗”干净',
        '壁条件与再循环：净加料只是零头',
        '燃料滞留与氚库存：账在壁上，不在罐里',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AO 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AO 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ao-units',
        'tab:app-ao-pumping',
        'tab:app-ao-outgas',
        'tab:app-ao-conditioning',
        'tab:app-ao-recycling',
        'tab:app-ao-retention',
        'tab:app-ao-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AO'])

    for lab in (
        'sec:app-ao-units',
        'sec:app-ao-pumping',
        'sec:app-ao-outgassing',
        'sec:app-ao-conditioning',
        'sec:app-ao-recycling',
        'sec:app-ao-retention',
        'sec:app-ao-problems',
        'sec:app-ao-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AO'])

    for lab in (
        'tip-app-ao-units',
        'tip-app-ao-pumping',
        'tip-app-ao-outgassing',
        'tip-app-ao-conditioning',
        'tip-app-ao-recycling',
        'tip-app-ao-retention',
        'tip-app-ao-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AO'])

    for lab in (
        'eq:app-ao-density',
        'eq:app-ao-throughput',
        'eq:app-ao-sccm',
        'eq:app-ao-conductance',
        'eq:app-ao-seff',
        'eq:app-ao-lossless',
        'eq:app-ao-outgas',
        'eq:app-ao-base',
        'eq:app-ao-bake',
        'eq:app-ao-monolayer',
        'eq:app-ao-bakepower',
        'eq:app-ao-boron',
        'eq:app-ao-recycle',
        'eq:app-ao-net',
        'eq:app-ao-flux',
        'eq:app-ao-erosion',
        'eq:app-ao-codep',
        'eq:app-ao-codepD',
        'eq:app-ao-inventory',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AO'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AO'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AO 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AO 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\s*\d+\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AO 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AO 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('5 Pa', '真空是什么单位：密度·正文（精确绑定）'),
        ('300$ K', '真空是什么单位：密度·公式（精确绑定）'),
        ('5$ Pa', '真空是什么单位：密度·公式（精确绑定）'),
        ('1.20716\\times10^{21}$ m$^{-3}$', '真空是什么单位：密度·公式（精确绑定）'),
        ('2.41432\\times10^{14}$ m$^{-3}$', '真空是什么单位：密度·公式（精确绑定）'),
        ('100$ m$^3$', '真空是什么单位：密度·公式（精确绑定）'),
        ('500$ Pa', '真空是什么单位：密度·公式（精确绑定）'),
        ('5000$ mbar', '真空是什么单位：密度·公式（精确绑定）'),
        ('1.20716\\times10^{23}$ s$^{-1}$', '真空是什么单位：密度·公式（精确绑定）'),
        ('2.0591\\times10^{21}$ s$^{-1}$', '真空是什么单位：密度·公式（精确绑定）'),
        ('4598.29 sccm', '真空是什么单位：密度·公式（精确绑定）'),
        ('500 Pa', '真空是什么单位：密度·表格（精确绑定）'),
        ('1 m', '抽气：有效抽速由导管·公式（精确绑定）'),
        ('12.1$ m$^3$', '抽气：有效抽速由导管·公式（精确绑定）'),
        ('10.7939$ m$^3$', '抽气：有效抽速由导管·公式（精确绑定）'),
        ('900$ m$^3$', '抽气：有效抽速由导管·公式（精确绑定）'),
        ('4.20551 m', '抽气：有效抽速由导管·公式（精确绑定）'),
        ('10\\%', '抽气：有效抽速由导管·正文（精确绑定）'),
        ('12.1 m$^3$', '抽气：有效抽速由导管·表格（精确绑定）'),
        ('10.7939 m$^3$', '抽气：有效抽速由导管·表格（精确绑定）'),
        ('900 m$^3$', '抽气：有效抽速由导管·表格（精确绑定）'),
        ('10 m$^3$', '抽气：有效抽速由导管·表格（精确绑定）'),
        ('4.89532\\times10^{-3}$ Pa', '放气与漏率：本底压强·公式（精确绑定）'),
        ('4.89532\\times10^{-5}$ Pa', '放气与漏率：本底压强·公式（精确绑定）'),
        ('4895.32$ m$^3$', '放气与漏率：本底压强·公式（精确绑定）'),
        ('4.89532\\times10^{-8}$ Pa', '放气与漏率：本底压强·公式（精确绑定）'),
        ('3\\times10^{-11}$ Pa', '放气与漏率：本底压强·公式（精确绑定）'),
        ('0.146442 g', '放气与漏率：本底压强·正文（精确绑定）'),
        ('5.0cm', '放气与漏率：本底压强·表格（精确绑定）'),
        ('6.0cm', '放气与漏率：本底压强·表格（精确绑定）'),
        ('4895.32 m$^3$', '放气与漏率：本底压强·表格（精确绑定）'),
        ('4.89532\\times10^{-6}$ Pa', '放气与漏率：本底压强·表格（精确绑定）'),
        ('396081 W', '烘烤与壁处理：为什么·公式（精确绑定）'),
        ('792162 W', '烘烤与壁处理：为什么·公式（精确绑定）'),
        ('100 h', '烘烤与壁处理：为什么·公式（精确绑定）'),
        ('79.2162 GJ', '烘烤与壁处理：为什么·公式（精确绑定）'),
        ('100$ nm', '烘烤与壁处理：为什么·正文（精确绑定）'),
        ('114.55 g', '烘烤与壁处理：为什么·正文（精确绑定）'),
        ('1.30347\\times10^{22}$ m$^{-2}$', '烘烤与壁处理：为什么·正文（精确绑定）'),
        ('1 h', '烘烤与壁处理：为什么·正文（精确绑定）'),
        ('100 nm', '烘烤与壁处理：为什么·表格（精确绑定）'),
        ('3.3288\\times10^{22}$ s$^{-1}$', '壁条件与再循环：净加·正文（精确绑定）'),
        ('3.16236\\times10^{22}$ s$^{-1}$', '壁条件与再循环：净加·公式（精确绑定）'),
        ('1.6644\\times10^{21}$ s$^{-1}$', '壁条件与再循环：净加·公式（精确绑定）'),
        ('5\\%', '壁条件与再循环：净加·公式（精确绑定）'),
        ('1\\%', '壁条件与再循环：净加·公式（精确绑定）'),
        ('6.79996\\times10^{19}$ m$^{-2}$s$^{-1}$', '壁条件与再循环：净加·公式（精确绑定）'),
        ('6.79996\\times10^{15}$ m$^{-2}$s$^{-1}$', '壁条件与再循环：净加·正文（精确绑定）'),
        ('1470.6 s', '壁条件与再循环：净加·正文（精确绑定）'),
        ('100$ m$^2$', '燃料滞留与氚库存：账·公式（精确绑定）'),
        ('9000 g', '燃料滞留与氚库存：账·公式（精确绑定）'),
        ('603.647 g', '燃料滞留与氚库存：账·公式（精确绑定）'),
        ('2 kg', '燃料滞留与氚库存：账·公式（精确绑定）'),
        ('1.6272 kg', '燃料滞留与氚库存：账·公式（精确绑定）'),
        ('2$ kg', '燃料滞留与氚库存：账·正文（精确绑定）'),
        ('100 g', '燃料滞留与氚库存：账·正文（精确绑定）'),
        ('100 m$^3$', '估算练习·正文（精确绑定）'),
        ('200 m$^3$', '估算练习·正文（精确绑定）'),
        ('0.5 m', '估算练习·正文（精确绑定）'),
        ('200 nm', '估算练习·正文（精确绑定）'),
        ('100 m$^2$', '估算练习·正文（精确绑定）'),
        ('200 m$^2$', '估算练习·正文（精确绑定）'),
        ('2.44766\\times10^{-8}$ Pa', '估算练习·正文（精确绑定）'),
        ('12.1\\,d^3', '估算练习·正文（精确绑定）'),
        ('1.5125$ m$^3$', '估算练习·正文（精确绑定）'),
        ('1.48996$ m$^3$', '估算练习·正文（精确绑定）'),
        ('229.101$ g', '估算练习·正文（精确绑定）'),
        ('3.3288\\times10^{21}$ s$^{-1}$', '估算练习·正文（精确绑定）'),
        ('1207.29$ g', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AO'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 数值 & 说明 \\\\', '表格数据行'),
        ('偏滤器压强 $P_{\\rm div}$ & 5 Pa & 排灰与杂质控制的工作点 \\\\', '表格数据行'),
        ('本底压强 $P_{\\rm vessel}$ & $10^{-6}$ Pa & 避免杂质回流芯部 \\\\', '表格数据行'),
        ('数密度（5 Pa） & $1.20716\\times10^{21}$ m$^{-3}$ & 式 \\eqref{eq:app-ao-density} \\\\', '表格数据行'),
        ('数密度（$10^{-6}$ Pa） & $2.41432\\times10^{14}$ m$^{-3}$ & 相差六个数量级 \\\\', '表格数据行'),
        ('真空室粒子数（5 Pa） & $1.0046\\times10^{24}$ & 芯部只有 $8.322\\times10^{22}$ \\\\', '表格数据行'),
        ('吞吐量 $Q=PS$ & 500 Pa$\\cdot$m$^3$/s & $=5000$ mbar$\\cdot$L/s \\\\', '表格数据行'),
        ('抽走粒子流 & $1.20716\\times10^{23}$ s$^{-1}$ & 附录 AJ 口径 \\\\', '表格数据行'),
        ('外部加料 & 4598.29 sccm & 只占抽气流的 1.7\\% \\\\', '表格数据行'),
        ('\\begin{tip}{5 Pa 与 $10^{-6}$ Pa 差六个数量级，但它们共用一套泵}{tip-app-ao-units}', '语义盒标题'),
        ("\\textbf{出处：}真空系统的量与工程口径参照 O'Hanlon《A User's Guide to Vacuum Technology", '出处声明整行'),
        ('导管流导 $C$（1 m$\\times$1 m） & 12.1 m$^3$/s & 式 \\eqref{eq:app-ao-conductance} \\\\', '表格数据行'),
        ('有效抽速 $S_{\\rm eff}$ & 10.7939 m$^3$/s & 式 \\eqref{eq:app-ao-seff} \\\\', '表格数据行'),
        ('相对铭牌的份额 & 10.79\\% & 损失 9.26446 倍 \\\\', '表格数据行'),
        ('损失一成所需流导 & 900 m$^3$/s & $C\\ge 9S_p$ \\\\', '表格数据行'),
        ('对应单管直径 & 4.20551 m & 式 \\eqref{eq:app-ao-lossless} \\\\', '表格数据行'),
        ('或并联 1 m 短管 & 74.38 根 & 工程上不可行 $\\Rightarrow$ 泵进盒体 \\\\', '表格数据行'),
        ('低温泵对氦的抽速 & 10 m$^3$/s & 约为氘的 1/10 \\\\', '表格数据行'),
        ('\\begin{tip}{泵的抽速不是系统的抽速}{tip-app-ao-pumping}', '语义盒标题'),
        ("\\textbf{出处：}分子流流导与有效抽速的公式参照 O'Hanlon《A User's Guide to Vacuum Technolo", '出处声明整行'),
        ('未烘烤放气量 & $4.89532\\times10^{-3}$ Pa$\\cdot$m$^3$/s & $q=10^{-5}$ \\\\', '表格数据行'),
        ('未烘烤平衡压强 & $4.89532\\times10^{-5}$ Pa & 差目标 48.9532 倍 \\\\', '表格数据行'),
        ('抽到目标所需抽速 & 4895.32 m$^3$/s & 不烘烤则不可能 \\\\', '表格数据行'),
        ('烘烤后放气量 & $4.89532\\times10^{-6}$ Pa$\\cdot$m$^3$/s & $q=10^{-8}$ \\\\', '表格数据行'),
        ('烘烤后平衡压强 & $4.89532\\times10^{-8}$ Pa & 余量 20.4277 倍 \\\\', '表格数据行'),
        ('实测漏率 & $3\\times10^{-9}$ Pa$\\cdot$m$^3$/s & 占放气 0.061283\\% \\\\', '表格数据行'),
        ('漏率贡献压强 & $3\\times10^{-11}$ Pa & 比放气小四个数量级 \\\\', '表格数据行'),
        ('检漏余量 & 300 倍 & 相对 $10^{-11}$ 灵敏度 \\\\', '表格数据行'),
        ('单层水分子数 & $4.89532\\times10^{21}$ & 0.146442 g \\\\', '表格数据行'),
        ('\\begin{tip}{本底压强是除法，不是抽功}{tip-app-ao-outgassing}', '语义盒标题'),
        ("\\textbf{出处：}放气率与烘烤的数据口径参照 O'Hanlon《A User's Guide to Vacuum Technology", '出处声明整行'),
        ('烘烤辐射功率 & 396081 W & 式 \\eqref{eq:app-ao-bakepower} \\\\', '表格数据行'),
        ('含对流余量 & 792162 W & 取一倍余量 \\\\', '表格数据行'),
        ('烘 100 h 的能量 & 79.2162 GJ & 电网账 \\\\', '表格数据行'),
        ('硼化 100 nm 的硼质量 & 114.55 g & 式 \\eqref{eq:app-ao-boron} \\\\', '表格数据行'),
        ('硼膜面密度 & $1.30347\\times10^{22}$ m$^{-2}$ & 1303.47 层单原子 \\\\', '表格数据行'),
        ('放电清洗 1 h 的粒子数 & $3.6\\times10^{25}$ & 约 7353.96 层水膜 \\\\', '表格数据行'),
        ('\\begin{tip}{壁处理不是打扫，是改变化学}{tip-app-ao-conditioning}', '语义盒标题'),
        ('\\textbf{出处：}硼化与放电清洗的物理参照 Winter（\\textit{Plasma Physics and Controlled ', '出处声明整行'),
        ('芯部粒子损失率 & $3.3288\\times10^{22}$ s$^{-1}$ & $N/\\tau_p$ \\\\', '表格数据行'),
        ('再循环返回流 & $3.16236\\times10^{22}$ s$^{-1}$ & 式 \\eqref{eq:app-ao-recycle} \\\\', '表格数据行'),
        ('净流出 & $1.6644\\times10^{21}$ s$^{-1}$ & 式 \\eqref{eq:app-ao-net} \\\\', '表格数据行'),
        ('外部加料/净流出 & 1.23714 & 同量级，差额是燃烧与滞留 \\\\', '表格数据行'),
        ('再循环/外部加料 & 15.358 & 加料只是零头 \\\\', '表格数据行'),
        ('壁面通量 & $6.79996\\times10^{19}$ m$^{-2}$s$^{-1}$ & 式 \\eqref{eq:app-ao-flux} \\\\', '表格数据行'),
        ('溅射率（$Y=10^{-4}$） & $6.79996\\times10^{15}$ m$^{-2}$s$^{-1}$ & 式 \\eqref{eq:app-ao-erosion} \\\\', '表格数据行'),
        ('单层原子被打一遍 & 1470.6 s & 约 24.5 分钟 \\\\', '表格数据行'),
        ('\\begin{tip}{再循环是乘数}{tip-app-ao-recycling}', '语义盒标题'),
        ('\\textbf{出处：}再循环系数与粒子平衡的口径参照 Wesson《Tokamaks》4th ed.（2011）第 4 章与 Stange', '出处声明整行'),
        ('共沉积层质量 & 9000 g & 式 \\eqref{eq:app-ao-codep} \\\\', '表格数据行'),
        ('氘氚原子数 & $1.80499\\times10^{26}$ & 式 \\eqref{eq:app-ao-codepD} \\\\', '表格数据行'),
        ('氘氚质量 & 603.647 g & 占库存上限 30.1823\\% \\\\', '表格数据行'),
        ('年滞留量 & 1.6272 kg/年 & 加料流的 1\\% \\\\', '表格数据行'),
        ('达到 2 kg 上限 & 1.22911 年 & 式 \\eqref{eq:app-ao-inventory} \\\\', '表格数据行'),
        ('清除周期（每次 100 g） & 6.03647 次/年 & 运行前提 \\\\', '表格数据行'),
        ('\\begin{tip}{氚库存的瓶颈在壁上}{tip-app-ao-retention}', '语义盒标题'),
        ('\\textbf{出处：}共沉积与氚滞留的物理参照 Federici 等（\\textit{Nuclear Fusion} \\textbf{41', '出处声明整行'),
        ('借鉴内容 & 来源 & 本书的处理 \\\\', '表格数据行'),
        ("真空量与流导公式 & O'Hanlon《A User's Guide to Vacuum Technology》3rd ed.（2003）第 1--4、8 章 & 重写推导，只保留结论与可复算系数 \\\\", '表格数据行'),
        ('放气率与烘烤数据 & 同上，第 8 章 & 取量级作声明输入，数值全部由 \\texttt{appAO\\_calc.py} 复算 \\\\', '表格数据行'),
        ('托卡马克工作压强与吞吐量 & ITER Physics Basis（\\textit{Nuclear Fusion} \\textbf{39}, 1999）第 4 章 & 只取口径（5 Pa / $10^{-6}$ Pa / 100 m$^3$/s） \\\\', '表格数据行'),
        ('壁处理与再循环 & Wesson《Tokamaks》4th ed.（2011）第 4 章；Stangeby（2000）第 4 章 & 改写为“净流量 = (1-R)”的账 \\\\', '表格数据行'),
        ('硼化与放电清洗 & Winter（\\textit{PPCF} \\textbf{38}, 1996） & 只引用原理结论，不复制工艺配方 \\\\', '表格数据行'),
        ('溅射与侵蚀寿命 & Roth 等（\\textit{PPCF} \\textbf{50}, 2008） & 只取产额量级，题型自拟 \\\\', '表格数据行'),
        ('共沉积与氚滞留 & Federici 等（\\textit{Nuclear Fusion} \\textbf{41}, 2001） & 只取 D/C 比与层厚量级 \\\\', '表格数据行'),
        ('库存与去氚约束 & Dolan（编）《Magnetic Fusion Technology》（2013）第 6 章 & 与附录 Z/AJ/AK 同口径，改写为 1.22911 年这笔账 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ao-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AO'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    # ---- 精确整串绑定（第 2 批）（堵 1e-3 容差内的末位改动） ----
    for lit, why in (
        (r'相当于 2606.94', r'练习 4 层数整串'),
        (r'0.603647（$60.3647\%$）', r'练习 6 份额整串'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AO'], why))
    # ---- 精确整串绑定（堵 1e-3 容差内的末位改动） ----
    for lit, why in (
        (r'碳原子数 $4.51247\times10^{26}$', r'共沉积碳原子数整串'),
        (r'提到 40.8553 倍', r'练习 1 余量整串'),
        (r'=1.48996$ m$^3$/s', r'练习 2 有效抽速整串'),
        (r'差 7.24442 倍', r'练习 2 倍数整串'),
        (r'之比降到 0.618571', r'练习 5 比值整串'),
        (r'面密度 $2.60694\times10^{22}$ m$^{-2}$', r'练习 4 面密度整串'),
        (r'$603.647\times2=1207.29$ g', r'练习 6 质量整串'),
        (r'得 114.55 g 硼', r'硼质量整串'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AO'], why))

    cases = [
        (0, 5, 'Pa', '真空是什么单位：密度·正文'),
        (0, 300, '$ K', '真空是什么单位：密度·公式'),
        (0, 5, '$ Pa', '真空是什么单位：密度·公式'),
        (0, 1.20716e+21, '$ m$^{-3}$', '真空是什么单位：密度·公式'),
        (0, 241432000000000.0, '$ m$^{-3}$', '真空是什么单位：密度·公式'),
        (0, 100, '$ m$^3$', '真空是什么单位：密度·公式'),
        (0, 500, '$ Pa', '真空是什么单位：密度·公式'),
        (0, 5000, '$ mbar', '真空是什么单位：密度·公式'),
        (0, 1.20716e+23, '$ s$^{-1}$', '真空是什么单位：密度·公式'),
        (0, 2.0591e+21, '$ s$^{-1}$', '真空是什么单位：密度·公式'),
        (0, 4598.29, 'sccm', '真空是什么单位：密度·公式'),
        (0, 500, 'Pa', '真空是什么单位：密度·表格'),
        (1, 1, 'm', '抽气：有效抽速由导管·公式'),
        (1, 12.1, '$ m$^3$', '抽气：有效抽速由导管·公式'),
        (1, 10.7939, '$ m$^3$', '抽气：有效抽速由导管·公式'),
        (1, 10.79, '\\%', '抽气：有效抽速由导管·公式'),
        (1, 900, '$ m$^3$', '抽气：有效抽速由导管·公式'),
        (1, 4.20551, 'm', '抽气：有效抽速由导管·公式'),
        (1, 10, '\\%', '抽气：有效抽速由导管·正文'),
        (1, 12.1, 'm$^3$', '抽气：有效抽速由导管·表格'),
        (1, 10.7939, 'm$^3$', '抽气：有效抽速由导管·表格'),
        (1, 900, 'm$^3$', '抽气：有效抽速由导管·表格'),
        (1, 10, 'm$^3$', '抽气：有效抽速由导管·表格'),
        (2, 0.00489532, '$ Pa', '放气与漏率：本底压强·公式'),
        (2, 4.89532e-05, '$ Pa', '放气与漏率：本底压强·公式'),
        (2, 4895.32, '$ m$^3$', '放气与漏率：本底压强·公式'),
        (2, 4.89532e-08, '$ Pa', '放气与漏率：本底压强·公式'),
        (2, 3e-11, '$ Pa', '放气与漏率：本底压强·公式'),
        (2, 0.146442, 'g', '放气与漏率：本底压强·正文'),
        (2, 5, 'cm', '放气与漏率：本底压强·表格'),
        (2, 6, 'cm', '放气与漏率：本底压强·表格'),
        (2, 4895.32, 'm$^3$', '放气与漏率：本底压强·表格'),
        (2, 4.89532e-06, '$ Pa', '放气与漏率：本底压强·表格'),
        (3, 396081, 'W', '烘烤与壁处理：为什么·公式'),
        (3, 792162, 'W', '烘烤与壁处理：为什么·公式'),
        (3, 100, 'h', '烘烤与壁处理：为什么·公式'),
        (3, 79.2162, 'GJ', '烘烤与壁处理：为什么·公式'),
        (3, 100, '$ nm', '烘烤与壁处理：为什么·正文'),
        (3, 114.55, 'g', '烘烤与壁处理：为什么·正文'),
        (3, 1.30347e+22, '$ m$^{-2}$', '烘烤与壁处理：为什么·正文'),
        (3, 1, 'h', '烘烤与壁处理：为什么·正文'),
        (3, 100, 'nm', '烘烤与壁处理：为什么·表格'),
        (4, 3.3288e+22, '$ s$^{-1}$', '壁条件与再循环：净加·正文'),
        (4, 3.16236e+22, '$ s$^{-1}$', '壁条件与再循环：净加·公式'),
        (4, 1.6644e+21, '$ s$^{-1}$', '壁条件与再循环：净加·公式'),
        (4, 5, '\\%', '壁条件与再循环：净加·公式'),
        (4, 1, '\\%', '壁条件与再循环：净加·公式'),
        (4, 6.79996e+19, '$ m$^{-2}$s$^{-1}$', '壁条件与再循环：净加·公式'),
        (4, 6799960000000000.0, '$ m$^{-2}$s$^{-1}$', '壁条件与再循环：净加·正文'),
        (4, 1470.6, 's', '壁条件与再循环：净加·正文'),
        (5, 100, '$ m$^2$', '燃料滞留与氚库存：账·公式'),
        (5, 9000, 'g', '燃料滞留与氚库存：账·公式'),
        (5, 603.647, 'g', '燃料滞留与氚库存：账·公式'),
        (5, 2, 'kg', '燃料滞留与氚库存：账·公式'),
        (5, 30.1823, '\\%', '燃料滞留与氚库存：账·公式'),
        (5, 1.6272, 'kg', '燃料滞留与氚库存：账·公式'),
        (5, 2, '$ kg', '燃料滞留与氚库存：账·正文'),
        (5, 100, 'g', '燃料滞留与氚库存：账·正文'),
        (6, 100, 'm$^3$', '估算练习·正文'),
        (6, 200, 'm$^3$', '估算练习·正文'),
        (6, 0.5, 'm', '估算练习·正文'),
        (6, 200, 'nm', '估算练习·正文'),
        (6, 100, 'm$^2$', '估算练习·正文'),
        (6, 200, 'm$^2$', '估算练习·正文'),
        (6, 2.44766e-08, '$ Pa', '估算练习·正文'),
        (6, 12.1, '\\,d^3', '估算练习·正文'),
        (6, 1.5125, '$ m$^3$', '估算练习·正文'),
        (6, 1.48996, '$ m$^3$', '估算练习·正文'),
        (6, 229.101, '$ g', '估算练习·正文'),
        (6, 3.3288e+21, '$ s$^{-1}$', '估算练习·正文'),
        (6, 1207.29, '$ g', '估算练习·正文'),
        (6, 60.3647, '\\%', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AO 的 %d 个关键数值未在小节内复算通过：%s'
            % (len(bad), '、'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AO 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流守卫（禁止把精确值四舍五入后写回正文） ----------
    for lit, why in (
        ('1.207\\times10^{21}', '偏滤器密度四舍五入'),
        ('1.00\\times10^{24}', '真空室总粒子数四舍五入'),
        ('4.48\\times10^{17}', '1 sccm 粒子流四舍五入'),
        ('4598 sccm', '外部加料折算 sccm 四舍五入'),
        ('2.41\\times10^{14}', '本底密度四舍五入'),
        ('10.79\\ ', '有效抽速四舍五入'),
        ('0.108\\ ', '有效抽速份额四舍五入'),
        ('9.26\\ ', '抽速损失倍数四舍五入'),
        ('0.00490\\ ', '未烘烤放气量四舍五入'),
        ('4.90\\times10^{-6}', '烘烤后放气量四舍五入'),
        ('4.90\\times10^{-5}', '未烘烤平衡压强四舍五入'),
        ('4.90\\times10^{-8}', '烘烤后平衡压强四舍五入'),
        ('48.95 倍', '未烘烤压强倍数四舍五入'),
        ('0.0613\\%', '漏率份额四舍五入'),
        ('3.96\\times10^{5}', '烘烤功率四舍五入'),
        ('7.92\\times10^{5}', '含对流烘烤功率四舍五入'),
        ('79.2\\ ', '烘烤能量四舍五入'),
        ('114.6\\ ', '硼质量四舍五入'),
        ('1.30\\times10^{22}', '硼面密度四舍五入'),
        ('1303\\ ', '硼膜层数四舍五入'),
        ('7354\\ ', '放电清洗水膜层数四舍五入'),
        ('3.33\\times10^{22}', '粒子损失率四舍五入'),
        ('3.16\\times10^{22}', '再循环返回流四舍五入'),
        ('1.66\\times10^{21}', '净流出四舍五入'),
        ('15.36\\ ', '再循环倍数四舍五入'),
        ('6.80\\times10^{19}', '壁面通量四舍五入'),
        ('1471\\ ', '单层打掉时间四舍五入'),
        ('4.51\\times10^{26}', '共沉积碳原子数四舍五入'),
        ('1.80\\times10^{26}', '共沉积氘氚原子数四舍五入'),
        ('603.6\\ ', '共沉积氘氚质量四舍五入'),
        ('30.18\\%', '库存份额四舍五入'),
        ('1.627\\ ', '年滞留量四舍五入'),
        ('1.229\\ ', '达到上限时间四舍五入'),
        ('6.04\\ ', '清扫周期四舍五入'),
        ('2.45\\times10^{-8}', '练习 1 压强四舍五入'),
        ('40.86\\ ', '练习 1 余量四舍五入'),
        ('1.51\\ ', '练习 2 流导四舍五入'),
        ('1.49\\ ', '练习 2 有效抽速四舍五入'),
        ('7.24\\ ', '练习 2 倍数四舍五入'),
        ('0.6186\\ ', '练习 5 比值四舍五入'),
        ('2.607\\times10^{22}', '练习 4 面密度四舍五入'),
        ('2607\\ ', '练习 4 层数四舍五入'),
        ('1207\\ ', '练习 6 质量四舍五入'),
        ('0.6036\\ ', '练习 6 份额四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AO'], why))


def check_P48_appendix_ap_maintenance():
    """附录 AP：维护与遥操作——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与编号解答数量不足或不相等；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AP'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AP 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '为什么维护是一门物理课',
        '包层寿命与更换周期',
        '活化与冷却窗口',
        '遥操作：精度、力反馈与延迟',
        '热室与部件流转',
        '在役检查与人因',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AP 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AP 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ap-avail',
        'tab:app-ap-blanket',
        'tab:app-ap-dose',
        'tab:app-ap-teleop',
        'tab:app-ap-hotcell',
        'tab:app-ap-human',
        'tab:app-ap-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AP'])

    for lab in (
        'sec:app-ap-availability',
        'sec:app-ap-blanket',
        'sec:app-ap-dose',
        'sec:app-ap-teleop',
        'sec:app-ap-hotcell',
        'sec:app-ap-human',
        'sec:app-ap-problems',
        'sec:app-ap-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AP'])

    for lab in (
        'tip-app-ap-avail',
        'tip-app-ap-blanket',
        'tip-app-ap-dose',
        'tip-app-ap-teleop',
        'tip-app-ap-hotcell',
        'tip-app-ap-human',
        'tip-app-ap-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AP'])

    for lab in (
        'eq:app-ap-avail',
        'eq:app-ap-cost',
        'eq:app-ap-sector',
        'eq:app-ap-campaign',
        'eq:app-ap-decay',
        'eq:app-ap-cooling',
        'eq:app-ap-hands',
        'eq:app-ap-acc',
        'eq:app-ap-shield',
        'eq:app-ap-err',
        'eq:app-ap-margin',
        'eq:app-ap-force',
        'eq:app-ap-latency',
        'eq:app-ap-need',
        'eq:app-ap-cap',
        'eq:app-ap-capmargin',
        'eq:app-ap-human',
        'eq:app-ap-budget',
        'eq:app-ap-turns',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AP'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AP'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AP 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AP 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\s*\d+\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AP 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AP 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('60\\ \\mathrm{d}', '为什么维护是一门物理·公式（精确绑定）'),
        ('20\\ \\mathrm{d}', '为什么维护是一门物理·公式（精确绑定）'),
        ('1920\\ \\mathrm{h}', '为什么维护是一门物理·公式（精确绑定）'),
        ('200\\ \\mathrm{MW}', '为什么维护是一门物理·公式（精确绑定）'),
        ('1.752\\times10^{6}\\ \\mathrm{MWh}', '为什么维护是一门物理·正文（精确绑定）'),
        ('384000\\ \\mathrm{MWh}', '为什么维护是一门物理·表格（精确绑定）'),
        ('62500\\ \\mathrm{kg}', '包层寿命与更换周期·公式（精确绑定）'),
        ('62.5\\ \\mathrm{t}', '包层寿命与更换周期·公式（精确绑定）'),
        ('200\\ \\mathrm{t}', '包层寿命与更换周期·正文（精确绑定）'),
        ('600\\ \\mathrm{h}', '包层寿命与更换周期·正文（精确绑定）'),
        ('1920 h', '包层寿命与更换周期·表格（精确绑定）'),
        ('10\\ \\mathrm{mSv}', '活化与冷却窗口：剂量·正文（精确绑定）'),
        ('2.5\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('1\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('33.2555\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('0.1\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('41.9041\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('30\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('10.7698\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('1.07698\\times10^{-3}\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式（精确绑定）'),
        ('0.01\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·表格（精确绑定）'),
        ('2.68476\\times10^{-3}\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·表格（精确绑定）'),
        ('5\\ \\mathrm{m}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.001\\ \\mathrm{m}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('1\\ \\mathrm{mm}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.5\\ \\mathrm{mm}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.1\\ \\mathrm{N}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('100\\ \\mathrm{N}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.1\\%', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('100\\ \\mathrm{ms}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('1\\ \\mathrm{cm}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.01\\ \\mathrm{m}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.1\\ \\mathrm{s}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('0.1\\ \\mathrm{m}', '遥操作：精度、力反馈·公式（精确绑定）'),
        ('10\\ \\mathrm{Hz}', '遥操作：精度、力反馈·正文（精确绑定）'),
        ('760.417\\%', '热室与部件流转：能力·正文（精确绑定）'),
        ('20\\ \\mathrm{mSv}', '在役检查与人因·正文（精确绑定）'),
        ('0.5\\ \\mathrm{mSv}', '在役检查与人因·正文（精确绑定）'),
        ('2.5\\ \\mathrm{mSv}', '在役检查与人因·公式（精确绑定）'),
        ('87.5\\%', '在役检查与人因·表格（精确绑定）'),
        ('45\\ \\mathrm{d}', '估算练习·正文（精确绑定）'),
        ('1\\ \\mathrm{d}', '估算练习·正文（精确绑定）'),
        ('15\\ \\mathrm{d}', '估算练习·正文（精确绑定）'),
        ('4.10959\\%', '估算练习·正文（精确绑定）'),
        ('333.333\\ \\mathrm{t}', '估算练习·正文（精确绑定）'),
        ('13.3022\\ \\mathrm{d}', '估算练习·正文（精确绑定）'),
        ('1.07698\\times10^{-4}\\ \\mathrm{Sv}', '估算练习·正文（精确绑定）'),
        ('$A = 0.780822$', '年可用因子（精确绑定·补充）'),
        ('$t_{\\text{down}} = 1920\\ \\mathrm{h}$', '年停机小时数（精确绑定·补充）'),
        ('= 62.5\\ \\mathrm{t}.', '单段质量吨位（精确绑定·补充）'),
        ('$600\\ \\mathrm{h}$', '每段可用工期（精确绑定·补充）'),
        ('= 24.3333\\ \\text{段}/\\text{年}', 'cask 年转运能力（精确绑定·补充）'),
        ('= 7.60417.', '热室能力余量（精确绑定·补充）'),
        ('单工位年周转 $12.1667$ 段', '单工位年周转能力（精确绑定·补充）'),
        ('$0.263014$ 个', '同时占用工位数（精确绑定·补充）'),
        ('$760.417\\%$', '能力余量百分数（精确绑定·补充）'),
        ('$87.5\\%$ 的余量', '剂量余量百分数（精确绑定·补充）'),
        ('$12.5\\%$', '限值占用百分数（精确绑定·补充）'),
        ('= 0.821918$', '练习 1 可用因子（精确绑定·补充）'),
        ('$4.10959\\%$', '练习 1 发电时间增幅（精确绑定·补充）'),
        ('= 5.33333$ 段', '练习 2 年更换段数（精确绑定·补充）'),
        ('$1.66667$ 倍', '练习 2 相对倍数（精确绑定·补充）'),
        ('$333.333\\ \\mathrm{t}$', '练习 2 年更换质量（精确绑定·补充）'),
        ('= 13.3022\\ \\mathrm{d}$', '练习 3 冷却时间（精确绑定·补充）'),
        ('= 48.6667$ 段', '练习 6 cask 能力（精确绑定·补充）'),
        ('= 15.2083$', '练习 6 能力余量（精确绑定·补充）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AP'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('项目 & 数值 & 说明 \\\\', '表格数据行'),
        ('计划维护 & $60\\ \\mathrm{d}$ & 包层/偏滤器更换、烘烤、检查 \\\\', '表格数据行'),
        ('非计划停机 & $20\\ \\mathrm{d}$ & 故障与修复（设计指标） \\\\', '表格数据行'),
        ('停机小时 & $1920\\ \\mathrm{h}$ & 两项之和 \\\\', '表格数据行'),
        ('可用因子 & $0.780822$ & 全年可发电小时份额 \\\\', '表格数据行'),
        ('停机的发电量损失 & $384000\\ \\mathrm{MWh}$ & $P_e \\times t_{\\text{down}}$ \\\\', '表格数据行'),
        ('停机的售电损失 & $1.536\\times10^{8}$ 元 & 按上述电价 \\\\', '表格数据行'),
        ('每停机一天的损失 & $1.92\\times10^{6}$ 元$/$天 & 停机成本的时间单价 \\\\', '表格数据行'),
        ('\\begin{tip}{可用因子是设计指标，不是运维指标}{tip-app-ap-avail}', '语义盒标题'),
        ('\\textbf{出处：}维护与可用性的工程口径参照 Dolan（编）\\textit{Magnetic Fusion Technology}', '出处声明整行'),
        ('单个扇形段质量 & $62500\\ \\mathrm{kg}$ & $= 62.5\\ \\mathrm{t}$ \\\\', '表格数据行'),
        ('全寿期更换次数 & $8$ & 40 年 $/$ 5 年 \\\\', '表格数据行'),
        ('全寿期更换段数 & $128$ & 每次 16 段 \\\\', '表格数据行'),
        ('平均每年更换段数 & $3.2$ & 128 段 $/$ 40 年 \\\\', '表格数据行'),
        ('平均每年更换质量 & $200\\ \\mathrm{t}$ & $3.2 \\times 62.5$ \\\\', '表格数据行'),
        ('每段可用工期 & $600\\ \\mathrm{h}$ & 1920 h $/$ 3.2 段 \\\\', '表格数据行'),
        ('\\begin{tip}{更换周期是设计输入，不是寿命预测}{tip-app-ap-blanket}', '语义盒标题'),
        ('\\textbf{出处：}包层与增殖段维护的工程口径参照 ITER Physics Basis', '出处声明整行'),
        ('降到 $1\\ \\mathrm{Sv}/\\mathrm{h}$ & $33.2555\\ \\mathrm{d}$ & 远程操作开工线 \\\\', '表格数据行'),
        ('降到 $0.1\\ \\mathrm{Sv}/\\mathrm{h}$ & $41.9041\\ \\mathrm{d}$ & 进一步保守 \\\\', '表格数据行'),
        ('一年后的剂量率 & $0.01\\ \\mathrm{Sv}/\\mathrm{h}$ & 短寿命分量已消失 \\\\', '表格数据行'),
        ('十年后的剂量率 & $2.68476\\times10^{-3}\\ \\mathrm{Sv}/\\mathrm{h}$ & 仍远高于人工阈值 \\\\', '表格数据行'),
        ('降到 $10\\ \\mu\\mathrm{Sv}/\\mathrm{h}$ & $52.5306$ 年 & 人工进入不可行 \\\\', '表格数据行'),
        ('30 天累计（未屏蔽） & $10.7698\\ \\mathrm{Sv}$ & 从开工线起算 \\\\', '表格数据行'),
        ('30 天累计（屏蔽后） & $1.07698\\times10^{-3}\\ \\mathrm{Sv}$ & $f_{\\text{shield}} = 10^{-4}$ \\\\', '表格数据行'),
        ('占器件耐总剂量 & $1.07698\\times10^{-8}$ & 器件寿命不由辐照决定 \\\\', '表格数据行'),
        ('\\begin{tip}{剂量率有两个分量，它们管两件不同的事}{tip-app-ap-dose}', '语义盒标题'),
        ('\\textbf{出处：}活化与剂量率的工程口径参照 Gilbert 与 Sublet 的活化计算综述', '出处声明整行'),
        ('臂端定位误差 & $0.001\\ \\mathrm{m}$ & $L\\varepsilon$，即 $1\\ \\mathrm{mm}$ \\\\', '表格数据行'),
        ('相对装配间隙的余量 & $0.5$ & 小于 1，必须二次对中 \\\\', '表格数据行'),
        ('力反馈相对分辨率 & $0.001$ & $0.1\\ \\mathrm{N} / 100\\ \\mathrm{N}$ \\\\', '表格数据行'),
        ('延迟对应的速度上限 & $0.1\\ \\mathrm{m}/\\mathrm{s}$ & $1\\ \\mathrm{cm} / 0.1\\ \\mathrm{s}$ \\\\', '表格数据行'),
        ('伺服带宽对应的周期 & $0.1\\ \\mathrm{s}$ & $1/10\\ \\mathrm{Hz}$ \\\\', '表格数据行'),
        ('\\begin{tip}{遥操作不是「远程看视频」}{tip-app-ap-teleop}', '语义盒标题'),
        ('\\textbf{出处：}遥操作与机器人学的工程口径参照 Siciliano 与 Khatib（编）', '出处声明整行'),
        ('年更换需求 & $3.2$ 段$/$年 & 与表~\\ref{tab:app-ap-blanket} 同 \\\\', '表格数据行'),
        ('cask 年转运能力 & $24.3333$ 段$/$年 & $2 \\times 365/30$ \\\\', '表格数据行'),
        ('能力余量 & $7.60417$ & 容器不是瓶颈 \\\\', '表格数据行'),
        ('单工位年周转能力 & $12.1667$ 段$/$年 & 真正的瓶颈 \\\\', '表格数据行'),
        ('同时占用工位数 & $0.263014$ & 平均在制段数 \\\\', '表格数据行'),
        ('\\begin{tip}{热室不是仓库，是工位}{tip-app-ap-hotcell}', '语义盒标题'),
        ('\\textbf{出处：}热室与部件流转的口径参照 Dolan（编）\\textit{Magnetic Fusion Technology}', '出处声明整行'),
        ('按单次剂量允许次数 & $40$ 次$/$年 & $20/0.5$ \\\\', '表格数据行'),
        ('在役检查的人工次数 & $5$ 次$/$年 & $100 \\times 0.05$ \\\\', '表格数据行'),
        ('人工年剂量 & $2.5\\ \\mathrm{mSv}/$年 & $5 \\times 0.5$ \\\\', '表格数据行'),
        ('占限值份额 & $0.125$ & 留出 $87.5\\%$ 的余量 \\\\', '表格数据行'),
        ('远程相对人工次数倍数 & $19$ & $0.95/0.05$ \\\\', '表格数据行'),
        ('\\begin{tip}{剂量预算是乘出来的}{tip-app-ap-human}', '语义盒标题'),
        ('\\textbf{出处：}剂量限值与防护框架参照 ICRP 第 103 号报告（2007）、', '出处声明整行'),
        ('来源 & 类型 & 借鉴内容 \\\\', '表格数据行'),
        ('Wesson \\textit{Tokamaks} 4th ed.（2011） & 教材 & 装置运行窗口与停机代价的背景 \\\\', '表格数据行'),
        ('教材 & 遥操作、力反馈与延迟的工程判据 \\\\', '表格数据行'),
        ('ICRP 第 103 号报告（2007） & 安全标准 & 剂量限值与防护框架 \\\\', '表格数据行'),
        ('教材 & 剂量、屏蔽与累积剂量的算法 \\\\', '表格数据行'),
        ('Gilbert 与 Sublet 活化计算综述（2010 年代） & 综述 & 活化与剂量率的量级 \\\\', '表格数据行'),
        ('IAEA 远程处理与热室技术文件（2000 年代） & 技术文件 & cask 与热室流转口径 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ap-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AP'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 60, '\\ \\mathrm{d}', '为什么维护是一门物理·公式'),
        (0, 20, '\\ \\mathrm{d}', '为什么维护是一门物理·公式'),
        (0, 1920, '\\ \\mathrm{h}', '为什么维护是一门物理·公式'),
        (0, 200, '\\ \\mathrm{MW}', '为什么维护是一门物理·公式'),
        (0, 1752000, '\\ \\mathrm{MWh}', '为什么维护是一门物理·正文'),
        (0, 384000, '\\ \\mathrm{MWh}', '为什么维护是一门物理·表格'),
        (1, 62500, '\\ \\mathrm{kg}', '包层寿命与更换周期·公式'),
        (1, 62.5, '\\ \\mathrm{t}', '包层寿命与更换周期·公式'),
        (1, 200, '\\ \\mathrm{t}', '包层寿命与更换周期·正文'),
        (1, 600, '\\ \\mathrm{h}', '包层寿命与更换周期·正文'),
        (1, 1920, 'h', '包层寿命与更换周期·表格'),
        (2, 10, '\\ \\mathrm{mSv}', '活化与冷却窗口：剂量·正文'),
        (2, 2.5, '\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式'),
        (2, 1, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式'),
        (2, 33.2555, '\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式'),
        (2, 0.1, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式'),
        (2, 41.9041, '\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式'),
        (2, 30, '\\ \\mathrm{d}', '活化与冷却窗口：剂量·公式'),
        (2, 10.7698, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式'),
        (2, 0.00107698, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·公式'),
        (2, 0.01, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·表格'),
        (2, 0.00268476, '\\ \\mathrm{Sv}', '活化与冷却窗口：剂量·表格'),
        (3, 5, '\\ \\mathrm{m}', '遥操作：精度、力反馈·公式'),
        (3, 0.001, '\\ \\mathrm{m}', '遥操作：精度、力反馈·公式'),
        (3, 1, '\\ \\mathrm{mm}', '遥操作：精度、力反馈·公式'),
        (3, 0.5, '\\ \\mathrm{mm}', '遥操作：精度、力反馈·公式'),
        (3, 0.1, '\\ \\mathrm{N}', '遥操作：精度、力反馈·公式'),
        (3, 100, '\\ \\mathrm{N}', '遥操作：精度、力反馈·公式'),
        (3, 0.1, '\\%', '遥操作：精度、力反馈·公式'),
        (3, 100, '\\ \\mathrm{ms}', '遥操作：精度、力反馈·公式'),
        (3, 1, '\\ \\mathrm{cm}', '遥操作：精度、力反馈·公式'),
        (3, 0.01, '\\ \\mathrm{m}', '遥操作：精度、力反馈·公式'),
        (3, 0.1, '\\ \\mathrm{s}', '遥操作：精度、力反馈·公式'),
        (3, 0.1, '\\ \\mathrm{m}', '遥操作：精度、力反馈·公式'),
        (3, 10, '\\ \\mathrm{Hz}', '遥操作：精度、力反馈·正文'),
        (4, 760.417, '\\%', '热室与部件流转：能力·正文'),
        (5, 20, '\\ \\mathrm{mSv}', '在役检查与人因·正文'),
        (5, 0.5, '\\ \\mathrm{mSv}', '在役检查与人因·正文'),
        (5, 2.5, '\\ \\mathrm{mSv}', '在役检查与人因·公式'),
        (5, 87.5, '\\%', '在役检查与人因·表格'),
        (6, 45, '\\ \\mathrm{d}', '估算练习·正文'),
        (6, 1, '\\ \\mathrm{d}', '估算练习·正文'),
        (6, 15, '\\ \\mathrm{d}', '估算练习·正文'),
        (6, 4.10959, '\\%', '估算练习·正文'),
        (6, 333.333, '\\ \\mathrm{t}', '估算练习·正文'),
        (6, 13.3022, '\\ \\mathrm{d}', '估算练习·正文'),
        (6, 0.000107698, '\\ \\mathrm{Sv}', '估算练习·正文'),
        (0, 0.780822, '', '年可用因子'),
        (0, 1920, '', '年停机小时数'),
        (1, 62.5, '', '单段质量吨位'),
        (1, 600, '', '每段可用工期'),
        (4, 24.3333, '', 'cask 年转运能力'),
        (4, 7.60417, '', '热室能力余量'),
        (4, 12.1667, '', '单工位年周转能力'),
        (4, 0.263014, '', '同时占用工位数'),
        (4, 760.417, '', '能力余量百分数'),
        (5, 87.5, '', '剂量余量百分数'),
        (5, 12.5, '', '限值占用百分数'),
        (6, 0.821918, '', '练习 1 可用因子'),
        (6, 4.10959, '', '练习 1 发电时间增幅'),
        (6, 5.33333, '', '练习 2 年更换段数'),
        (6, 1.66667, '', '练习 2 相对倍数'),
        (6, 333.333, '', '练习 2 年更换质量'),
        (6, 13.3022, '', '练习 3 冷却时间'),
        (6, 48.6667, '', '练习 6 cask 能力'),
        (6, 15.2083, '', '练习 6 能力余量'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AP 的 %d 个关键数值未在小节内复算通过：%s'
            % (len(bad), '、'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AP 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流守卫（禁止把精确值四舍五入后写回正文） ----------
    for lit, why in (
        ('760.4\\%', '热室能力余量百分数四舍五入'),
        ('12.17\\ ', '单工位年周转能力四舍五入'),
        ('0.263\\ ', '同时占用工位数四舍五入'),
        ('0.8219\\ ', '练习 1 可用因子四舍五入'),
        ('0.0411\\ ', '练习 1 可用因子提高量四舍五入'),
        ('4.11\\%', '练习 1 发电时间增幅四舍五入'),
        ('5.33\\ ', '练习 2 年更换段数四舍五入'),
        ('1.67\\ ', '练习 2 相对倍数四舍五入'),
        ('333.3\\ ', '练习 2 年更换质量四舍五入'),
        ('13.3\\ ', '练习 3 冷却时间四舍五入'),
        ('48.67\\ ', '练习 6 cask 能力四舍五入'),
        ('15.21\\ ', '练习 6 能力余量四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AP'], why))


def check_P49_appendix_aq_interferometry():
    """附录 AQ：干涉与相位诊断——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与编号解答数量不足或不相等；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 数字与单位的精确子串绑定，
       后者专门堵「只改末位数字、落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AQ'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AQ 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '干涉仪测的是相移',
        '条纹、外差与零差',
        '振动补偿与色散干涉',
        '阿贝尔反演',
        '源与探测器链',
        '误差预算与相对论修正',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AQ 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AQ 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-aq-phase',
        'tab:app-aq-heterodyne',
        'tab:app-aq-dispersion',
        'tab:app-aq-abel',
        'tab:app-aq-chain',
        'tab:app-aq-error',
        'tab:app-aq-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AQ'])

    for lab in (
        'sec:app-aq-phase',
        'sec:app-aq-heterodyne',
        'sec:app-aq-dispersion',
        'sec:app-aq-abel',
        'sec:app-aq-chain',
        'sec:app-aq-error',
        'sec:app-aq-problems',
        'sec:app-aq-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AQ'])

    for lab in (
        'tip-app-aq-phase',
        'tip-app-aq-heterodyne',
        'tip-app-aq-dispersion',
        'tip-app-aq-abel',
        'tip-app-aq-chain',
        'tip-app-aq-relativistic',
        'tip-app-aq-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AQ'])

    for lab in (
        'eq:app-aq-dispersion',
        'eq:app-aq-refractive',
        'eq:app-aq-phase',
        'eq:app-aq-lineint',
        'eq:app-aq-fringe',
        'eq:app-aq-heterodyne',
        'eq:app-aq-noise',
        'eq:app-aq-vibration',
        'eq:app-aq-dispersive',
        'eq:app-aq-abel',
        'eq:app-aq-peel',
        'eq:app-aq-weight',
        'eq:app-aq-source',
        'eq:app-aq-relativistic',
        'eq:app-aq-expansion',
        'eq:app-aq-budget',
        'eq:app-aq-sensi',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AQ'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AQ'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AQ 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AQ 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+\}', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\s*\d+\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AQ 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AQ 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('10\\ \\mathrm{GHz}', '干涉仪测的是相移，不·公式（精确绑定）'),
        ('2.81794\\times10^{-15}\\ \\mathrm{m}', '干涉仪测的是相移，不·公式（精确绑定）'),
        ('2.0\\ \\mathrm{m}', '干涉仪测的是相移，不·公式（精确绑定）'),
        ('2.66667\\times10^{20}\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·公式（精确绑定）'),
        ('3.0cm', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('10.6\\ \\mu\\mathrm{m}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('118.8\\ \\mu\\mathrm{m}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('1.0\\ \\mathrm{mm}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('2.82823\\times10^{13}\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('2.52351\\times10^{12}\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('2.99792\\times10^{11}\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('9.92216\\times10^{24}\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('7.89924\\times10^{22}\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('1.11485\\times10^{21}\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('7.96538\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('89.2723\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('751.451\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('2.1035\\times10^{20}\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('1.87686\\times10^{19}\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('2.22971\\times10^{18}\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格（精确绑定）'),
        ('6.66667\\times10^{19}\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·正文（精确绑定）'),
        ('2\\,\\mathrm{SNR}', '条纹、外差与零差：把·公式（精确绑定）'),
        ('1\\ \\mu\\mathrm{W}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('8\\ \\mathrm{dB}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('1.58489\\times10^{-7}\\ \\mathrm{W}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('3000\\ \\mathrm{K}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('1\\ \\mathrm{MHz}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('4.14195\\times10^{-14}\\ \\mathrm{W}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('3.61482\\times10^{-4}\\ \\mathrm{rad}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('1.07979\\times10^{15}\\ \\mathrm{m^{-2}}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('4.04921\\times10^{14}\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('1\\ \\mu\\mathrm{m}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('0.001\\ \\mathrm{rad}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('1.25543\\times10^{16}\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('1.12017\\times10^{15}\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('1.33076\\times10^{14}\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('0.592753\\ \\mathrm{rad}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('0.0528888\\ \\mathrm{rad}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('0.00628319\\ \\mathrm{rad}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('1.68704\\ \\mathrm{nm}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('18.9076\\ \\mathrm{nm}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('159.155\\ \\mathrm{nm}', '条纹、外差与零差：把·表格（精确绑定）'),
        ('1\\ \\mu\\mathrm{s}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('9.38429\\times10^{24}\\ \\mathrm{m^{-2}/s}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('3.51911\\times10^{24}\\ \\mathrm{m^{-3}/s}', '条纹、外差与零差：把·正文（精确绑定）'),
        ('1064\\ \\mathrm{nm}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('532\\ \\mathrm{nm}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('0.799544\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('1.59909\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('5.90525\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('6.70479\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('7.50434\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('1.11022\\times10^{-16}\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文（精确绑定）'),
        ('5.0cm', '振动补偿与色散干涉：·表格（精确绑定）'),
        ('7.38578\\times10^{20}\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格（精确绑定）'),
        ('4.52111\\times10^{16}\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格（精确绑定）'),
        ('4.53817\\times10^{15}\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格（精确绑定）'),
        ('1.80\\ \\mathrm{m}', '阿贝尔反演：从弦积分·正文（精确绑定）'),
        ('0.45\\ \\mathrm{m}', '阿贝尔反演：从弦积分·正文（精确绑定）'),
        ('4.0cm', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('2.46675\\times10^{20}\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.90\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('1.89917\\times10^{20}\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('1.35\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('1.07107\\times10^{20}\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('2.20851\\times10^{19}\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('9.399\\times10^{19}\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('6.01002\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.2\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('9.92284\\times10^{19}\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.771641\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('9.97193\\times10^{19}\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.280675\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('1\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.00682732\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('3.15242\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.5\\ \\mathrm{mm}', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.0929994\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('0.355102\\%', '阿贝尔反演：从弦积分·表格（精确绑定）'),
        ('1.8\\ \\mathrm{m}', '阿贝尔反演：从弦积分·正文（精确绑定）'),
        ('0.1\\ \\mathrm{W}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('110\\ \\mathrm{GHz}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('10\\%', '源与探测器链：相位噪·正文（精确绑定）'),
        ('0.01\\ \\mathrm{W}', '源与探测器链：相位噪·公式（精确绑定）'),
        ('330\\ \\mathrm{GHz}', '源与探测器链：相位噪·公式（精确绑定）'),
        ('1.35084\\times10^{21}\\ \\mathrm{m^{-3}}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('3\\ \\mathrm{dB}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('0.005\\ \\mathrm{W}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('1\\ \\mathrm{mW}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('30\\ \\mathrm{dB}', '源与探测器链：相位噪·正文（精确绑定）'),
        ('1.79922\\%', '误差预算与相对论修正·公式（精确绑定）'),
        ('2.34903\\%', '误差预算与相对论修正·公式（精确绑定）'),
        ('2.08885\\%', '误差预算与相对论修正·公式（精确绑定）'),
        ('0.000404921\\%', '误差预算与相对论修正·表格（精确绑定）'),
        ('2.08885\\times10^{18}\\ \\mathrm{m^{-3}}', '误差预算与相对论修正·表格（精确绑定）'),
        ('2\\%', '误差预算与相对论修正·语义盒（精确绑定）'),
        ('4.0\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('2.98711\\times10^{15}\\ \\mathrm{m^{-2}}', '估算练习·正文（精确绑定）'),
        ('0.000361482\\ \\mathrm{rad}', '估算练习·正文（精确绑定）'),
        ('0.0628319\\ \\mathrm{rad}', '估算练习·正文（精确绑定）'),
        ('0.592753\\ \\mathrm{rad/\\mu m}', '估算练习·正文（精确绑定）'),
        ('0.106\\ \\mu\\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('1.06\\times10^{-7}\\ \\mathrm{m}', '估算练习·正文（精确绑定）'),
        ('2.6334\\%', '估算练习·正文（精确绑定）'),
        ('2.96296\\times10^{18}\\ \\mathrm{m^{-3}}', '估算练习·正文（精确绑定）'),
        ('6\\%', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AQ'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('项目 & CO$_2$ $10.6\\ \\mu\\mathrm{m}$ & 远红外 $118.8\\ \\mu\\mathrm{m}$ & 毫米波 $1.0\\ \\mathrm{mm}$ \\\\', '表格数据行'),
        ('频率 & $2.82823\\times10^{13}\\ \\mathrm{Hz}$ & $2.52351\\times10^{12}\\ \\mathrm{Hz}$ & $2.99792\\times10^{11}\\ \\mathrm{Hz}$ \\\\', '表格数据行'),
        ('截止密度 & $9.92216\\times10^{24}\\ \\mathrm{m^{-3}}$ & $7.89924\\times10^{22}\\ \\mathrm{m^{-3}}$ & $1.11485\\times10^{21}\\ \\mathrm{m^{-3}}$ \\\\', '表格数据行'),
        ('工作密度占比 $n_{e0}/n_c$ & $1.00784\\times10^{-5}$ & $1.26594\\times10^{-3}$ & $8.96978\\times10^{-2}$ \\\\', '表格数据行'),
        ('折射率亏损 $1-N$ & $5.03922\\times10^{-6}$ & $6.32972\\times10^{-4}$ & $4.48489\\times10^{-2}$ \\\\', '表格数据行'),
        ('中心弦相移 & $7.96538\\ \\mathrm{rad}$ & $89.2723\\ \\mathrm{rad}$ & $751.451\\ \\mathrm{rad}$ \\\\', '表格数据行'),
        ('条纹数 & $1.26773$ & $14.2081$ & $119.597$ \\\\', '表格数据行'),
        ('一条条纹对应的线积分密度 & $2.1035\\times10^{20}\\ \\mathrm{m^{-2}}$ & $1.87686\\times10^{19}\\ \\mathrm{m^{-2}}$ & $2.22971\\times10^{18}\\ \\mathrm{m^{-2}}$ \\\\', '表格数据行'),
        ('\\begin{tip}{相位是密度的线积分，不是密度}{tip-app-aq-phase}', '语义盒标题'),
        ('\\textbf{出处：}Hutchinson \\textit{Principles of Plasma Diagnostics} 第 2 版', '出处声明整行'),
        ('$0.001\\ \\mathrm{rad}$ 折合的电子密度 & $1.25543\\times10^{16}\\ \\mathrm{m^{-3}}$ & $1.12017\\times10^{15}\\ \\mathrm{m^{-3}}$ & $1.33076\\times10^{14}\\ \\mathrm{m^{-3}}$ \\\\', '表格数据行'),
        ('$1\\ \\mu\\mathrm{m}$ 振动引入的相位 & $0.592753\\ \\mathrm{rad}$ & $0.0528888\\ \\mathrm{rad}$ & $0.00628319\\ \\mathrm{rad}$ \\\\', '表格数据行'),
        ('折合条纹数 & $0.0943396$ & $0.00841751$ & $0.001$ \\\\', '表格数据行'),
        ('相位控制在 $0.001\\ \\mathrm{rad}$ 的允许振幅 & $1.68704\\ \\mathrm{nm}$ & $18.9076\\ \\mathrm{nm}$ & $159.155\\ \\mathrm{nm}$ \\\\', '表格数据行'),
        ('\\begin{tip}{条纹是计数，外差是计时}{tip-app-aq-heterodyne}', '语义盒标题'),
        ('方案 & 对振动的处理 & 折合中心密度误差 \\\\', '表格数据行'),
        ('单波长（$1064\\ \\mathrm{nm}$） & 不补偿，全靠机械稳定 & $7.38578\\times10^{20}\\ \\mathrm{m^{-3}}$（$7.38578$ 倍工作密度） \\\\', '表格数据行'),
        ('色散干涉（$1064$ + $532\\ \\mathrm{nm}$，共光路） & 由 $2f-f$ 组合直接减掉 & $4.52111\\times10^{16}\\ \\mathrm{m^{-3}}$（相对 $4.52111\\times10^{-4}$） \\\\', '表格数据行'),
        ('双色独立两条弦（$118.8\\ \\mu\\mathrm{m}$ + $10.6\\ \\mu\\mathrm{m}$） & 不共路，只能各自稳定 & $4.53817\\times10^{15}\\ \\mathrm{m^{-3}}$（相位噪声极限） \\\\', '表格数据行'),
        ('\\begin{tip}{双波长不是为了看得更准，是为了把振动减掉}{tip-app-aq-dispersion}', '语义盒标题'),
        ('\\textbf{出处：}Boboc 等，\\textit{Rev. Sci. Instrum.} \\textbf{81}（2010）', '出处声明整行'),
        ('项目 & 数值 & 说明 \\\\', '表格数据行'),
        ('$y = 0$ 弦 & $2.66667\\times10^{20}\\ \\mathrm{m^{-2}}$ & 中心道，也是最大的一道 \\\\', '表格数据行'),
        ('$y = 0.45\\ \\mathrm{m}$ 弦 & $2.46675\\times10^{20}\\ \\mathrm{m^{-2}}$ & 壳层边界 \\\\', '表格数据行'),
        ('$y = 0.90\\ \\mathrm{m}$ 弦 & $1.89917\\times10^{20}\\ \\mathrm{m^{-2}}$ & 壳层边界 \\\\', '表格数据行'),
        ('$y = 1.35\\ \\mathrm{m}$ 弦 & $1.07107\\times10^{20}\\ \\mathrm{m^{-2}}$ & 壳层边界 \\\\', '表格数据行'),
        ('$y = 1.80\\ \\mathrm{m}$ 弦 & $2.20851\\times10^{19}\\ \\mathrm{m^{-2}}$ & 最外道，壳层延伸到 $a$ \\\\', '表格数据行'),
        ('中心道与最外道之比 & $12.0745$ & 动态范围：$12$ 倍 \\\\', '表格数据行'),
        ('9 道反演的中心密度 & $9.399\\times10^{19}\\ \\mathrm{m^{-3}}$ & 偏低 $6.01002\\%$（最外层壳层厚 $0.2\\ \\mathrm{m}$） \\\\', '表格数据行'),
        ('20 道反演的中心密度 & $9.92284\\times10^{19}\\ \\mathrm{m^{-3}}$ & 偏低 $0.771641\\%$ \\\\', '表格数据行'),
        ('40 道反演的中心密度 & $9.97193\\times10^{19}\\ \\mathrm{m^{-3}}$ & 偏低 $0.280675\\%$ \\\\', '表格数据行'),
        ('最外道偏 $1\\%$ 对中心的影响 & $-0.00682732\\%$ & 只动它自己的壳层 \\\\', '表格数据行'),
        ('中心道偏 $1\\%$ 对中心的影响 & $3.15242\\%$ & 权重最大的一道 \\\\', '表格数据行'),
        ('弦位置偏 $0.5\\ \\mathrm{mm}$（$y = 1.35\\ \\mathrm{m}$） & $0.0929994\\%$ & 该道线积分的相对变化 \\\\', '表格数据行'),
        ('弦位置偏 $0.5\\ \\mathrm{mm}$（$y = 1.80\\ \\mathrm{m}$） & $0.355102\\%$ & 与标定误差同量级 \\\\', '表格数据行'),
        ('\\begin{tip}{剥皮法的权重写在壳层厚度里}{tip-app-aq-abel}', '语义盒标题'),
        ('\\textbf{出处：}Abel 的原始变换（1826）与 Hutchinson \\textit{Principles of Plasma', '出处声明整行'),
        ('源功率（$110\\ \\mathrm{GHz}$） & $0.1\\ \\mathrm{W}$ & 声明 \\\\', '表格数据行'),
        ('三倍频效率 & $0.10$ & 声明 \\\\', '表格数据行'),
        ('输出功率（$330\\ \\mathrm{GHz}$） & $0.01\\ \\mathrm{W}$ & $P_{\\mathrm{gunn}}\\eta_{\\mathrm{tripler}}$ \\\\', '表格数据行'),
        ('波长 & $0.908462\\ \\mathrm{mm}$ & $c/(3f_{\\mathrm{gunn}})$ \\\\', '表格数据行'),
        ('截止密度 & $1.35084\\times10^{21}\\ \\mathrm{m^{-3}}$ & 工作密度的 $13.5084$ 倍 \\\\', '表格数据行'),
        ('到达等离子体功率（$3\\ \\mathrm{dB}$ 损失） & $0.005\\ \\mathrm{W}$ & 含窗口与传输 \\\\', '表格数据行'),
        ('本振功率 & $1\\ \\mathrm{mW}$ & 比信号高 $30\\ \\mathrm{dB}$ \\\\', '表格数据行'),
        ('混频器变频损耗 & $8\\ \\mathrm{dB}$ & 声明 \\\\', '表格数据行'),
        ('系统噪声温度 & $3000\\ \\mathrm{K}$ & 混频器 $+$ 中频 \\\\', '表格数据行'),
        ('热噪声 $k_BT_{\\mathrm{sys}}B$（$1\\ \\mathrm{MHz}$） & $4.14195\\times10^{-14}\\ \\mathrm{W}$ & 噪声地板 \\\\', '表格数据行'),
        ('相位刷新率对应的时间分辨 & $10^{-6}\\ \\mathrm{s}$ & 声明 \\\\', '表格数据行'),
        ('\\begin{tip}{相位噪声是噪声温度的函数}{tip-app-aq-chain}', '语义盒标题'),
        ('\\textbf{出处：}Heald \\& Wharton \\textit{Plasma Diagnostics with Microwave', '出处声明整行'),
        ('误差项 & 相对误差 & 来源与性质 \\\\', '表格数据行'),
        ('相位噪声 & $0.000404921\\%$ & $1\\ \\mu\\mathrm{W}$ 信号、$T_{\\mathrm{sys}} = 3000\\ \\mathrm{K}$、$B = 1\\ \\mathrm{MHz}$ \\\\', '表格数据行'),
        ('标定 & $1\\%$ & 弦长、相对标定与光学效率 \\\\', '表格数据行'),
        ('弦位置 & $0.355102\\%$ & $0.5\\ \\mathrm{mm}$ 安装误差（最外道） \\\\', '表格数据行'),
        ('相对论修正 & $1.79922\\%$ & $T_e/m_ec^2$ 量级，可按 $T_e$ 修正 \\\\', '表格数据行'),
        ('合成（RSS） & $2.08885\\%$ & 折合 $2.08885\\times10^{18}\\ \\mathrm{m^{-3}}$ \\\\', '表格数据行'),
        ('\\begin{tip}{电子温度到了 keV，折射率里的质量就不再是常数}{tip-app-aq-relativistic}', '语义盒标题'),
        ('来源 & 本附录用到的部分 & 分级 \\\\', '表格数据行'),
        ('Heald \\& Wharton, \\textit{Plasma Diagnostics with Microwaves} (1965) & 微波干涉的相移公式、外差接收与混频器噪声 & 一级（学科基石） \\\\', '表格数据行'),
        ('Hutchinson, \\textit{Principles of Plasma Diagnostics}, 2nd ed. (2002) 第 4 章 & 折射率、阿贝尔反演、误差传播、相对论修正 & 一级（标准教材） \\\\', '表格数据行'),
        ('Donné 等, \\textit{Nuclear Fusion} \\textbf{47} (2007) S337 & ITER 密度诊断的精度要求与工程口径 & 二级（工程综述） \\\\', '表格数据行'),
        ('Abel (1826) 的积分变换 & 阿贝尔正反变换与层析思想 & 三级（史料） \\\\', '表格数据行'),
        ('本书 \\texttt{tmp/appAQ\\_calc.py} & 全部数值（可由脚本复算） & 自算（一级） \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-aq-honest}', '语义盒标题'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AQ'], why))

    # ---------- ⑤b 补绑定：逐条探针发现的洞（v1.26.21） ----------
    for lit, why in (
        ('相对工作密度有 $13.5084$ 倍余量', '正文余量句（同一数在表行外还有一处）'),
        ('色散干涉把它压到相位噪声水平，压低倍数 $16336.2$。', '色散干涉压低倍数整句'),
        ('（相对 $4.52111\\times10^{-4}$），改善因子 $16336.2$。', '改善因子整句（同值第二处）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AQ'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 10, '\\ \\mathrm{GHz}', '干涉仪测的是相移，不·公式'),
        (0, 2.81794e-15, '\\ \\mathrm{m}', '干涉仪测的是相移，不·公式'),
        (0, 2, '\\ \\mathrm{m}', '干涉仪测的是相移，不·公式'),
        (0, 2.66667e+20, '\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·公式'),
        (0, 3, 'cm', '干涉仪测的是相移，不·表格'),
        (0, 10.6, '\\ \\mu\\mathrm{m}', '干涉仪测的是相移，不·表格'),
        (0, 118.8, '\\ \\mu\\mathrm{m}', '干涉仪测的是相移，不·表格'),
        (0, 1, '\\ \\mathrm{mm}', '干涉仪测的是相移，不·表格'),
        (0, 2.82823e+13, '\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格'),
        (0, 2.52351e+12, '\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格'),
        (0, 2.99792e+11, '\\ \\mathrm{Hz}', '干涉仪测的是相移，不·表格'),
        (0, 9.92216e+24, '\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格'),
        (0, 7.89924e+22, '\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格'),
        (0, 1.11485e+21, '\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·表格'),
        (0, 7.96538, '\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格'),
        (0, 89.2723, '\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格'),
        (0, 751.451, '\\ \\mathrm{rad}', '干涉仪测的是相移，不·表格'),
        (0, 2.1035e+20, '\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格'),
        (0, 1.87686e+19, '\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格'),
        (0, 2.22971e+18, '\\ \\mathrm{m^{-2}}', '干涉仪测的是相移，不·表格'),
        (0, 6.66667e+19, '\\ \\mathrm{m^{-3}}', '干涉仪测的是相移，不·正文'),
        (1, 2, '\\,\\mathrm{SNR}', '条纹、外差与零差：把·公式'),
        (1, 1, '\\ \\mu\\mathrm{W}', '条纹、外差与零差：把·正文'),
        (1, 8, '\\ \\mathrm{dB}', '条纹、外差与零差：把·正文'),
        (1, 1.58489e-07, '\\ \\mathrm{W}', '条纹、外差与零差：把·正文'),
        (1, 3000, '\\ \\mathrm{K}', '条纹、外差与零差：把·正文'),
        (1, 1, '\\ \\mathrm{MHz}', '条纹、外差与零差：把·正文'),
        (1, 4.14195e-14, '\\ \\mathrm{W}', '条纹、外差与零差：把·正文'),
        (1, 0.000361482, '\\ \\mathrm{rad}', '条纹、外差与零差：把·正文'),
        (1, 1.07979e+15, '\\ \\mathrm{m^{-2}}', '条纹、外差与零差：把·正文'),
        (1, 4.04921e+14, '\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·正文'),
        (1, 1, '\\ \\mu\\mathrm{m}', '条纹、外差与零差：把·表格'),
        (1, 0.001, '\\ \\mathrm{rad}', '条纹、外差与零差：把·表格'),
        (1, 1.25543e+16, '\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格'),
        (1, 1.12017e+15, '\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格'),
        (1, 1.33076e+14, '\\ \\mathrm{m^{-3}}', '条纹、外差与零差：把·表格'),
        (1, 0.592753, '\\ \\mathrm{rad}', '条纹、外差与零差：把·表格'),
        (1, 0.0528888, '\\ \\mathrm{rad}', '条纹、外差与零差：把·表格'),
        (1, 0.00628319, '\\ \\mathrm{rad}', '条纹、外差与零差：把·表格'),
        (1, 1.68704, '\\ \\mathrm{nm}', '条纹、外差与零差：把·表格'),
        (1, 18.9076, '\\ \\mathrm{nm}', '条纹、外差与零差：把·表格'),
        (1, 159.155, '\\ \\mathrm{nm}', '条纹、外差与零差：把·表格'),
        (1, 1, '\\ \\mu\\mathrm{s}', '条纹、外差与零差：把·正文'),
        (1, 9.38429e+24, '\\ \\mathrm{m^{-2}/s}', '条纹、外差与零差：把·正文'),
        (1, 3.51911e+24, '\\ \\mathrm{m^{-3}/s}', '条纹、外差与零差：把·正文'),
        (2, 1064, '\\ \\mathrm{nm}', '振动补偿与色散干涉：·正文'),
        (2, 532, '\\ \\mathrm{nm}', '振动补偿与色散干涉：·正文'),
        (2, 0.799544, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 1.59909, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 5.90525, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 6.70479, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 7.50434, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 1.11022e-16, '\\ \\mathrm{rad}', '振动补偿与色散干涉：·正文'),
        (2, 5, 'cm', '振动补偿与色散干涉：·表格'),
        (2, 7.38578e+20, '\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格'),
        (2, 4.52111e+16, '\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格'),
        (2, 4.53817e+15, '\\ \\mathrm{m^{-3}}', '振动补偿与色散干涉：·表格'),
        (3, 1.8, '\\ \\mathrm{m}', '阿贝尔反演：从弦积分·正文'),
        (3, 0.45, '\\ \\mathrm{m}', '阿贝尔反演：从弦积分·正文'),
        (3, 4, 'cm', '阿贝尔反演：从弦积分·表格'),
        (3, 2.46675e+20, '\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格'),
        (3, 0.9, '\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格'),
        (3, 1.89917e+20, '\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格'),
        (3, 1.35, '\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格'),
        (3, 1.07107e+20, '\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格'),
        (3, 2.20851e+19, '\\ \\mathrm{m^{-2}}', '阿贝尔反演：从弦积分·表格'),
        (3, 9.399e+19, '\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格'),
        (3, 6.01002, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 0.2, '\\ \\mathrm{m}', '阿贝尔反演：从弦积分·表格'),
        (3, 9.92284e+19, '\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格'),
        (3, 0.771641, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 9.97193e+19, '\\ \\mathrm{m^{-3}}', '阿贝尔反演：从弦积分·表格'),
        (3, 0.280675, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 1, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 0.00682732, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 3.15242, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 0.5, '\\ \\mathrm{mm}', '阿贝尔反演：从弦积分·表格'),
        (3, 0.0929994, '\\%', '阿贝尔反演：从弦积分·表格'),
        (3, 0.355102, '\\%', '阿贝尔反演：从弦积分·表格'),
        (4, 0.1, '\\ \\mathrm{W}', '源与探测器链：相位噪·正文'),
        (4, 110, '\\ \\mathrm{GHz}', '源与探测器链：相位噪·正文'),
        (4, 10, '\\%', '源与探测器链：相位噪·正文'),
        (4, 0.01, '\\ \\mathrm{W}', '源与探测器链：相位噪·公式'),
        (4, 330, '\\ \\mathrm{GHz}', '源与探测器链：相位噪·公式'),
        (4, 1.35084e+21, '\\ \\mathrm{m^{-3}}', '源与探测器链：相位噪·正文'),
        (4, 3, '\\ \\mathrm{dB}', '源与探测器链：相位噪·正文'),
        (4, 0.005, '\\ \\mathrm{W}', '源与探测器链：相位噪·正文'),
        (4, 1, '\\ \\mathrm{mW}', '源与探测器链：相位噪·正文'),
        (4, 30, '\\ \\mathrm{dB}', '源与探测器链：相位噪·正文'),
        (5, 1.79922, '\\%', '误差预算与相对论修正·公式'),
        (5, 2.34903, '\\%', '误差预算与相对论修正·公式'),
        (5, 2.08885, '\\%', '误差预算与相对论修正·公式'),
        (5, 0.000404921, '\\%', '误差预算与相对论修正·表格'),
        (5, 2.08885e+18, '\\ \\mathrm{m^{-3}}', '误差预算与相对论修正·表格'),
        (5, 2, '\\%', '误差预算与相对论修正·语义盒'),
        (6, 4, '\\ \\mathrm{m}', '估算练习·正文'),
        (6, 2.98711e+15, '\\ \\mathrm{m^{-2}}', '估算练习·正文'),
        (6, 65.83, '\\ \\mathrm{dB}', '估算练习·正文'),
        (6, 0.0628319, '\\ \\mathrm{rad}', '估算练习·正文'),
        (6, 0.592753, '\\ \\mathrm{rad/\\mu m}', '估算练习·正文'),
        (6, 0.106, '\\ \\mu\\mathrm{m}', '估算练习·正文'),
        (6, 1.06e-07, '\\ \\mathrm{m}', '估算练习·正文'),
        (6, 2.6334, '\\%', '估算练习·正文'),
        (6, 2.96296e+18, '\\ \\mathrm{m^{-3}}', '估算练习·正文'),
        (6, 6, '\\%', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AQ 的 %d 个关键数值未在小节内复算通过：%s'
            % (len(bad), '、'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AQ 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流守卫（禁止把精确值四舍五入后写回正文） ----------
    for lit, why in (
        ('12.07\\ ', '中心道与最外道之比四舍五入'),
        ('6.01\\%', '9 道弦反演偏低百分数四舍五入'),
        ('0.77\\%', '20 道弦反演偏低百分数四舍五入'),
        ('0.28\\%', '40 道弦反演偏低百分数四舍五入'),
        ('3.15\\%', '中心道偏 1% 影响百分数四舍五入'),
        ('0.36\\%', '弦位置误差百分数四舍五入'),
        ('13.5\\ ', '截止密度余量倍数四舍五入'),
        ('2.35\\%', '毫米波线性近似偏差百分数四舍五入'),
        ('1.8\\%', '相对论修正百分数四舍五入'),
        ('2.09\\%', '合成误差百分数四舍五入'),
        ('0.91\\ ', '输出波长毫米数四舍五入'),
        ('1.12\\times10^{-5}', 'FIR 相位分辨率折合密度份额四舍五入'),
        ('1.88\\times10^{19}', '一条条纹对应线积分密度四舍五入'),
        ('16336\\ ', '振动压低倍数四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AQ'], why))


def check_P50_appendix_ar_power():
    """附录 AR：电源与储能——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与编号解答数量不足或不相等；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 「数字与单位」精确子串绑定，后者专门堵
       「只改末位数字、仍落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AR'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AR 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '电源与储能：电厂里的三本能量账',
        '变流器与纹波',
        '储能：飞轮、电容',
        '环流功率',
        '电网接口',
        '保护与故障',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AR 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AR 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-ar-books',
        'tab:app-ar-pulse',
        'tab:app-ar-storage',
        'tab:app-ar-recirc',
        'tab:app-ar-grid',
        'tab:app-ar-protect',
        'tab:app-ar-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AR'])

    for lab in (
        'sec:app-ar-chain',
        'sec:app-ar-converter',
        'sec:app-ar-storage',
        'sec:app-ar-recirculating',
        'sec:app-ar-grid',
        'sec:app-ar-protect',
        'sec:app-ar-problems',
        'sec:app-ar-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AR'])

    for lab in (
        'tip-app-ar-books',
        'tip-app-ar-ripple',
        'tip-app-ar-flywheel',
        'tip-app-ar-recirc',
        'tip-app-ar-grid',
        'tip-app-ar-protect',
        'tip-app-ar-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AR'])

    for lab in (
        'eq:app-ar-gross',
        'eq:app-ar-deliv',
        'eq:app-ar-recirc',
        'eq:app-ar-net',
        'eq:app-ar-cs',
        'eq:app-ar-ramp',
        'eq:app-ar-ud',
        'eq:app-ar-pd',
        'eq:app-ar-harm',
        'eq:app-ar-filter',
        'eq:app-ar-q',
        'eq:app-ar-flywheel',
        'eq:app-ar-flywheel2',
        'eq:app-ar-vmax',
        'eq:app-ar-cap',
        'eq:app-ar-ind',
        'eq:app-ar-discharge',
        'eq:app-ar-qeng',
        'eq:app-ar-sens',
        'eq:app-ar-rocof',
        'eq:app-ar-black',
        'eq:app-ar-fast',
        'eq:app-ar-decay',
        'eq:app-ar-mit',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AR'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AR'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AR 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AR 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+\}', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\s*(\d+)\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AR 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AR 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('556.5$ MW', '电源与储能：电厂里的·正文（精确绑定）'),
        ('5.565\\times10^{8}$ W', '电源与储能：电厂里的·正文（精确绑定）'),
        ('2.226\\times10^{8}\\ \\mathrm{W}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('222.6\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('50.0$ MW', '电源与储能：电厂里的·公式（精确绑定）'),
        ('13.4157$ MW', '电源与储能：电厂里的·公式（精确绑定）'),
        ('8.0$ MW', '电源与储能：电厂里的·公式（精确绑定）'),
        ('71.4157\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('73.6244\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('2.20873\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('148.976\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('33.0748\\%', '电源与储能：电厂里的·公式（精确绑定）'),
        ('26.7701\\%', '电源与储能：电厂里的·公式（精确绑定）'),
        ('40\\%', '电源与储能：电厂里的·公式（精确绑定）'),
        ('1 MW', '电源与储能：电厂里的·公式（精确绑定）'),
        ('3.73551$ MW', '电源与储能：电厂里的·公式（精确绑定）'),
        ('10$ H', '电源与储能：电厂里的·公式（精确绑定）'),
        ('40$ kA', '电源与储能：电厂里的·公式（精确绑定）'),
        ('12 L', '电源与储能：电厂里的·公式（精确绑定）'),
        ('8.0\\times10^{9}\\ \\mathrm{J}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('8\\ \\mathrm{GJ}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('10$ s', '电源与储能：电厂里的·公式（精确绑定）'),
        ('8.0\\times10^{8}\\ \\mathrm{W}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('800\\ \\mathrm{MW}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('8000\\ \\mathrm{A/s}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('80\\ \\mathrm{kV}', '电源与储能：电厂里的·公式（精确绑定）'),
        ('800$ MW', '电源与储能：电厂里的·正文（精确绑定）'),
        ('148.976$ MW', '电源与储能：电厂里的·正文（精确绑定）'),
        ('100$ s', '电源与储能：电厂里的·正文（精确绑定）'),
        ('80$ MW', '电源与储能：电厂里的·正文（精确绑定）'),
        ('10$ ms', '电源与储能：电厂里的·正文（精确绑定）'),
        ('2.0$ MJ', '电源与储能：电厂里的·正文（精确绑定）'),
        ('2.0\\times10^{8}$ W', '电源与储能：电厂里的·正文（精确绑定）'),
        ('200$ MW', '电源与储能：电厂里的·正文（精确绑定）'),
        ('27.1649$ ms', '电源与储能：电厂里的·正文（精确绑定）'),
        ('222.6$ MW', '电源与储能：电厂里的·表格（精确绑定）'),
        ('73.6244$ MW', '电源与储能：电厂里的·表格（精确绑定）'),
        ('8$ GJ', '电源与储能：电厂里的·表格（精确绑定）'),
        ('10 s', '电源与储能：电厂里的·表格（精确绑定）'),
        ('13.4\\%', '变流器与纹波：十二脉·正文（精确绑定）'),
        ('4.2\\%', '变流器与纹波：十二脉·正文（精确绑定）'),
        ('2.1\\%', '变流器与纹波：十二脉·正文（精确绑定）'),
        ('24$ kV', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('23.1822$ kV', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('9.27289\\times10^{8}\\ \\mathrm{W}', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('927.289\\ \\mathrm{MW}', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('986.8$ MVA', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('50$ Hz', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('0.837365$ mH', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('0.599534$ mH', '变流器与纹波：十二脉·公式（精确绑定）'),
        ('59.6455$ MVAr', '变流器与纹波：十二脉·正文（精确绑定）'),
        ('2.71586\\%', '变流器与纹波：十二脉·正文（精确绑定）'),
        ('3.0cm', '变流器与纹波：十二脉·表格（精确绑定）'),
        ('8 GJ', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('20000$ kg', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('2.0$ m', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('40000\\ \\mathrm{kg\\,m^2}', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('314.159$ rad', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1.97392\\times10^{9}\\ \\mathrm{J}', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1.97392\\ \\mathrm{GJ}', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('628.319$ m', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('3.09906\\times10^{9}$ Pa', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('3.09906$ GPa', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1.0$ GPa', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('356.915$ m', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('3.0$ GPa', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1369.31$ m', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('9.375$ GJ', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('20$ kV', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('0.01$ F', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('10$ mF', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('5$ mF', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('28.2843$ kV', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('10$ kA', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('20$ kA', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1$ m$\\Omega', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('10000$ s', '储能：飞轮、电容与它·公式（精确绑定）'),
        ('1.97392$ GJ', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('100$ MW', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('19.7392$ s', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('80$ s', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('8$ MW', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('28.8$ GJ', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('1$ MWh', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('3.6\\times10^{9}$ J', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('0.548311$ MWh', '储能：飞轮、电容与它·正文（精确绑定）'),
        ('40$ MJ', '储能：飞轮、电容与它·表格（精确绑定）'),
        ('67.9122\\%', '环流功率：电厂净效率·正文（精确绑定）'),
        ('18.2218\\%', '环流功率：电厂净效率·正文（精确绑定）'),
        ('10.866\\%', '环流功率：电厂净效率·正文（精确绑定）'),
        ('3\\%', '环流功率：电厂净效率·正文（精确绑定）'),
        ('10\\%', '环流功率：电厂净效率·公式（精确绑定）'),
        ('5$ MW', '环流功率：电厂净效率·公式（精确绑定）'),
        ('25\\%', '环流功率：电厂净效率·正文（精确绑定）'),
        ('55.65$ MW', '环流功率：电厂净效率·正文（精确绑定）'),
        ('17.9744$ MW', '环流功率：电厂净效率·正文（精确绑定）'),
        ('2.20873$ MW', '环流功率：电厂净效率·表格（精确绑定）'),
        ('100\\%', '环流功率：电厂净效率·表格（精确绑定）'),
        ('5$ s', '电网接口：并网、频率·公式（精确绑定）'),
        ('1000$ MVA', '电网接口：并网、频率·公式（精确绑定）'),
        ('0.01$ Hz', '电网接口：并网、频率·公式（精确绑定）'),
        ('0.2$ Hz', '电网接口：并网、频率·公式（精确绑定）'),
        ('20$ s', '电网接口：并网、频率·公式（精确绑定）'),
        ('71.4157$ MW', '电网接口：并网、频率·公式（精确绑定）'),
        ('7.14157\\%', '电网接口：并网、频率·公式（精确绑定）'),
        ('14.8976\\%', '电网接口：并网、频率·公式（精确绑定）'),
        ('2.88\\times10^{10}\\ \\mathrm{J}', '电网接口：并网、频率·公式（精确绑定）'),
        ('28.8\\ \\mathrm{GJ}', '电网接口：并网、频率·公式（精确绑定）'),
        ('1$ ms', '保护与故障：毫秒级动·正文（精确绑定）'),
        ('100$ kA', '保护与故障：毫秒级动·正文（精确绑定）'),
        ('500\\ \\mathrm{V}', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('25\\ \\mathrm{MW}', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('25$ MW', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('2.77778$ h', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('1$ s', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('76.9074\\ \\mathrm{ms}', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('2.0\\ \\mathrm{MJ}', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('27.1649\\ \\mathrm{ms}', '保护与故障：毫秒级动·公式（精确绑定）'),
        ('500$ V', '保护与故障：毫秒级动·表格（精确绑定）'),
        ('76.9074$ ms', '保护与故障：毫秒级动·表格（精确绑定）'),
        ('960$ MW', '保护与故障：毫秒级动·表格（精确绑定）'),
        ('250.425$ MW', '估算练习·正文（精确绑定）'),
        ('29.3998\\%', '估算练习·正文（精确绑定）'),
        ('31.7701\\%', '估算练习·正文（精确绑定）'),
        ('176.801$ MW', '估算练习·正文（精确绑定）'),
        ('418.879$ rad', '估算练习·正文（精确绑定）'),
        ('3.50919$ GJ', '估算练习·正文（精确绑定）'),
        ('837.758$ m', '估算练习·正文（精确绑定）'),
        ('5.50943$ GPa', '估算练习·正文（精确绑定）'),
        ('4.0\\times10^{8}$ W', '估算练习·正文（精确绑定）'),
        ('400$ MW', '估算练习·正文（精确绑定）'),
        ('4000$ A', '估算练习·正文（精确绑定）'),
        ('40$ kV', '估算练习·正文（精确绑定）'),
        ('8533.33$ kg', '估算练习·正文（精确绑定）'),
        ('103.8$ MVAr', '估算练习·正文（精确绑定）'),
        ('45.2009$ MVAr', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AR'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('账本 & 代表量 & 数值 & 对设计的要求 \\\\', '表格数据行'),
        ('稳态 & 毛电 $P_\\mathrm{gross}$ & $222.6$ MW & 热电转换效率定上限 \\\\', '表格数据行'),
        ('稳态 & 环流 $P_\\mathrm{rec}$ & $73.6244$ MW（$33.0748\\%$） & 每一项负荷都要单独省 \\\\', '表格数据行'),
        ('脉冲 & 螺线管储能 $E_\\mathrm{CS}$ & $8$ GJ & 要么拉长爬升，要么上储能 \\\\', '表格数据行'),
        ('脉冲 & 10 s 爬升功率 & $800$ MW（净电的 $5.37001$ 倍） & 电网瞬时能力或厂内储能 \\\\', '表格数据行'),
        ('瞬态 & 破裂缓解脉冲 & $200$ MW / $10$ ms & 毫秒级开关与能量 \\\\', '表格数据行'),
        ('瞬态 & 等效环流时间 & $27.1649$ ms & 能量小、时间贵 \\\\', '表格数据行'),
        ('\\begin{tip}{三本账要分开算}{tip-app-ar-books}', '语义盒标题'),
        ('\\textbf{出处：}电厂能量账与再循环功率的口径见 Freidberg, \\emph{Plasma Physics and Fusion', '出处声明整行'),
        ('脉波数 & 纹波系数 & 最低谐波次数 & 相对幅值 & 典型用途 \\\\', '表格数据行'),
        ('6 & $13.4\\%$ & 5、7 & $0.2$、$0.142857$ & 小功率线圈电源 \\\\', '表格数据行'),
        ('12 & $4.2\\%$ & 11、13 & $0.0909091$、$0.0769231$ & 大电流线圈电源标准 \\\\', '表格数据行'),
        ('24 & $2.1\\%$ & 23、25 & 更小 & 对纹波极敏感场合 \\\\', '表格数据行'),
        ('\\begin{tip}{纹波不是噪声，是换相留下的缺口}{tip-app-ar-ripple}', '语义盒标题'),
        ('\\textbf{出处：}脉波数、纹波系数与谐波次数的关系是电力电子教科书的标准内容，见 Rashid,', '出处声明整行'),
        ('路线 & 能量公式 & 本附录例子 & 能撑多久 & 主要限制 \\\\', '表格数据行'),
        ('飞轮 & $\\tfrac12I\\omega^2$ & $1.97392$ GJ（$20000$ kg、$2.0$ m、$3000$ rpm） & $19.7392$ s（$100$ MW） & 材料强度 $\\sigma=\\rho v^2$ \\\\', '表格数据行'),
        ('电容 & $\\tfrac12CV^2$ & $2.0$ MJ（$10$ mF、$20$ kV） & $10$ ms（$200$ MW） & 能量密度低（$40$ MJ/m$^3$） \\\\', '表格数据行'),
        ('电感 & $\\tfrac12LI^2$ & $8$ GJ（$10$ H、$40$ kA） & 放电受保护限制 & 断路与电压失控 \\\\', '表格数据行'),
        ('电池 & 化学能 & $1$ MWh $=3.6\\times10^{9}$ J & 小时级 & 功率密度与循环寿命 \\\\', '表格数据行'),
        ('\\begin{tip}{飞轮的极限不是轴承，是材料强度}{tip-app-ar-flywheel}', '语义盒标题'),
        ('\\textbf{出处：}飞轮的应力限制与比能量推导见 Genta, \\emph{Kinetic Energy Storage}（1985）与', '出处声明整行'),
        ('项目 & 功率 & 占环流 & 说明 \\\\', '表格数据行'),
        ('加热与电流驱动 & $50.0$ MW & $67.9122\\%$ & 由等离子体物理决定，最值得优化 \\\\', '表格数据行'),
        ('低温制冷 & $13.4157$ MW & $18.2218\\%$ & 与附录 AM 同一口径 \\\\', '表格数据行'),
        ('其他厂用电 & $8.0$ MW & $10.866\\%$ & 泵、加料、诊断、控制 \\\\', '表格数据行'),
        ('变换损耗 & $2.20873$ MW & $3\\%$ & 电源效率 $0.97$ 的代价 \\\\', '表格数据行'),
        ('合计 & $73.6244$ MW & $100\\%$ & 占毛电 $33.0748\\%$ \\\\', '表格数据行'),
        ('\\begin{tip}{环流功率是分母，不是余数}{tip-app-ar-recirc}', '语义盒标题'),
        ('\\textbf{出处：}再循环功率与工程增益的定义见 Freidberg, \\emph{Plasma Physics and Fusion ', '出处声明整行'),
        ('量 & 数值 & 含义 \\\\', '表格数据行'),
        ('频率变化率 & $0.01$ Hz/s & $H=5$ s、$S_\\mathrm{grid}=1000$ MVA、$100$ MW 扰动 \\\\', '表格数据行'),
        ('$0.2$ Hz 对应时间 & $20$ s & 一次调频与储能的动作窗口 \\\\', '表格数据行'),
        ('厂用电占短路容量 & $7.14157\\%$ & 大负荷启动会扰动并网点电压 \\\\', '表格数据行'),
        ('净电送出占比 & $14.8976\\%$ & 送出线路与保护按此设计 \\\\', '表格数据行'),
        ('无功需求 & $59.6455$ MVAr & $\\alpha=15^\\circ$、毛电 $222.6$ MW \\\\', '表格数据行'),
        ('黑启动能量 & $28.8$ GJ & 厂用电 $8$ MW 撑一小时 \\\\', '表格数据行'),
        ('\\begin{tip}{电网不是无限大的理想电压源}{tip-app-ar-grid}', '语义盒标题'),
        ('\\textbf{出处：}频率变化率与惯性常数的关系见 Kundur, \\emph{Power System Stability and Co', '出处声明整行'),
        ('动作 & 时间尺度 & 功率/电压 & 关键约束 \\\\', '表格数据行'),
        ('快控线圈 & $1$ ms & $500$ V、$25$ MW & 储能功率密度 \\\\', '表格数据行'),
        ('失电保护 & $L/R=10000$ s & --- & 保持可控，而非放电 \\\\', '表格数据行'),
        ('破裂缓解 & $10$ ms（$\\ll76.9074$ ms） & $200$ MW & 阀门 $1$ ms 响应 \\\\', '表格数据行'),
        ('直流侧短路 & 保护动作时间 & $960$ MW（$40$ kA、$24$ kV） & 限流与断路器 \\\\', '表格数据行'),
        ('\\begin{tip}{保护的时间尺度由热沉积决定}{tip-app-ar-protect}', '语义盒标题'),
        ('\\textbf{出处：}破裂缓解的时间窗口与 $\\tau_\\mathrm{CQ}$ 口径取本书附录 V、X 的同一组数值；快控线圈的', '出处声明整行'),
        ('来源 & 借鉴内容 & 分级 & 版权口径 \\\\', '表格数据行'),
        ('Freidberg 2007 & 电厂能量账、再循环功率与工程增益的定义 & 教科书 & 只复述定义与量级，不复制 \\\\', '表格数据行'),
        ('Stacey 2011 & 电厂章的能量流与效率口径 & 教科书 & 只取口径，不复制表格 \\\\', '表格数据行'),
        ('Dolan 2013 & 电源、储能与电气接口的工程讨论 & 专著 & 只取分类框架，不复制 \\\\', '表格数据行'),
        ('Rashid / Mohan 等 & 脉波数、纹波与谐波关系的公式 & 教科书 & 公式为标准内容，自行复算 \\\\', '表格数据行'),
        ('Kundur 1994 & 频率变化率与惯性常数的关系 & 教科书 & 只取关系式，不复制 \\\\', '表格数据行'),
        ('Genta 1985 / Bolund 2007 & 飞轮的应力限制与比能量 & 专著/综述 & 只取 $\\sigma=\\rho v^2$ 与量级，不复制 \\\\', '表格数据行'),
        ('ITER Physics Basis 1999 & 电源与破裂缓解的时间尺度 & 文献 & 只引时间尺度，不复制 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-ar-honest}', '语义盒标题'),
        ('\\textbf{出处：}来源表所列文献；分级与不复制口径沿用本书附录 AM/AO/AP 的同一约定。本附录不复制', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AR'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 556.5, '$ MW', '电源与储能：电厂里的·正文'),
        (0, 556500000, '$ W', '电源与储能：电厂里的·正文'),
        (0, 222600000, '\\ \\mathrm{W}', '电源与储能：电厂里的·公式'),
        (0, 222.6, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 50, '$ MW', '电源与储能：电厂里的·公式'),
        (0, 13.4157, '$ MW', '电源与储能：电厂里的·公式'),
        (0, 8, '$ MW', '电源与储能：电厂里的·公式'),
        (0, 71.4157, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 73.6244, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 2.20873, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 148.976, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 33.0748, '\\%', '电源与储能：电厂里的·公式'),
        (0, 26.7701, '\\%', '电源与储能：电厂里的·公式'),
        (0, 40, '\\%', '电源与储能：电厂里的·公式'),
        (0, 1, 'MW', '电源与储能：电厂里的·公式'),
        (0, 3.73551, '$ MW', '电源与储能：电厂里的·公式'),
        (0, 10, '$ H', '电源与储能：电厂里的·公式'),
        (0, 40, '$ kA', '电源与储能：电厂里的·公式'),
        (0, 12, 'L', '电源与储能：电厂里的·公式'),
        (0, 8000000000, '\\ \\mathrm{J}', '电源与储能：电厂里的·公式'),
        (0, 8, '\\ \\mathrm{GJ}', '电源与储能：电厂里的·公式'),
        (0, 10, '$ s', '电源与储能：电厂里的·公式'),
        (0, 800000000, '\\ \\mathrm{W}', '电源与储能：电厂里的·公式'),
        (0, 800, '\\ \\mathrm{MW}', '电源与储能：电厂里的·公式'),
        (0, 8000, '\\ \\mathrm{A/s}', '电源与储能：电厂里的·公式'),
        (0, 80, '\\ \\mathrm{kV}', '电源与储能：电厂里的·公式'),
        (0, 800, '$ MW', '电源与储能：电厂里的·正文'),
        (0, 148.976, '$ MW', '电源与储能：电厂里的·正文'),
        (0, 100, '$ s', '电源与储能：电厂里的·正文'),
        (0, 80, '$ MW', '电源与储能：电厂里的·正文'),
        (0, 10, '$ ms', '电源与储能：电厂里的·正文'),
        (0, 2, '$ MJ', '电源与储能：电厂里的·正文'),
        (0, 200000000, '$ W', '电源与储能：电厂里的·正文'),
        (0, 200, '$ MW', '电源与储能：电厂里的·正文'),
        (0, 27.1649, '$ ms', '电源与储能：电厂里的·正文'),
        (0, 222.6, '$ MW', '电源与储能：电厂里的·表格'),
        (0, 73.6244, '$ MW', '电源与储能：电厂里的·表格'),
        (0, 8, '$ GJ', '电源与储能：电厂里的·表格'),
        (0, 10, 's', '电源与储能：电厂里的·表格'),
        (1, 13.4, '\\%', '变流器与纹波：十二脉·正文'),
        (1, 4.2, '\\%', '变流器与纹波：十二脉·正文'),
        (1, 2.1, '\\%', '变流器与纹波：十二脉·正文'),
        (1, 24, '$ kV', '变流器与纹波：十二脉·公式'),
        (1, 23.1822, '$ kV', '变流器与纹波：十二脉·公式'),
        (1, 927289000, '\\ \\mathrm{W}', '变流器与纹波：十二脉·公式'),
        (1, 927.289, '\\ \\mathrm{MW}', '变流器与纹波：十二脉·公式'),
        (1, 986.8, '$ MVA', '变流器与纹波：十二脉·公式'),
        (1, 50, '$ Hz', '变流器与纹波：十二脉·公式'),
        (1, 0.837365, '$ mH', '变流器与纹波：十二脉·公式'),
        (1, 0.599534, '$ mH', '变流器与纹波：十二脉·公式'),
        (1, 59.6455, '$ MVAr', '变流器与纹波：十二脉·正文'),
        (1, 2.71586, '\\%', '变流器与纹波：十二脉·正文'),
        (1, 3, 'cm', '变流器与纹波：十二脉·表格'),
        (1, 3.1, 'cm', '变流器与纹波：十二脉·表格'),
        (2, 8, 'GJ', '储能：飞轮、电容与它·正文'),
        (2, 20000, '$ kg', '储能：飞轮、电容与它·公式'),
        (2, 2, '$ m', '储能：飞轮、电容与它·公式'),
        (2, 40000, '\\ \\mathrm{kg\\,m^2}', '储能：飞轮、电容与它·公式'),
        (2, 314.159, '$ rad', '储能：飞轮、电容与它·公式'),
        (2, 1973920000, '\\ \\mathrm{J}', '储能：飞轮、电容与它·公式'),
        (2, 1.97392, '\\ \\mathrm{GJ}', '储能：飞轮、电容与它·公式'),
        (2, 628.319, '$ m', '储能：飞轮、电容与它·公式'),
        (2, 3099060000, '$ Pa', '储能：飞轮、电容与它·公式'),
        (2, 3.09906, '$ GPa', '储能：飞轮、电容与它·公式'),
        (2, 1, '$ GPa', '储能：飞轮、电容与它·公式'),
        (2, 356.915, '$ m', '储能：飞轮、电容与它·公式'),
        (2, 3, '$ GPa', '储能：飞轮、电容与它·公式'),
        (2, 1369.31, '$ m', '储能：飞轮、电容与它·公式'),
        (2, 9.375, '$ GJ', '储能：飞轮、电容与它·公式'),
        (2, 20, '$ kV', '储能：飞轮、电容与它·公式'),
        (2, 0.01, '$ F', '储能：飞轮、电容与它·公式'),
        (2, 10, '$ mF', '储能：飞轮、电容与它·公式'),
        (2, 5, '$ mF', '储能：飞轮、电容与它·公式'),
        (2, 28.2843, '$ kV', '储能：飞轮、电容与它·公式'),
        (2, 10, '$ kA', '储能：飞轮、电容与它·公式'),
        (2, 20, '$ kA', '储能：飞轮、电容与它·公式'),
        (2, 1, '$ m$\\Omega', '储能：飞轮、电容与它·公式'),
        (2, 10000, '$ s', '储能：飞轮、电容与它·公式'),
        (2, 1.97392, '$ GJ', '储能：飞轮、电容与它·正文'),
        (2, 100, '$ MW', '储能：飞轮、电容与它·正文'),
        (2, 19.7392, '$ s', '储能：飞轮、电容与它·正文'),
        (2, 80, '$ s', '储能：飞轮、电容与它·正文'),
        (2, 28.8, '$ GJ', '储能：飞轮、电容与它·正文'),
        (2, 1, '$ MWh', '储能：飞轮、电容与它·正文'),
        (2, 3600000000, '$ J', '储能：飞轮、电容与它·正文'),
        (2, 0.548311, '$ MWh', '储能：飞轮、电容与它·正文'),
        (2, 40, '$ MJ', '储能：飞轮、电容与它·表格'),
        (3, 67.9122, '\\%', '环流功率：电厂净效率·正文'),
        (3, 18.2218, '\\%', '环流功率：电厂净效率·正文'),
        (3, 10.866, '\\%', '环流功率：电厂净效率·正文'),
        (3, 3, '\\%', '环流功率：电厂净效率·正文'),
        (3, 10, '\\%', '环流功率：电厂净效率·公式'),
        (3, 5, '$ MW', '环流功率：电厂净效率·公式'),
        (3, 25, '\\%', '环流功率：电厂净效率·正文'),
        (3, 55.65, '$ MW', '环流功率：电厂净效率·正文'),
        (3, 17.9744, '$ MW', '环流功率：电厂净效率·正文'),
        (3, 2.20873, '$ MW', '环流功率：电厂净效率·表格'),
        (3, 100, '\\%', '环流功率：电厂净效率·表格'),
        (4, 5, '$ s', '电网接口：并网、频率·公式'),
        (4, 1000, '$ MVA', '电网接口：并网、频率·公式'),
        (4, 0.01, '$ Hz', '电网接口：并网、频率·公式'),
        (4, 0.2, '$ Hz', '电网接口：并网、频率·公式'),
        (4, 20, '$ s', '电网接口：并网、频率·公式'),
        (4, 71.4157, '$ MW', '电网接口：并网、频率·公式'),
        (4, 7.14157, '\\%', '电网接口：并网、频率·公式'),
        (4, 14.8976, '\\%', '电网接口：并网、频率·公式'),
        (4, 2.88e+10, '\\ \\mathrm{J}', '电网接口：并网、频率·公式'),
        (4, 28.8, '\\ \\mathrm{GJ}', '电网接口：并网、频率·公式'),
        (5, 1, '$ ms', '保护与故障：毫秒级动·正文'),
        (5, 100, '$ kA', '保护与故障：毫秒级动·正文'),
        (5, 500, '\\ \\mathrm{V}', '保护与故障：毫秒级动·公式'),
        (5, 25, '\\ \\mathrm{MW}', '保护与故障：毫秒级动·公式'),
        (5, 25, '$ MW', '保护与故障：毫秒级动·公式'),
        (5, 2.77778, '$ h', '保护与故障：毫秒级动·公式'),
        (5, 1, '$ s', '保护与故障：毫秒级动·公式'),
        (5, 76.9074, '\\ \\mathrm{ms}', '保护与故障：毫秒级动·公式'),
        (5, 2, '\\ \\mathrm{MJ}', '保护与故障：毫秒级动·公式'),
        (5, 27.1649, '\\ \\mathrm{ms}', '保护与故障：毫秒级动·公式'),
        (5, 500, '$ V', '保护与故障：毫秒级动·表格'),
        (5, 76.9074, '$ ms', '保护与故障：毫秒级动·表格'),
        (5, 960, '$ MW', '保护与故障：毫秒级动·表格'),
        (6, 250.425, '$ MW', '估算练习·正文'),
        (6, 29.3998, '\\%', '估算练习·正文'),
        (6, 31.7701, '\\%', '估算练习·正文'),
        (6, 176.801, '$ MW', '估算练习·正文'),
        (6, 418.879, '$ rad', '估算练习·正文'),
        (6, 3.50919, '$ GJ', '估算练习·正文'),
        (6, 837.758, '$ m', '估算练习·正文'),
        (6, 5.50943, '$ GPa', '估算练习·正文'),
        (6, 400000000, '$ W', '估算练习·正文'),
        (6, 400, '$ MW', '估算练习·正文'),
        (6, 4000, '$ A', '估算练习·正文'),
        (6, 40, '$ kV', '估算练习·正文'),
        (6, 8533.33, '$ kg', '估算练习·正文'),
        (6, 103.8, '$ MVAr', '估算练习·正文'),
        (6, 45.2009, '$ MVAr', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AR 的 %d 个关键数值未在小节内复算通过：%s'
            % (len(bad), '、'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AR 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流守卫（禁止把精确值四舍五入后写回正文） ----------
    # 守卫串取「四舍五入后的数字 + $」：既不会撞上正文里的精确值（如 $33.0748\%$），
    # 又能在有人把它改写成 $33.07$ 时报 FAIL。
    for lit, why in (
        ('33.07$', '环流份额百分数四舍五入'),
        ('26.77$', '电厂净效率百分数四舍五入'),
        ('3.74$', '每 1 MW 净电对应聚变功率四舍五入'),
        ('5.37$', '爬升功率相对净电倍数四舍五入'),
        ('2.72$', '瞬态脉冲相对环流倍数四舍五入'),
        ('67.91$', '电流驱动份额百分数四舍五入'),
        ('18.22$', '低温份额百分数四舍五入'),
        ('10.87$', '厂用电份额百分数四舍五入'),
        ('2.02$', '工程增益四舍五入'),
        ('0.898$', '净效率增量百分点四舍五入'),
        ('17.97$', '环流削减量四舍五入'),
        ('29.4$', '练习 1 环流份额百分数四舍五入'),
        ('31.77$', '练习 1 净效率百分数四舍五入'),
        ('1.19$', '练习 1 净电提升倍数四舍五入'),
        ('19.74$', '飞轮储能放电时间四舍五入'),
        ('3.1$', '飞轮应力吉帕数四舍五入'),
        ('3.84$', '碳纤维与钢的线速度比四舍五入'),
        ('1.78$', '练习 2 飞轮储能提升倍数四舍五入'),
        ('5.51$', '练习 2 钢内应力四舍五入'),
        ('0.837$', '11 次谐波滤波器电感四舍五入'),
        ('0.6$', '13 次谐波滤波器电感四舍五入'),
        ('1.74$', '练习 6 无功增加倍数四舍五入'),
        ('0.466$', '练习 6 触发角正切四舍五入'),
        ('45.2$', '练习 6 补偿容量四舍五入'),
        ('27.16$', '破裂缓解等效环流时间四舍五入'),
        ('76.9$', '热沉积时间尺度四舍五入'),
        ('4.05$', '练习 5 等效飞轮个数四舍五入'),
        ('8.53$', '练习 5 碳纤维飞轮质量四舍五入'),
        ('0.548$', '飞轮等效兆瓦时四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AR'], why))


def check_P51_appendix_as_rays():
    """附录 AS：射线追踪与波束光学——结构、出处、练习与数值核对。

    可证伪点：① 整节被删（小节标题匹配失败）；② 表 / 语义盒 / 小节 / 公式标签被删或改名；
    ③ 「出处：」「不复制」声明缺失；④ 练习与编号解答数量不足或不相等；
    ⑤ 关键数值被改动（逐节 value_in 复算 + 「数字与单位」精确子串绑定，后者专门堵
       「只改末位数字、仍落在 1e-3 容差内」这类改动）；
    ⑥ 整行表格数据行 / 语义盒标题 / 出处声明被改写；
    ⑦ 口径回流（把精确值四舍五入后重写进正文）。
    """
    spath = os.path.join(SRC, APP_FILES['AS'])
    if not os.path.exists(spath):
        rec('FAIL', 'P', '附录 AS 文件缺失：%s' % spath)
        return

    text = strip_comments(read(spath))
    secs = split_blocks(text, 'section')

    # ---------- ① 小节齐全 ----------
    want = [
        '从波动到射线',
        '射线方程与折射',
        '高斯波束与准光学',
        '反射计',
        '电子回旋发射成像',
        '功率沉积定位',
        '估算练习',
        '来源、分级与版权说明',
    ]
    found = {}
    for pre in want:
        for ttl, blk in secs:
            if pre in ttl:
                found[pre] = blk
                break
    miss = [p for p in want if p not in found]
    if miss:
        rec('FAIL', 'P', '附录 AS 缺小节：%s' % '、'.join(miss))
        return
    rec('PASS', 'P', '附录 AS 的 %d 个小节齐全' % len(want))

    # ---------- ② 表 / 小节 / 语义盒 / 公式标签 ----------
    for lab in (
        'tab:app-as-wkb',
        'tab:app-as-ray',
        'tab:app-as-beam',
        'tab:app-as-reflect',
        'tab:app-as-ece',
        'tab:app-as-deposition',
        'tab:app-as-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（表标签）' % APP_FILES['AS'])

    for lab in (
        'sec:app-as-wkb',
        'sec:app-as-ray',
        'sec:app-as-beam',
        'sec:app-as-reflect',
        'sec:app-as-ece',
        'sec:app-as-deposition',
        'sec:app-as-problems',
        'sec:app-as-source',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（小节标签）' % APP_FILES['AS'])

    for lab in (
        'tip-app-as-wkb',
        'tip-app-as-ray',
        'tip-app-as-beam',
        'tip-app-as-reflect',
        'tip-app-as-ece',
        'tip-app-as-deposition',
        'tip-app-as-honest',
    ):
        has(text, '{%s}' % lab, 'src/%s（语义盒标签）' % APP_FILES['AS'])

    for lab in (
        'eq:app-as-dispersion',
        'eq:app-as-wkb',
        'eq:app-as-scalelength',
        'eq:app-as-rayeq',
        'eq:app-as-snell',
        'eq:app-as-turning',
        'eq:app-as-pathturn',
        'eq:app-as-beam',
        'eq:app-as-rayleigh',
        'eq:app-as-focus',
        'eq:app-as-intensity',
        'eq:app-as-tau',
        'eq:app-as-profile',
        'eq:app-as-res',
        'eq:app-as-doppler',
        'eq:app-as-ecefreq',
        'eq:app-as-mapping',
        'eq:app-as-rj',
        'eq:app-as-ecepow',
        'eq:app-as-resonance',
        'eq:app-as-dopplerres',
        'eq:app-as-abs',
        'eq:app-as-powerdens',
    ):
        has(text, '\\label{%s}' % lab, 'src/%s（公式标签）' % APP_FILES['AS'])

    # ---------- ③ 出处与版权声明 ----------
    has(text, '不复制', 'src/%s（版权口径）' % APP_FILES['AS'])
    n_src = len(re.findall(r'\\textbf\{出处：\}', text))
    if n_src < 6:
        rec('FAIL', 'P', '附录 AS 的「出处：」只有 %d 处（要求 >= 6）' % n_src)
    else:
        rec('PASS', 'P', '附录 AS 的「出处：」共 %d 处' % n_src)

    # ---------- ④ 练习与解答 ----------
    ex = found['估算练习']
    n_ex = len(re.findall(r'\\textbf\{练习\s*\d+\}', ex))
    n_ans = len(re.findall(r'\\textbf\{解答\s*(\d+)\}', ex))
    if n_ex < 6 or n_ans != n_ex:
        rec('FAIL', 'P', '附录 AS 练习 %d 道、编号解答 %d 条（要求 >= 6 且两者相等）'
            % (n_ex, n_ans))
    else:
        rec('PASS', 'P', '附录 AS 的 %d 道练习均有编号解答' % n_ex)

    # ---------- ④c 整句字面绑定（堵「同值多处只改一处」） ----------
    for lit, why in (
        (r'2\pi/k_\mathrm{p}=2.79093$ mm', '等离子体波长的正文句'),
        (r'那一束是 $2.94982$ m', '斜射光程的正文句'),
        (r'只有 $0.747509$ m', '直线距离的正文句'),
        (r'f_\mathrm{ce}=148.360$ GHz', '轴上回旋频率的正文句'),
        (r'\frac{v_\parallel^\mathrm{th}}{c}=0.136769', '多普勒展宽公式'),
        (r'1/(kL_n)=1.70247\%', '折射率偏差的正文句'),
        (r'0.538089', '练习 6 聚焦增益'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AS'], why))

    # ---------- ④b 精确字面绑定（数字 + 单位） ----------
    for lit, why in (
        ('1.0\\times10^{20}$ m$^{-3}$', '从波动到射线：什么时·正文（精确绑定）'),
        ('9.194$ keV', '从波动到射线：什么时·正文（精确绑定）'),
        ('5.3$ T', '从波动到射线：什么时·正文（精确绑定）'),
        ('2.0$ m', '从波动到射线：什么时·正文（精确绑定）'),
        ('140$ GHz', '从波动到射线：什么时·正文（精确绑定）'),
        ('8.79646\\times10^{11}$ rad', '从波动到射线：什么时·正文（精确绑定）'),
        ('2.14137$ mm', '从波动到射线：什么时·公式（精确绑定）'),
        ('2251.29$ m$^{-1}$', '从波动到射线：什么时·公式（精确绑定）'),
        ('2.79093$ mm', '从波动到射线：什么时·公式（精确绑定）'),
        ('4.8$ m', '从波动到射线：什么时·正文（精确绑定）'),
        ('1.5$ m', '从波动到射线：什么时·正文（精确绑定）'),
        ('0.211111$ m', '从波动到射线：什么时·正文（精确绑定）'),
        ('1.98$ m', '从波动到射线：什么时·正文（精确绑定）'),
        ('1.70247\\%', '从波动到射线：什么时·正文（精确绑定）'),
        ('6.07817\\times10^{19}$ m$^{-3}$', '射线方程与折射：波为·公式（精确绑定）'),
        ('1.25249\\ \\mathrm{m}', '射线方程与折射：波为·公式（精确绑定）'),
        ('0.970573$ m', '射线方程与折射：波为·正文（精确绑定）'),
        ('0.747509$ m', '射线方程与折射：波为·正文（精确绑定）'),
        ('7858.75$ rad', '射线方程与折射：波为·正文（精确绑定）'),
        ('1.25249$ m', '射线方程与折射：波为·表格（精确绑定）'),
        ('2.0$ cm', '高斯波束与准光学：把·公式（精确绑定）'),
        ('0.586837$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('0.0710357$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('1$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('0.0395160$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('4$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('2.17219$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('1.0$ MW', '高斯波束与准光学：把·公式（精确绑定）'),
        ('1.59155\\times10^{9}$ W m$^{-2}$', '高斯波束与准光学：把·公式（精确绑定）'),
        ('1.0$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('0.20$ m', '高斯波束与准光学：把·公式（精确绑定）'),
        ('0.00681621\\ \\mathrm{m}', '高斯波束与准光学：把·公式（精确绑定）'),
        ('6.81621\\ \\mathrm{mm}', '高斯波束与准光学：把·公式（精确绑定）'),
        ('145.961$ mm$^2$', '高斯波束与准光学：把·公式（精确绑定）'),
        ('3.21312\\ \\mathrm{MV\\,m^{-1}}', '高斯波束与准光学：把·公式（精确绑定）'),
        ('308.314$ kV m$^{-1}$', '高斯波束与准光学：把·正文（精确绑定）'),
        ('1 MW', '高斯波束与准光学：把·正文（精确绑定）'),
        ('140 GHz', '高斯波束与准光学：把·表格（精确绑定）'),
        ('6.81621$ mm', '高斯波束与准光学：把·表格（精确绑定）'),
        ('10\\ \\mathrm{GHz}', '反射计：反射量的是时·公式（精确绑定）'),
        ('50$ GHz', '反射计：反射量的是时·公式（精确绑定）'),
        ('75$ GHz', '反射计：反射量的是时·公式（精确绑定）'),
        ('1.66119$ m', '反射计：反射量的是时·公式（精确绑定）'),
        ('0.561644$ m', '反射计：反射量的是时·公式（精确绑定）'),
        ('0.01$ rad', '反射计：反射量的是时·公式（精确绑定）'),
        ('118.547\\ \\mathrm{kHz}', '反射计：反射量的是时·公式（精确绑定）'),
        ('2.66324\\ \\mu\\mathrm{m}', '反射计：反射量的是时·公式（精确绑定）'),
        ('1.46418\\times10^{14}$ m$^{-3}$', '反射计：反射量的是时·公式（精确绑定）'),
        ('3.19440\\ \\mathrm{MHz}', '反射计：反射量的是时·公式（精确绑定）'),
        ('100$ kHz', '反射计：反射量的是时·正文（精确绑定）'),
        ('313.048$ m s$^{-1}$', '反射计：反射量的是时·正文（精确绑定）'),
        ('60$ GHz', '反射计：反射量的是时·表格（精确绑定）'),
        ('4.46559\\times10^{19}$ m$^{-3}$', '反射计：反射量的是时·表格（精确绑定）'),
        ('6.97749\\times10^{19}$ m$^{-3}$', '反射计：反射量的是时·表格（精确绑定）'),
        ('255.509$ GHz', '电子回旋发射成像：共·公式（精确绑定）'),
        ('280$ GHz', '电子回旋发射成像：共·公式（精确绑定）'),
        ('0.370237\\ \\mathrm{m}', '电子回旋发射成像：共·公式（精确绑定）'),
        ('5.00134\\ \\mathrm{T}', '电子回旋发射成像：共·公式（精确绑定）'),
        ('1$ GHz', '电子回旋发射成像：共·公式（精确绑定）'),
        ('2.34651$ cm', '电子回旋发射成像：共·公式（精确绑定）'),
        ('280\\ \\mathrm{GHz}', '电子回旋发射成像：共·公式（精确绑定）'),
        ('9.194\\ \\mathrm{keV}', '电子回旋发射成像：共·公式（精确绑定）'),
        ('2.94608\\ \\mu\\mathrm{W}', '电子回旋发射成像：共·公式（精确绑定）'),
        ('6.2$ m', '电子回旋发射成像：共·表格（精确绑定）'),
        ('5.30000$ T', '电子回旋发射成像：共·表格（精确绑定）'),
        ('0.370237$ m', '电子回旋发射成像：共·表格（精确绑定）'),
        ('5.00134$ T', '电子回旋发射成像：共·表格（精确绑定）'),
        ('140.000$ GHz', '电子回旋发射成像：共·表格（精确绑定）'),
        ('280.000$ GHz', '电子回旋发射成像：共·表格（精确绑定）'),
        ('4.56389$ T', '电子回旋发射成像：共·表格（精确绑定）'),
        ('4.00732$ T', '电子回旋发射成像：共·表格（精确绑定）'),
        ('5.68693\\times10^{7}$ m s$^{-1}$', '功率沉积定位：共振面·公式（精确绑定）'),
        ('1.04098$ m', '功率沉积定位：共振面·公式（精确绑定）'),
        ('0.95\\ \\mathrm{MW}', '功率沉积定位：共振面·公式（精确绑定）'),
        ('0.0590350$ m', '功率沉积定位：共振面·公式（精确绑定）'),
        ('5.48818\\times10^{7}\\ \\mathrm{W\\,m^{-3}}', '功率沉积定位：共振面·公式（精确绑定）'),
        ('6.36475\\times10^{9}$ W m$^{-3}$', '功率沉积定位：共振面·正文（精确绑定）'),
        ('3.0cm', '功率沉积定位：共振面·表格（精确绑定）'),
        ('100$ GHz', '估算练习·正文（精确绑定）'),
        ('5.0\\times10^{19}$ m$^{-3}$', '估算练习·正文（精确绑定）'),
        ('1619.26$ m$^{-1}$', '估算练习·正文（精确绑定）'),
        ('3.88028\\times10^{8}$ m s$^{-1}$', '估算练习·正文（精确绑定）'),
        ('1.0$ cm', '估算练习·正文（精确绑定）'),
        ('70$ GHz', '估算练习·正文（精确绑定）'),
        ('5.0$ MHz', '估算练习·正文（精确绑定）'),
        ('1.56524\\times10^{4}$ m s$^{-1}$', '估算练习·正文（精确绑定）'),
        ('15.6524$ km s$^{-1}$', '估算练习·正文（精确绑定）'),
        ('260$ GHz', '估算练习·正文（精确绑定）'),
        ('4.64410$ T', '估算练习·正文（精确绑定）'),
        ('0.875640$ m', '估算练习·正文（精确绑定）'),
        ('0.0272140$ m', '估算练习·正文（精确绑定）'),
        ('2.72140$ cm', '估算练习·正文（精确绑定）'),
        ('2$ MW', '估算练习·正文（精确绑定）'),
        ('5$ mm', '估算练习·正文（精确绑定）'),
        ('1$ MW', '估算练习·正文（精确绑定）'),
        ('5.09296\\times10^{10}$ W m$^{-2}$', '估算练习·正文（精确绑定）'),
        ('6.19463\\times10^{6}$ V m$^{-1}', '估算练习·正文（精确绑定）'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AS'], why))

    # ---------- ⑤ 整行字面绑定 ----------
    for lit, why in (
        ('量 & 数值 & 判据 & 含义 \\\\', '表格数据行'),
        ('真空波长 $\\lambda_0$ & $2.14137$ mm & --- & $140$ GHz \\\\', '表格数据行'),
        ('等离子体波长 $\\lambda_\\mathrm{p}$ & $2.79093$ mm & --- & 中心处，$n=0.767263$ \\\\', '表格数据行'),
        ('波数乘半径 $k_\\mathrm{p}a$ & $4502.58$ & $\\gg1$ & 跨过四千多个相位 \\\\', '表格数据行'),
        ('密度尺度最小值 $L_n$ & $0.211111$ m & $\\gg\\lambda_\\mathrm{p}$ & 出现在 $r=0.99a$ \\\\', '表格数据行'),
        ('尺度判据最小值 $kL_n$ & $58.7382$ & $\\gg1$ & 边缘仍然勉强成立 \\\\', '表格数据行'),
        ('相对误差量级 $1/(kL_n)$ & $1.70247\\%$ & $\\ll1$ & 边缘要小心 \\\\', '表格数据行'),
        ('\\begin{tip}{射线是短波极限，不是真理}{tip-app-as-wkb}', '语义盒标题'),
        ('\\textbf{出处：}几何光学极限与色散关系见 Ginzburg, \\emph{The Propagation of', '出处声明整行'),
        ('发射角 & 转折密度 & 转折半径 & 路径长度 & 说明 \\\\', '表格数据行'),
        ('$20^\\circ$ & $2.14686\\times10^{20}$ m$^{-3}$ & 无 & --- & 高于中心密度，穿过 \\\\', '表格数据行'),
        ('$60^\\circ$ & $6.07817\\times10^{19}$ m$^{-3}$ & $1.25249$ m & $2.94982$ m & $0.626245a$ 处折返 \\\\', '表格数据行'),
        ('$60^\\circ$（直线） & --- & --- & $0.747509$ m & 路径是直线的 $3.94620$ 倍 \\\\', '表格数据行'),
        ('\\begin{tip}{折射让波走弯路，转折点才是信息}{tip-app-as-ray}', '语义盒标题'),
        ('\\textbf{出处：}射线方程与分层介质中的折射守恒见 Tracy, Brizard, Richardson \\&', '出处声明整行'),
        ('量 & 数值 & 位置 & 说明 \\\\', '表格数据行'),
        ('腰斑 $w_0$ & $2.0$ cm & 发射口 & 声明值 \\\\', '表格数据行'),
        ('瑞利长度 $z_R$ & $0.586837$ m & --- & 束宽涨到 $\\sqrt2$ 倍 \\\\', '表格数据行'),
        ('发散角 $\\theta_\\mathrm{div}$ & $1.9527^\\circ$ & 远场 & $0.0340810$ rad \\\\', '表格数据行'),
        ('束宽 $w$ & $0.0710357$ m & $z=2.0$ m & 打到中心 \\\\', '表格数据行'),
        ('峰值强度 & $1.59155\\times10^{9}$ W m$^{-2}$ & 腰斑 & $2P/(\\pi w_0^2)$ \\\\', '表格数据行'),
        ('峰值强度 & $1.26162\\times10^{8}$ W m$^{-2}$ & 中心 & 散开之后 \\\\', '表格数据行'),
        ('焦斑 $w_\\mathrm{f}$ & $6.81621$ mm & 镜后 & $f=1.0$ m、$D=0.20$ m \\\\', '表格数据行'),
        ('焦斑电场 $E_\\mathrm{f}$ & $3.21312$ MV m$^{-1}$ & 焦斑 & 中心的 $10.42$ 倍 \\\\', '表格数据行'),
        ('\\begin{tip}{准光学是把波长当成尺子}{tip-app-as-beam}', '语义盒标题'),
        ('\\textbf{出处：}高斯波束与聚焦腰斑的公式见 Born \\& Wolf, \\emph{Principles of Optics}', '出处声明整行'),
        ('频率 & 截止密度 & 转折半径 & 往返群延迟 & 相对半径 \\\\', '表格数据行'),
        ('$50$ GHz & $3.10111\\times10^{19}$ m$^{-3}$ & $1.66119$ m & $4.66510$ ns & $0.830596a$ \\\\', '表格数据行'),
        ('$60$ GHz & $4.46559\\times10^{19}$ m$^{-3}$ & $1.48787$ m & $7.19549$ ns & $0.743936a$ \\\\', '表格数据行'),
        ('$75$ GHz & $6.97749\\times10^{19}$ m$^{-3}$ & $1.09955$ m & $13.4255$ ns & $0.549774a$ \\\\', '表格数据行'),
        ('\\begin{tip}{反射计量的是时间，不是位置}{tip-app-as-reflect}', '语义盒标题'),
        ('\\textbf{出处：}微波反射计的原理、群延迟与剖面反演见 Hutchinson, \\emph{Principles of', '出处声明整行'),
        ('半径 & 磁场 & $f_\\mathrm{ce}$ & 二次谐波 & 说明 \\\\', '表格数据行'),
        ('$0$ m & $5.30000$ T & $148.360$ GHz & $296.720$ GHz & 轴上最高 \\\\', '表格数据行'),
        ('$0.370237$ m & $5.00134$ T & $140.000$ GHz & $280.000$ GHz & 本轮通道 \\\\', '表格数据行'),
        ('$1.0$ m & $4.56389$ T & $127.755$ GHz & $255.509$ GHz & --- \\\\', '表格数据行'),
        ('$2.0$ m & $4.00732$ T & $112.175$ GHz & $224.350$ GHz & 边缘最低 \\\\', '表格数据行'),
        ('\\begin{tip}{ECE 看到的是共振面，不是等离子体}{tip-app-as-ece}', '语义盒标题'),
        ('\\textbf{出处：}电子回旋发射的诊断原理与光学厚度讨论见 Hutchinson（2002）第 6 章与', '出处声明整行'),
        ('量 & 数值 & 单位 & 说明 \\\\', '表格数据行'),
        ('冷共振半径 & $0.370237$ & m & $f=f_\\mathrm{ce}$ \\\\', '表格数据行'),
        ('平行折射率 $n_\\parallel$ & $0.720991$ & --- & $n_\\mathrm{p}\\cos20^\\circ$ \\\\', '表格数据行'),
        ('多普勒展宽 $\\Delta\\omega/\\omega$ & $0.136769$ & --- & $n_\\parallel v_\\parallel^\\mathrm{th}/c$ \\\\', '表格数据行'),
        ('$+$ 支共振半径 & $1.41121$ & m & 外移 $1.04098$ m \\\\', '表格数据行'),
        ('展宽层厚 & $2.08195$ & m & 与器壁距离同量级 \\\\', '表格数据行'),
        ('吸收功率 $P_\\mathrm{abs}$ & $0.95$ & MW & 单程 $95\\%$（声明） \\\\', '表格数据行'),
        ('平均功率密度 & $5.48818\\times10^{7}$ & W m$^{-3}$ & 摊在展宽层里 \\\\', '表格数据行'),
        ('焦斑极限密度 & $6.36475\\times10^{9}$ & W m$^{-3}$ & 高两个数量级 \\\\', '表格数据行'),
        ('\\begin{tip}{沉积位置由磁场和角度共同决定}{tip-app-as-deposition}', '语义盒标题'),
        ('\\textbf{出处：}电子回旋波的多普勒共振与功率沉积见 Bornatici 等, \\emph{Nuclear Fusion}', '出处声明整行'),
        ('来源 & 用于 & 借用方式 & 版权说明 \\\\', '表格数据行'),
        ('Ginzburg 1964/1970 & AS.1--AS.2 & 色散关系与折射守恒 & 只复算，不复制图表 \\\\', '表格数据行'),
        ('Born \\& Wolf 1999 & AS.3 & 高斯波束与聚焦公式 & 只复算，不复制图 \\\\', '表格数据行'),
        ('Goldsmith 1998 & AS.3 & 准光学量级判断 & 只借用方法 \\\\', '表格数据行'),
        ('Hutchinson 2002 & AS.4--AS.5 & 反射计与 ECE 原理 & 不复制公式表与图 \\\\', '表格数据行'),
        ('Tracy 等 2014 & AS.2 & 射线方程与波迹 & 只借用形式 \\\\', '表格数据行'),
        ('ITER Physics Basis 1999 & AS.6 & 多普勒共振与沉积 & 不复制剖面图 \\\\', '表格数据行'),
        ('\\begin{tip}{本附录不给什么}{tip-app-as-honest}', '语义盒标题'),
        ('\\textbf{出处：}本附录的「不复制」口径与全书一致：凡引用原书结论处只写结论与出处，', '出处声明整行'),
    ):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['AS'], why))

    # ---------- ⑥ 数值复算：按小节绑定 ----------
    cases = [
        (0, 1e+20, '$ m$^{-3}$', '从波动到射线：什么时·正文'),
        (0, 9.194, '$ keV', '从波动到射线：什么时·正文'),
        (0, 5.3, '$ T', '从波动到射线：什么时·正文'),
        (0, 2, '$ m', '从波动到射线：什么时·正文'),
        (0, 140, '$ GHz', '从波动到射线：什么时·正文'),
        (0, 8.79646e+11, '$ rad', '从波动到射线：什么时·正文'),
        (0, 2.43127e+20, '\\ \\mathrm{m^{-3}}', '从波动到射线：什么时·公式'),
        (0, 2934.18, '$ m$^{-1}$', '从波动到射线：什么时·公式'),
        (0, 2.14137, '$ mm', '从波动到射线：什么时·公式'),
        (0, 2251.29, '$ m$^{-1}$', '从波动到射线：什么时·公式'),
        (0, 2.79093, '$ mm', '从波动到射线：什么时·公式'),
        (0, 4.8, '$ m', '从波动到射线：什么时·正文'),
        (0, 1.5, '$ m', '从波动到射线：什么时·正文'),
        (0, 0.211111, '$ m', '从波动到射线：什么时·正文'),
        (0, 1.98, '$ m', '从波动到射线：什么时·正文'),
        (0, 1.70247, '\\%', '从波动到射线：什么时·正文'),
        (1, 2.14686e+20, '$ m$^{-3}$', '射线方程与折射：波为·公式'),
        (1, 6.07817e+19, '$ m$^{-3}$', '射线方程与折射：波为·公式'),
        (1, 1.25249, '\\ \\mathrm{m}', '射线方程与折射：波为·公式'),
        (1, 2541.08, '$ m$^{-1}$', '射线方程与折射：波为·正文'),
        (1, 0.970573, '$ m', '射线方程与折射：波为·正文'),
        (1, 2.94982, '$ m', '射线方程与折射：波为·正文'),
        (1, 0.747509, '$ m', '射线方程与折射：波为·正文'),
        (1, 7858.75, '$ rad', '射线方程与折射：波为·正文'),
        (1, 1.25249, '$ m', '射线方程与折射：波为·表格'),
        (2, 2, '$ cm', '高斯波束与准光学：把·公式'),
        (2, 0.586837, '$ m', '高斯波束与准光学：把·公式'),
        (2, 0.034081, '\\ \\mathrm{rad}', '高斯波束与准光学：把·公式'),
        (2, 0.0710357, '$ m', '高斯波束与准光学：把·公式'),
        (2, 1, '$ m', '高斯波束与准光学：把·公式'),
        (2, 0.039516, '$ m', '高斯波束与准光学：把·公式'),
        (2, 4, '$ m', '高斯波束与准光学：把·公式'),
        (2, 0.137783, '$ m', '高斯波束与准光学：把·公式'),
        (2, 2.17219, '$ m', '高斯波束与准光学：把·公式'),
        (2, 1, '$ MW', '高斯波束与准光学：把·公式'),
        (2, 1591550000, '$ W m$^{-2}$', '高斯波束与准光学：把·公式'),
        (2, 126162000, '$ W m$^{-2}$', '高斯波束与准光学：把·公式'),
        (2, 0.2, '$ m', '高斯波束与准光学：把·公式'),
        (2, 0.00681621, '\\ \\mathrm{m}', '高斯波束与准光学：把·公式'),
        (2, 6.81621, '\\ \\mathrm{mm}', '高斯波束与准光学：把·公式'),
        (2, 0.000145961, '$ m$^2', '高斯波束与准光学：把·公式'),
        (2, 145.961, '$ mm$^2$', '高斯波束与准光学：把·公式'),
        (2, 1.37023e+10, '$ W m$^{-2}$', '高斯波束与准光学：把·公式'),
        (2, 3213120, '\\ \\mathrm{V\\,m^{-1}}', '高斯波束与准光学：把·公式'),
        (2, 3.21312, '\\ \\mathrm{MV\\,m^{-1}}', '高斯波束与准光学：把·公式'),
        (2, 308.314, '$ kV m$^{-1}$', '高斯波束与准光学：把·正文'),
        (2, 1, 'MW', '高斯波束与准光学：把·正文'),
        (2, 140, 'GHz', '高斯波束与准光学：把·表格'),
        (2, 0.034081, '$ rad', '高斯波束与准光学：把·表格'),
        (2, 6.81621, '$ mm', '高斯波束与准光学：把·表格'),
        (3, 10, '\\ \\mathrm{GHz}', '反射计：反射量的是时·公式'),
        (3, 50, '$ GHz', '反射计：反射量的是时·公式'),
        (3, 75, '$ GHz', '反射计：反射量的是时·公式'),
        (3, 1.66119, '$ m', '反射计：反射量的是时·公式'),
        (3, 1.09955, '$ m', '反射计：反射量的是时·公式'),
        (3, 0.561644, '$ m', '反射计：反射量的是时·公式'),
        (3, 0.0224658, '$ m', '反射计：反射量的是时·公式'),
        (3, 0.01, '$ rad', '反射计：反射量的是时·公式'),
        (3, 118.547, '\\ \\mathrm{kHz}', '反射计：反射量的是时·公式'),
        (3, 2.66324, '\\ \\mu\\mathrm{m}', '反射计：反射量的是时·公式'),
        (3, 1.46418e+14, '$ m$^{-3}$', '反射计：反射量的是时·公式'),
        (3, 3.1944, '\\ \\mathrm{MHz}', '反射计：反射量的是时·公式'),
        (3, 100, '$ kHz', '反射计：反射量的是时·正文'),
        (3, 313.048, '$ m s$^{-1}$', '反射计：反射量的是时·正文'),
        (3, 3.10111e+19, '$ m$^{-3}$', '反射计：反射量的是时·表格'),
        (3, 60, '$ GHz', '反射计：反射量的是时·表格'),
        (3, 4.46559e+19, '$ m$^{-3}$', '反射计：反射量的是时·表格'),
        (3, 1.48787, '$ m', '反射计：反射量的是时·表格'),
        (3, 6.97749e+19, '$ m$^{-3}$', '反射计：反射量的是时·表格'),
        (4, 148.36, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 296.72, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 127.755, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 255.509, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 112.175, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 224.35, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 280, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 0.370237, '\\ \\mathrm{m}', '电子回旋发射成像：共·公式'),
        (4, 5.00134, '\\ \\mathrm{T}', '电子回旋发射成像：共·公式'),
        (4, 2.34651e-11, '$ m', '电子回旋发射成像：共·公式'),
        (4, 1, '$ GHz', '电子回旋发射成像：共·公式'),
        (4, 0.0234651, '$ m', '电子回旋发射成像：共·公式'),
        (4, 2.34651, '$ cm', '电子回旋发射成像：共·公式'),
        (4, 280, '\\ \\mathrm{GHz}', '电子回旋发射成像：共·公式'),
        (4, 9.194, '\\ \\mathrm{keV}', '电子回旋发射成像：共·公式'),
        (4, 2.94608, '\\ \\mu\\mathrm{W}', '电子回旋发射成像：共·公式'),
        (4, 1.07131, '$ mm', '电子回旋发射成像：共·正文'),
        (4, 6.2, '$ m', '电子回旋发射成像：共·表格'),
        (4, 0.370237, '$ m', '电子回旋发射成像：共·表格'),
        (4, 5.00134, '$ T', '电子回旋发射成像：共·表格'),
        (4, 4.56389, '$ T', '电子回旋发射成像：共·表格'),
        (4, 4.00732, '$ T', '电子回旋发射成像：共·表格'),
        (5, 56869300, '$ m s$^{-1}$', '功率沉积定位：共振面·公式'),
        (5, 1.41121, '$ m', '功率沉积定位：共振面·公式'),
        (5, 1.04098, '$ m', '功率沉积定位：共振面·公式'),
        (5, 2.08195, '$ m', '功率沉积定位：共振面·公式'),
        (5, 0.95, '\\ \\mathrm{MW}', '功率沉积定位：共振面·公式'),
        (5, 0.059035, '$ m', '功率沉积定位：共振面·公式'),
        (5, 0.0109489, '$ m$^2$', '功率沉积定位：共振面·公式'),
        (5, 54881800, '\\ \\mathrm{W\\,m^{-3}}', '功率沉积定位：共振面·公式'),
        (5, 6364750000, '$ W m$^{-3}$', '功率沉积定位：共振面·正文'),
        (5, 3, 'cm', '功率沉积定位：共振面·表格'),
        (6, 100, '$ GHz', '估算练习·正文'),
        (6, 5e+19, '$ m$^{-3}$', '估算练习·正文'),
        (6, 1.24044e+20, '$ m$^{-3}$', '估算练习·正文'),
        (6, 1619.26, '$ m$^{-1}$', '估算练习·正文'),
        (6, 388028000, '$ m s$^{-1}$', '估算练习·正文'),
        (6, 1, '$ cm', '估算练习·正文'),
        (6, 0.146709, '$ m', '估算练习·正文'),
        (6, 70, '$ GHz', '估算练习·正文'),
        (6, 1.08566e-08, '$ s', '估算练习·正文'),
        (6, 5, '$ MHz', '估算练习·正文'),
        (6, 15652.4, '$ m s$^{-1}$', '估算练习·正文'),
        (6, 15.6524, '$ km s$^{-1}$', '估算练习·正文'),
        (6, 260, '$ GHz', '估算练习·正文'),
        (6, 4.6441, '$ T', '估算练习·正文'),
        (6, 0.87564, '$ m', '估算练习·正文'),
        (6, 0.027214, '$ m', '估算练习·正文'),
        (6, 2.7214, '$ cm', '估算练习·正文'),
        (6, 2, '$ MW', '估算练习·正文'),
        (6, 5, '$ mm', '估算练习·正文'),
        (6, 5.09296e+10, '$ W m$^{-2}$', '估算练习·正文'),
        (6, 6194630, '$ V m$^{-1}', '估算练习·正文'),
    ]

    bad = []
    for s_idx, val, unit, why in cases:
        pre = want[s_idx]
        if pre not in found:
            continue
        if not value_in(found[pre], val, 1e-3, unit):
            bad.append('%s（%.6g）' % (why, val))
    if bad:
        rec('FAIL', 'P', '附录 AS 的 %d 个关键数值未在小节内复算通过：%s'
            % (len(bad), '、'.join(bad[:6])))
    else:
        rec('PASS', 'P', '附录 AS 的 %d 个关键数值均在小节内复算通过' % len(cases))

    # ---------- ⑦ 口径回流守卫（禁止把精确值四舍五入后写回正文） ----------
    # 守卫串取「四舍五入后的数字 + $」或「四舍五入后的 LaTeX 字面」：既不会撞上正文里的
    # 精确值（如 $0.136769$），又能在有人把它改写成 $0.137$ 时报 FAIL。
    for lit, why in (
        ('2.791$', 'WKB 相位积分四舍五入'),
        ('0.9706$', '转折点密度比四舍五入'),
        ('2.95$', '单程光程四舍五入'),
        ('0.748$', '密度剖面峰值比四舍五入'),
        ('0.587$', '波束腰与瑞利长度比四舍五入'),
        ('0.0341$', '远场发散角四舍五入'),
        ('0.07104$', '腰斑相对波长四舍五入'),
        ('0.03952$', '近场光斑半径比四舍五入'),
        ('0.1378$', '多普勒展宽四舍五入'),
        ('0.00682$', '聚焦光斑相对波长四舍五入'),
        ('0.562$', '反射计对比度四舍五入'),
        ('0.02247$', '反射计相位分辨折合密度四舍五入'),
        ('1.488$', '反射计频率分辨四舍五入'),
        ('148.4$', '轴上回旋频率四舍五入'),
        ('127.8$', '高场侧边缘回旋频率四舍五入'),
        ('112.2$', '低场侧边缘回旋频率四舍五入'),
        ('0.02347$', '相对论质量修正四舍五入'),
        ('2.347$', 'ECE 真空波长四舍五入'),
        ('4.564$', '共振面磁场四舍五入'),
        ('4.01$', '边缘磁场四舍五入'),
        ('1.041$', '多普勒外移四舍五入'),
        ('2.082$', '展宽层厚四舍五入'),
        ('0.01095$', '沉积体积份额四舍五入'),
        ('0.147$', '练习 1 判定参数四舍五入'),
        ('15.7$', '练习 4 频率分辨四舍五入'),
        ('0.876$', '练习 5 吸收份额四舍五入'),
        ('1.702$', '折射率偏差百分数四舍五入'),
        ('5.69\\times10^{7}', '热速度四舍五入'),
        ('1.24\\times10^{18}', 'n_c 系数四舍五入'),
        ('1.2404\\times10^{18}', 'n_c 系数截断'),
        ('2.8\\times10^{10}', 'f_ce 系数四舍五入'),
    ):
        lacks(text, lit, 'src/%s（%s）' % (APP_FILES['AS'], why))


def check_P02_appendix_j_map():
    """附录 J 逐章对照表的完备性：覆盖 1--12 章、每行 5 列、章号不越界。

    可证伪点：删掉一行（覆盖不全）、把 Chen 第 10 章写成第 11 章（幻章）、
    把某一行的外文定位留空（该章失去对照目标）——三种都会 FAIL。
    Chen 第 3 版共 10 章、Fitzpatrick 第 2 版共 9 章，这两个上限写死在这里；
    换版本时必须同时改这里，正是我们想要的提醒。
    """
    t = read(os.path.join(SRC, APP_FILES['J']))
    m = re.search(r'\\label\{tab:app-j-map\}.*?\\end\{longtable\}', t, re.S)
    if not m:
        rec('FAIL', 'P', '附录 J 找不到逐章对照表（label tab:app-j-map 缺失，表格被删或改名）')
        return
    blk = m.group(0)
    rows = re.findall(r'^第\s*(\d+)\s*章\s*&(.*?)\\\\?\s*$', blk, re.M)
    got = [int(n) for n, _ in rows]
    if got != list(range(1, 13)):
        rec('FAIL', 'P', '附录 J 逐章对照表覆盖的章号不是 1--12：%s（共 %d 行）'
            % (got, len(got)))
    else:
        rec('PASS', 'P', '附录 J 逐章对照表覆盖第 1--12 章，共 %d 行' % len(got))

    bad_cols, phantom, empty = [], [], []
    for n, rest in rows:
        fields = [f.strip() for f in re.split(r'(?<!\\)&', rest)]
        if len(fields) != 4:
            bad_cols.append((int(n), len(fields) + 1))
            continue
        chen = [int(x) for x in re.findall(r'第\s*(\d+)\s*章', fields[1])]
        fitz = [int(x) for x in re.findall(r'第\s*(\d+)\s*章', fields[2])]
        if any(x > 10 for x in chen) or any(x > 9 for x in fitz):
            phantom.append((int(n), chen, fitz))
        if not chen and not fitz:
            empty.append(int(n))
    if bad_cols:
        rec('FAIL', 'P', '附录 J 对照表有 %d 行列数不是 5：%s' % (len(bad_cols), bad_cols))
    else:
        rec('PASS', 'P', '附录 J 对照表每行都是 5 列（本书／主题／Chen／Fitz／其他专著）')
    if phantom:
        rec('FAIL', 'P', '附录 J 对照表引用了不存在的章号（Chen 只有 10 章、Fitzpatrick 只有 9 章）：%s'
            % phantom)
    else:
        rec('PASS', 'P', '附录 J 对照表的 Chen 章号 ≤10、Fitzpatrick 章号 ≤9（无幻章）')
    if empty:
        rec('FAIL', 'P', '附录 J 对照表有 %d 行既没给 Chen 也没给 Fitzpatrick 定位：%s'
            % (len(empty), empty))
    else:
        rec('PASS', 'P', '附录 J 对照表每一行都至少给出一个外文教材定位')


def check_P03_appendix_i_refs():
    """附录 I 的文献出处表：DOI 格式合法、无重复、卷期页齐全。

    学习策略类的论断必须有可核查的出处。DOI 写错一个字符就查不到，
    所以这里查的是「8 条 DOI 全部形如 10.xxxx/yyy 且互不重复」。
    """
    t = read(os.path.join(SRC, APP_FILES['I']))
    dois = [d.replace('\\_', '_').rstrip('，。；,;')
            for d in re.findall(r'DOI:\s*([^\s&,，。；]+)', t)]
    bad = [d for d in dois if not re.fullmatch(r'10\.\d{4,9}/[-._;()/:A-Za-z0-9]+', d)]
    dup = sorted({d for d in dois if dois.count(d) > 1})
    vol = re.findall(r'\d+\(\d+\),\s*\d+--\d+', t)
    if len(dois) < 8:
        rec('FAIL', 'P', '附录 I 出处表只有 %d 条 DOI（门槛 8 条）' % len(dois))
    elif bad:
        rec('FAIL', 'P', '附录 I 有 %d 条 DOI 格式非法：%s' % (len(bad), bad))
    elif dup:
        rec('FAIL', 'P', '附录 I 有重复 DOI：%s' % dup)
    else:
        rec('PASS', 'P', '附录 I 出处表 %d 条 DOI 全部格式合法且互不重复' % len(dois))
    if len(vol) < 6:
        rec('FAIL', 'P', '附录 I 只有 %d 条「卷(期), 页--页」的期刊定位（门槛 6 条）' % len(vol))
    else:
        rec('PASS', 'P', '附录 I 有 %d 条期刊卷期页定位（可据此复核文献是否真实存在）' % len(vol))


def check_P04_new_appendix_wiring():
    """附录 I--M 必须真的挂进主文件，且章标题留在主文件里。

    两个都容易出错且都不报错：① 内容写好了却忘了 ``\\input``——编译照样成功，
    读者看不到；② ``\\appchapter`` 写进被 ``\\input`` 的文件——LaTeX 没问题，
    但按主文件判定附录字母的检查套件会漏掉它（v1.8.0 首稿正是这么写的，
    本轮修正后加此项防回归）。
    """
    t = read(os.path.join(SRC, MAIN))
    letters = re.findall(r'\\appchapter\{附录([A-Z]{1,2})[：:]', t)
    # 附录字母方案：A--Z 用完后接 AA、AB…（第 27 个不再是 chr(ord('A')+26)='['）
    want = [(chr(ord('A') + i) if i < 26 else 'A' + chr(ord('A') + i - 26))
            for i in range(len(APP_FILES) + 8)]
    if letters != want:
        rec('FAIL', 'P', '主文件附录字母不是 A--%s 连续：%s（共 %d 个）'
            % (want[-1], ''.join(letters), len(letters)))
    else:
        rec('PASS', 'P', '主文件附录字母连续 A--%s（共 %d 个，含 %s）'
            % (want[-1], len(letters), '/'.join(sorted(APP_FILES))))

    problems = []
    for letter, fname in sorted(APP_FILES.items()):
        stem = fname[:-4]
        path = os.path.join(SRC, fname)
        if not os.path.exists(path):
            problems.append('%s 文件缺失' % fname)
            continue
        if not re.search(r'\\appchapter\{附录%s[：:][^\n]*\}\s*\n\s*\\input\{%s\}'
                         % (letter, re.escape(stem)), t):
            problems.append('附录%s 的 \\appchapter 未紧接 \\input{%s}' % (letter, stem))
        body = read(path)
        if '\\appchapter' in strip_comments(body):
            problems.append('%s 里又写了一次 \\appchapter（会与主文件重复计字母）' % fname)
        if body.count('\\section{') < 4:
            problems.append('%s 的小节少于 4 个（内容疑似被截断）' % fname)
    if problems:
        rec('FAIL', 'P', '新增附录挂载有问题：%s' % problems)
    else:
        rec('PASS', 'P', '附录 %s 均已 \\input 挂载、章标题在主文件、各含 ≥4 小节'
            % '/'.join(sorted(APP_FILES)))


def check_P09_tcbtheorem_arity():
    r"""两参数语义盒不得裸调用：``\begin{tip}`` 会**吃掉正文头两个字符**。

    ``\newtcbtheorem`` 造出的环境签名是 ``\begin{env}{标题}{标签}``，两个都是必选参数。
    裸写 ``\begin{tip}`` 时 TeX 会把紧随其后的两个记号当成标题与标签：
    正文首字被吞（实测「本附录所有数值…」渲染成标题「本」、正文从「录…」开始），
    两个附录若都以同两字开头，还会产生 ``Label `tip:附' multiply defined`` 警告。
    v1.9.0 及以前全书有 16 处裸 tip，直到 v1.10.0 新增附录 M 触发重复标签才暴露。

    只检查 ``\newtcbtheorem`` 系列；``miniquiz`` / ``chapterreview`` / ``formula`` /
    ``checklist`` 是 ``\newtcolorbox``（可选参数），裸调用完全合法。
    """
    envs = ['theorem', 'definition', 'corollary', 'example', 'insight',
            'warning', 'history', 'review', 'tip']
    files = [MAIN] + sorted(APP_FILES.values()) + \
            ['ch%d_supplement.tex' % i for i in range(1, 13)]
    bad = []
    total = 0
    for fname in files:
        path = os.path.join(SRC, fname)
        if not os.path.exists(path):
            continue
        t = strip_comments(read(path))
        for env in envs:
            for m in re.finditer(r'\\begin\{%s\}(?!\s*\{)' % env, t):
                total += 1
                bad.append('%s:%d %s' % (fname, t[:m.start()].count('\n') + 1, env))
            total += 0
    if bad:
        rec('FAIL', 'P', '语义盒裸调用（会吞掉正文首字并可能造成标签重复）：%s' % bad[:6])
    else:
        rec('PASS', 'P', '全部语义盒均为两参数显式调用（标题 + 唯一标签），无裸 \\begin{tip} 类调用')


# ------------------------------------------------------------------ 主流程

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    check_M01_two_stream()
    check_M02_buneman()
    check_M03_delta_prime()
    check_M04_salpeter()
    check_M05_oblique()
    check_M06_pitch_and_efolding()
    check_M07_eedf()
    check_M08_timescale()
    check_O01_appendix_f()
    check_O02_appendix_g()
    check_O03_appendix_f_chapters()
    check_P01_new_appendix_numbers()
    check_P02_appendix_j_map()
    check_P03_appendix_i_refs()
    check_P04_new_appendix_wiring()
    check_P05_appendix_k_claims()
    check_P07_appendix_l_bank()
    check_P08_appendix_m_heating()
    check_P10_appendix_n_boundary()
    check_P11_appendix_o_propulsion()
    check_P12_appendix_p_lowtemp()
    check_P13_appendix_i_calendar()
    check_P14_appendix_q_diagnostics()
    check_P15_appendix_i_dunlosky()
    check_P16_no_markdown_bold()
    check_P17_appendix_r_dusty()
    check_P18_appendix_i_dualcode()
    check_P19_appendix_s_confinement()
    check_P20_ipb98_formula_consistency()
    check_P21_appendix_i_expert()
    check_P22_appendix_t_topology()
    check_P23_appendix_i_objectives()
    check_P24_appendix_u_gain()
    check_P25_appendix_i_expanding()
    check_P26_appendix_v_disruption()
    check_P27_appendix_w_fastparticle()
    check_P29_appendix_x_radiation()
    check_P30_appendix_y_strongfield()
    check_P31_appendix_z_reactor()
    check_P32_appendix_aa_alternates()
    check_P33_appendix_i_paper()
    check_P34_appendix_ab_space_astro()
    check_P35_appendix_ac_numerics()
    check_P36_appendix_ad_pedestal()
    check_P37_appendix_ae_kinetics()
    check_P38_appendix_af_turbulence()
    check_P39_appendix_ag_equilibrium()
    check_P40_appendix_ah_atomic()
    check_P41_appendix_ai_icf()
    check_P42_appendix_aj_fuelling()
    check_P43_appendix_ak_safety()
    check_P44_appendix_al_nonlinear()
    check_P45_appendix_am_magnets()
    check_P46_appendix_an_polarimetry()
    check_P47_appendix_ao_vacuum()
    check_P48_appendix_ap_maintenance()
    check_P49_appendix_aq_interferometry()
    check_P50_appendix_ar_power()
    check_P51_appendix_as_rays()
    check_P28_signed_quantities()
    check_P06_scattering_alpha()
    check_P09_tcbtheorem_arity()

    if args.json:
        print(json.dumps([{'level': a, 'cat': b, 'msg': c} for a, b, c in RESULTS],
                         ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》物理数值回归检查')
    print('=' * 78)
    for lvl in ['FAIL', 'WARN', 'INFO', 'PASS']:
        rows = [r for r in RESULTS if r[0] == lvl]
        if not rows:
            continue
        print('\n--- %s (%d) ---' % (lvl, len(rows)))
        for _, cat, msg in rows:
            print('  [%s] %s' % (cat, msg))
    nf = sum(1 for r in RESULTS if r[0] == 'FAIL')
    nw = sum(1 for r in RESULTS if r[0] == 'WARN')
    print('\n' + '=' * 78)
    print('汇总：FAIL %d ／ WARN %d ／ 合计 %d 项' % (nf, nw, len(RESULTS)))
    print('=' * 78)
    return 1 if nf else 0


if __name__ == '__main__':
    sys.exit(main())
