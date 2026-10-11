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
            (r'w=0.0698\ \mathrm{m}$、$\tau_E=2.4762\ \mathrm{s}', '练习 4 给定的 w 与 τ_E')):
        has(text, lit, 'src/%s（%s）' % (APP_FILES['T'], why))

    # 口径回流防护（这些写法一旦出现即为 FAIL）
    for lit, why in (
            (r'\chi_\parallel=\frac{3.16\,T_e\,\tau_e}{m_i}', '平行电子热导被写成离子质量'),
            (r'\tau_e=\frac{3.44\times10^{5}\,T_e^{1/2}', '电子碰撞时间丢了 3/2 次幂'),
            (r'w=2\sqrt{\frac{\tilde\Psi}{B_\theta\,s}}', '岛宽系数从 4 退回 2'),
            (r'q(\rho)=q_0+(q_{95}-q_0)\,\rho,', 'q 剖面从 ρ² 退回 ρ'),
            (r'r_{\mathrm{inv}}=a\frac{1-q_0}{q_{95}-q_0}', '反转半径丢掉根号'),
            ("\\Psi_h''(r_{mn})=B_\\theta\\,s^2", '螺旋磁通二阶导数写成 s²')):
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
    letters = re.findall(r'\\appchapter\{附录([A-Z])[：:]', t)
    want = [chr(ord('A') + i) for i in range(len(APP_FILES) + 8)]
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
