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
}

# CODATA 2018 常数（与附录 D 的单位换算表同一组值）
E_CHG = 1.602176634e-19
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
M_P = 1.67262192369e-27
M_D = 2.013553212745 * M_P        # 氘核质量（含其电子）
M_AR = 39.948 * M_P               # 氩原子质量（≈ Ar+）


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
