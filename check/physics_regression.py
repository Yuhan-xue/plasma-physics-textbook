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

    main = strip_comments(read(os.path.join(SRC, MAIN)))
    has(main, r'F(\eps) = \frac{2\sqrt{2 m_e \eps}}{A e^3}', 'src/%s（Druyvesteyn 公式）' % MAIN)
    has(main, r'\int_0^\infty F(\eps)\diff\eps=n_e', 'src/%s（EEDF 归一化）' % MAIN)
    has(main, '没有"固定偏差"这回事', 'src/%s（EEDF 警告框）' % MAIN)
    has(main, r'\diff v/\diff\eps=1/\sqrt{2m_e\eps}', 'src/%s（EEDF 雅可比）' % MAIN)


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
