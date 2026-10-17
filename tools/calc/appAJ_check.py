#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""附录 AJ 数值交叉核对：正文里每个数字都必须能在 tmp/appAJ_calc.txt（canon）或 DECLARED 里找到。
只读。用法：python tmp/appAJ_check.py"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, 'src', 'appendix_aj_fuelling.tex')
CANON = os.path.join(ROOT, 'tmp', 'appAJ_calc.txt')

# 声明级 / 文献级 / 单位换算 / 正文四舍五入（不参与 canon 比对）
DECLARED = {
    # 声明输入（脚本 [DECLARED] 段）
    1.0e20, 9194.0, 5.3, 2.0, 6.2, 832.2, 556.5, 111.3, 17.6, 3.7, 2.5, 0.95,
    0.1, 0.6, 0.05, 4.0e-4, 0.01, 100.0, 300.0, 12.32, 1.5, 255.0,
    5.00827e-27, 4.15135e-27,
    # 单位换算与正文四舍五入（canon 用 SI 打印，正文写 ms / 整数百分比）
    0.4, 84.0, 4.4, 0.646087, 1.26992, 12.0, 0.9, 0.81, 0.03, 5.0, 10.0,
    468.57, 2.0, 1000.0, 2500.0, 3.0, 4.0, 6.0,
    # 正文里的整数百分比 / 中间乘积 / 物理常数
    19.0, 98.0, 60.0, 15.0, 95.0, 16.8, 6.18573, 0.912256, 1.5e-3, 81.0,
    5.0e18, 3.0e18, 5.565e8, 1.10844e7, 9.88481e-7,
    1.602176634e-13, 1.380649e-23, 1.602176634e-19,
    # 量级与结构性数字
    0.0, 1.0, 10.0, 20.0, 40.0, 200.0, 1000.0, 0.2, 0.5, 0.75, 0.19,
    # 文献年份 / 卷期页
    1978, 1735, 21, 1995, 35, 657, 2007, 47, 203, 1999, 39, 2137, 1992, 32,
    611, 2011, 2000, 2002, 300, 91, 1988, 28, 2199, 4, 13, 14, 1, 2, 3, 6,
}


def nums_canon():
    t = io.open(CANON, encoding='utf-8', errors='replace').read()
    out = set()
    for m in re.finditer(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', t):
        try:
            out.add(float(m.group(0)))
        except ValueError:
            pass
    return out


def main():
    raw = io.open(TEX, encoding='utf-8', newline='').read()
    if '\r' in raw:
        print('[ABORT] 不是 LF')
        return 1
    txt = '\n'.join(l for l in raw.split('\n') if not l.lstrip().startswith('%'))
    cut = txt.index('\\section{来源、分级与版权说明}')
    body = txt[:cut]
    body = re.sub(r'p\{[^}]*\}', ' ', body)
    body = re.sub(r'\\label\{[^}]*\}', ' ', body)
    body = re.sub(r'\\(?:eqref|ref)\{[^}]*\}', ' ', body)
    sci = [(float(a) * 10.0 ** int(b), m.group(0)) for m, a, b in
           ((m, m.group(1), m.group(2)) for m in
            re.finditer(r'([0-9]+(?:\.[0-9]+)?)\s*\\times\s*10\^\{?(-?\d+)\}?', body))]
    body_nosci = re.sub(r'[0-9]+(?:\.[0-9]+)?\s*\\times\s*10\^\{?-?\d+\}?', ' ', body)
    plain = [(float(m.group(0)), m.group(0)) for m in
             re.finditer(r'(?<![0-9.eE])([0-9]+(?:\.[0-9]+)?)(?![0-9])', body_nosci)]
    canon = nums_canon()
    print('canon 数字 %d 个；正文科学计数法 %d 处；裸数字 %d 个'
          % (len(canon), len(sci), len(plain)))

    def found(v):
        for d in DECLARED:
            if d == 0:
                if v == 0:
                    return True
                continue
            if abs(v - d) <= 1e-3 * abs(d):
                return True
        for c in canon:
            if c == 0:
                if v == 0:
                    return True
                continue
            if abs(v - c) <= 1e-3 * abs(c):
                return True
        return False

    bad = []
    for v, s in sci + plain:
        if not found(v):
            bad.append((v, s))
    seen = set()
    uniq = []
    for v, s in bad:
        if s in seen:
            continue
        seen.add(s)
        uniq.append((v, s))
    print('未在 canon/DECLARED 中找到：%d' % len(uniq))
    for v, s in sorted(uniq, key=lambda x: -x[0])[:40]:
        print('  %-24s %s' % (s, v))
    return 0


if __name__ == '__main__':
    sys.exit(main())
