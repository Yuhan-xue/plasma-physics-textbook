# -*- coding: utf-8 -*-
"""附录 AG 交叉核对：正文里的每个数字都要能在 canon 或 DECLARED 里找到。"""
import io
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
TEX = os.path.join(ROOT, 'src', 'appendix_ag_equilibrium.tex')
CANON = os.path.join(ROOT, 'tmp', 'appAG_calc.txt')

# 声明级输入与跨附录引用值（正文里写明是声明值/别处口径）
DECLARED = {
    1.7, 0.4, 1.0, 1.0, 2.5, 3.0, 20.0, 2.8, 3.5, 1.8,  # 声明与文献量级
    1.31084, 0.819273, 1.2, 2.54715, 0.85,                # 附录 AD / AE 引用
    1.2566e-6, 1.5e7, 8e19, 5.5, 1.6,                     # 常数写法与练习输入
    157.0, 213.0,                                          # 文献页码/卷号
    9194.0, 1e20, 6.2, 2.0, 5.3, 15.0, 0.5,               # 工作点输入
    832.2, 489.532, 4.0, 12.0, 95.0, 2.2, 3.1,            # 几何与文字性数
    1984.0, 1999.0, 2137.0, 1958.0, 1966.0, 1973.0, 1964.0, 1971.0,
    1985.0, 1983.0, 1996.0, 2014.0, 2011.0, 1968.0, 1981.0, 1997.0,
    31.0, 2.0, 6.0, 8.0, 11.0, 12.0, 14.0, 18.0, 26.0, 39.0, 97.0, 209.0,
    1611.0, 3553.0, 219.0, 671.0, 103.0, 151.0, 190.0, 5.0, 4.0,
}


def ensure_canon():
    if os.path.exists(CANON):
        return
    calc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'appAG_calc.py')
    r = subprocess.run([sys.executable, '-X', 'utf8', calc], capture_output=True,
                       text=True, encoding='utf-8')
    io.open(CANON, 'w', encoding='utf-8').write(r.stdout)
    print('已生成 ' + CANON)


def canon_numbers():
    if not os.path.exists(CANON):
        print('缺少 %s：先跑 python -X utf8 tools/calc/appAG_calc.py > tmp/appAG_calc.txt' % CANON)
        sys.exit(2)
    t = io.open(CANON, encoding='utf-8').read()
    vals = set()
    for m in re.finditer(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', t):
        try:
            vals.add(float(m.group(0)))
        except ValueError:
            pass
    return vals


def strip_tex(text):
    text = re.sub(r'(?<!\\)%.*', '', text)
    text = re.sub(r'p\{[^}]*\}', '', text)
    text = re.sub(r'\\label\{[^}]*\}', '', text)
    text = re.sub(r'\\(?:eqref|ref|cite)\{[^}]*\}', '', text)
    text = re.sub(r'\\includegraphics(?:\[[^\]]*\])?\{[^}]*\}', '', text)
    return text


def prose_values(text):
    out = []
    for m in re.finditer(r'(\d+(?:\.\d+)?)\\times10\^\{?(-?\d+)\}?', text):
        out.append((float(m.group(1)) * 10 ** int(m.group(2)), m.group(0)))
    tmp = re.sub(r'\d+(?:\.\d+)?\\times10\^\{?-?\d+\}?', ' ', text)
    for m in re.finditer(r'(?<![\d.])(\d+\.\d+|\d{2,})(?:[eE]([-+]?\d+))?(?![\d.])', tmp):
        v = float(m.group(1)) * (10 ** int(m.group(2)) if m.group(2) else 1.0)
        out.append((v, m.group(0)))
    return out


def main():
    ensure_canon()
    raw = io.open(TEX, 'rb').read()
    t = strip_tex(raw.decode('utf-8'))
    blocks = re.split(r'\\section\{', t)
    cnum = canon_numbers()
    bad = []
    seen = set()
    total = 0
    for b in blocks[1:]:
        title = b.split('}')[0]
        if '来源' in title:
            continue
        for v, s in prose_values(b):
            total += 1
            if v in seen:
                continue
            seen.add(v)
            if v in DECLARED:
                continue
            if any(abs(v - c) <= 1e-3 * max(abs(c), 1e-30) for c in cnum):
                continue
            bad.append((title[:14], v, s))
    print('字节数 %d LF=%s 末尾换行=%s' % (len(raw), b'\r' not in raw, raw.endswith(b'\n')))
    print('扫描数字 %d 处（去重 %d）' % (total, len(seen)))
    print('未在 calc/DECLARED 中找到：%d' % len(bad))
    for title, v, s in bad:
        print('   %-16s %-22s <- %s' % (title, repr(v), s))
    return 0 if not bad else 1


if __name__ == '__main__':
    raise SystemExit(main())
