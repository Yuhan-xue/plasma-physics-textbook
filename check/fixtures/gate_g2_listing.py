#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可证伪夹具：附录 G.2「印刷输出 == 代码实跑输出」
=========================================================================

为什么需要它：`check/physics_regression.py` 的 O02 只做静态断言（源码里是否
出现 `success`/`status` 等字样），无法证明书里印刷的那段实测输出**真的**是这段
代码跑出来的。本夹具把 G.2 的 Verbatim 代码块抽出来实跑，与紧随其后的印刷输出
块逐行比对——这是「输出不是手写的」唯一硬证据。

三段判定：
  1. extract   : 在 `\\section{德拜屏蔽电势的数值求解}` 与 G.3 节之间必须恰好
                 有 2 个 Verbatim 块（第一个=代码，第二个=印刷输出）。
  2. reproduce : 用 `MPLBACKEND=Agg` 实跑代码，stdout（去尾空格、丢空行）必须
                 与印刷输出**逐行完全相同**。
  3. comparator: 把印刷输出篡改一位数字 / 删掉一行后，比对器必须报出差异
                 ——证明第 2 步不是「永远通过」的空检查。

严格逐行相等，不做数值容差：印刷输出块旁就写着实测环境（NumPy 2.5.3 /
SciPy 1.18.1）。若本机版本不同导致末位数字漂移，那正是需要人看一眼的信号，
不应被容差悄悄吞掉。

用法：
    python check/fixtures/gate_g2_listing.py           # 需要 numpy/scipy/matplotlib
    python check/fixtures/gate_g2_listing.py --keep    # 保留生成的 .py 与 PNG

工作目录固定在仓库内 `tmp/g2-listing-fixture/`（代码里的 `plt.savefig` 会往
当前目录写 PNG；系统临时目录在受限环境不可写）。控制台只打印 ASCII 结论行，
避免 Windows CP936 下 UnicodeEncodeError。
退出码：0 = 全部符合预期；1 = 有不符预期的项。
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MAIN_TEX = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.tex')
OUTDIR = os.path.join(ROOT, 'tmp', 'g2-listing-fixture')

SEC_G2 = r'\section{德拜屏蔽电势的数值求解}'
SEC_G3 = r'\section{朗缪尔波色散曲线绘制}'
BEGIN = r'\begin{Verbatim}'
END = r'\end{Verbatim}'
SCRIPT_REL = os.path.join('code', 'appendix_g', 'debye_bvp.py')
POINTER = r'code/appendix\_g/debye\_bvp.py'


def extract_blocks():
    """返回 printed_lines；结构不符时抛 ValueError。

    结构（v1.26.52 起）：G.2 节里恰好 2 个 Verbatim 块——① 正文骨架（≤20 行，
    代码本体已移入仓库）② 印刷输出；正文必须给出脚本路径指针。
    """
    text = open(MAIN_TEX, encoding='utf-8').read()
    i0 = text.find(SEC_G2)
    i1 = text.find(SEC_G3)
    if i0 < 0 or i1 < 0 or i1 <= i0:
        raise ValueError('G.2/G.3 section markers not found in main tex')
    region = text[i0:i1]
    if POINTER not in region:
        raise ValueError('G.2 does not point at %s' % SCRIPT_REL)
    blocks = re.findall(
        re.escape(BEGIN) + r'[^\n]*\n(.*?)' + re.escape(END),
        region, re.S)
    if len(blocks) != 2:
        raise ValueError(
            'expected exactly 2 Verbatim blocks (skeleton + printed output) in G.2, got %d'
            % len(blocks))
    skeleton = [ln for ln in blocks[0].splitlines() if ln.strip()]
    if len(skeleton) > 20:
        raise ValueError('G.2 skeleton is %d lines (limit 20)' % len(skeleton))
    return blocks[1].splitlines()


def norm(lines):
    """去行尾空白、丢空行——只比较内容，不比较排版留白。"""
    return [ln.rstrip() for ln in lines if ln.strip()]


def compare(expected, actual):
    """返回 (ok, detail, diff)：逐行严格相等。"""
    diff = []
    for i in range(max(len(expected), len(actual))):
        e = expected[i] if i < len(expected) else '<缺行>'
        a = actual[i] if i < len(actual) else '<缺行>'
        if e != a:
            diff.append((i + 1, e, a))
    detail = '%d/%d lines' % (len(actual), len(expected))
    return (not diff), detail, diff


def run_listing(keep):
    """把仓库脚本复制进工作目录后实跑（脚本里的 savefig 落在工作目录）。"""
    os.makedirs(OUTDIR, exist_ok=True)
    src_path = os.path.join(ROOT, SCRIPT_REL)
    if not os.path.exists(src_path):
        return 127, ['<缺少脚本：%s>' % SCRIPT_REL]
    code_path = os.path.join(OUTDIR, 'g2_debye_bvp.py')
    shutil.copyfile(src_path, code_path)
    env = dict(os.environ)
    env['MPLBACKEND'] = 'Agg'
    env['PYTHONIOENCODING'] = 'utf-8'
    proc = subprocess.run(
        [sys.executable, code_path], cwd=OUTDIR, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = proc.stdout.decode('utf-8', 'replace')
    return proc.returncode, out.splitlines()


def versions():
    parts = []
    for mod in ('numpy', 'scipy', 'matplotlib'):
        try:
            parts.append('%s %s' % (mod, __import__(mod).__version__))
        except Exception:
            parts.append('%s ?' % mod)
    return ' / '.join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--keep', action='store_true', help='保留生成的 .py 与 PNG')
    args = ap.parse_args()

    for mod in ('numpy', 'scipy', 'matplotlib'):
        try:
            __import__(mod)
        except ImportError:
            print('[SKIP] %s not installed; fixture not executed' % mod)
            return 0

    bad = 0

    try:
        printed_lines = extract_blocks()
    except ValueError as exc:
        print('[X] extract: %s' % exc)
        return 1
    printed = norm(printed_lines)
    print('[OK] extract: printed_output=%d lines (script: %s)'
          % (len(printed), SCRIPT_REL))

    rc, actual = run_listing(args.keep)
    actual = norm(actual)
    ok, detail, diff = compare(printed, actual)
    print('[%s] reproduce: exit=%d, match=%s (%s)'
          % ('OK' if ok else 'X', rc, detail, versions()))
    if not ok:
        bad += 1
        for lineno, e, a in diff[:5]:
            print('      line %d: printed=%r actual=%r' % (lineno, e, a))

    tampered = list(printed)
    if tampered:
        tampered[0] = re.sub(r'\d', '0', tampered[0], count=len(re.findall(r'\d', tampered[0])))
        if tampered[0] == printed[0]:
            tampered[0] = printed[0] + ' 篡改'
    neg_ok, neg_detail, neg_diff = compare(tampered, actual)
    print('[%s] comparator: tampered-number detected=%s (%d shown)'
          % ('OK' if not neg_ok else 'X', not neg_ok, len(neg_diff)))
    if neg_ok:
        bad += 1

    dropped = printed[:-1] if printed else []
    drop_ok, drop_detail, drop_diff = compare(dropped, actual)
    print('[%s] comparator: dropped-line detected=%s (%d shown)'
          % ('OK' if not drop_ok else 'X', not drop_ok, len(drop_diff)))
    if drop_ok:
        bad += 1

    if not args.keep:
        for name in ('g2_debye_bvp.py', 'debye_shielding.png'):
            path = os.path.join(OUTDIR, name)
            if os.path.exists(path):
                os.remove(path)
        if os.path.isdir(OUTDIR) and not os.listdir(OUTDIR):
            os.rmdir(OUTDIR)

    print('fixture verdict: %s' % ('PASS' if bad == 0 else 'FAIL (%d)' % bad))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
