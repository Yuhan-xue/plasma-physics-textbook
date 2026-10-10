#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可证伪夹具：check/pdf_outline.py（第六套「书签／目录／引用一致性守卫」）
=====================================================================

为什么需要它：本项目纪律是「检查通过必须配负例」——一个从不失败的检查
等于没有检查。本夹具在**副本**上篡改 PDF 结构，验证守卫确实会 FAIL：

  1. title : 把某条书签标题改成不存在的名字 → 期望 O3a FAIL（找不到同名项）
  2. order : 把首尾两条书签页码对调造成倒序   → 期望 O2  FAIL（倒退）
  3. empty : 清空全部书签                      → 期望 O1  FAIL（没有任何书签）

不加「页码越界」反例：PyMuPDF 的 `set_toc()` 自身就拒绝写入超出页数的行
（`ValueError: row 0: page number out of range`），无法通过该路径构造出越界
目录；守卫里的越界分支保留为纵深防御，但不由本夹具覆盖。

用法：
    python check/fixtures/gate_pdf_outline.py          # 需要 PyMuPDF
    python check/fixtures/gate_pdf_outline.py --keep   # 保留篡改副本

输出目录固定在仓库内 `tmp/outline-fixture/`（系统临时目录在受限环境不可写），
控制台只打印 ASCII 结论行，避免 Windows CP936 下 UnicodeEncodeError。
退出码：0 = 全部符合预期；1 = 有反例未按预期失败。
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC_PDF = os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf')
CHECKER = os.path.join(ROOT, 'check', 'pdf_outline.py')
OUTDIR = os.path.join(ROOT, 'tmp', 'outline-fixture')


def run_checker(pdf):
    """跑守卫并返回 (exit_code, [{'level','cat','msg'}])。"""
    proc = subprocess.run(
        [sys.executable, CHECKER, '--pdf', pdf, '--json'],
        cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    data = []
    try:
        data = json.loads(proc.stdout.decode('utf-8', 'replace') or '[]')
    except ValueError:
        data = []
    return proc.returncode, data


def has_needle(data, level, needle):
    return any(r.get('level') == level and needle in r.get('msg', '') for r in data)


def tamper(kind, dst):
    import pymupdf
    doc = pymupdf.open(SRC_PDF)
    toc = doc.get_toc()
    if kind == 'title':
        toc[0][1] = '篡改标题XYZ不存在'
    elif kind == 'order':
        toc[0][2], toc[-1][2] = toc[-1][2], toc[0][2]
    elif kind == 'empty':
        toc = []
    doc.set_toc(toc)
    doc.save(dst, garbage=0, deflate=True)
    doc.close()


CASES = [
    ('title', 1, 'FAIL', 'O3a', '找不到同名项'),
    ('order', 1, 'FAIL', 'O2', '倒退'),
    ('empty', 1, 'FAIL', 'O1', '没有任何书签'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--keep', action='store_true', help='保留篡改副本')
    args = ap.parse_args()

    try:
        import pymupdf  # noqa: F401
    except ImportError:
        print('[SKIP] pymupdf not installed; fixture not executed')
        return 0

    os.makedirs(OUTDIR, exist_ok=True)
    bad = 0
    for kind, want_exit, level, cat, needle in CASES:
        dst = os.path.join(OUTDIR, 'tamper_%s.pdf' % kind)
        tamper(kind, dst)
        code, data = run_checker(dst)
        ok = (code == want_exit) and has_needle(data, level, needle)
        cats = sorted({r.get('cat') for r in data if r.get('level') == 'FAIL'})
        print('[%s] %-6s exit=%d want=%d needle=%s FAIL_cats=%s'
              % ('OK' if ok else 'X', kind, code, want_exit, needle, cats))
        if not ok:
            bad += 1
            for r in data:
                if r.get('level') == 'FAIL':
                    print('      %s %s' % (r.get('cat'), r.get('msg')[:100]))
        if not args.keep:
            os.remove(dst)

    print('fixture verdict: %s' % ('PASS' if bad == 0 else 'FAIL (%d)' % bad))
    if not args.keep and os.path.isdir(OUTDIR) and not os.listdir(OUTDIR):
        os.rmdir(OUTDIR)
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
