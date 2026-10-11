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
  4. tail  : 把目录**最后一页**的链接删到只剩 1 个（目录短尾页的真实情形，
             链接数 1 < MIN_TOC_LINKS）        → 期望 O5 WARN 且页跨度仍含该页
             这条锁住 `detect_toc_pages` 的「短尾页吸收」：老实现会把该页判成正文页，
             于是 O5 报出的跨度少一页（`第 2--N-1 页`），needle 对不上 → 夹具失败。

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

sys.path.insert(0, os.path.join(ROOT, 'check'))
from pdf_outline import MIN_TOC_LINKS, TOC_FIRST      # noqa: E402  与守卫共用同一口径


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
    """在副本上植入缺陷；返回该用例额外的 needle 参数（没有则 None）。"""
    import pymupdf
    doc = pymupdf.open(SRC_PDF)
    toc = doc.get_toc()
    if kind == 'title':
        toc[0][1] = '篡改标题XYZ不存在'
    elif kind == 'order':
        toc[0][2], toc[-1][2] = toc[-1][2], toc[0][2]
    elif kind == 'empty':
        toc = []
    elif kind == 'toc-tail':
        # 先按守卫的口径找印刷目录页前缀（第 2 页起连续「链接数 ≥12」），
        # 再把最后一页的链接删到只剩 1 个 —— 复现「目录短尾页」。
        n = 0
        for pno in range(TOC_FIRST, doc.page_count):
            if len(doc[pno].get_links()) >= MIN_TOC_LINKS:
                n += 1
            else:
                break
        last = TOC_FIRST + n - 1                      # 0-based：目录末页
        page = doc[last]
        for lk in page.get_links()[1:]:
            page.delete_link(lk)
        doc.set_toc(toc)
        doc.save(dst, garbage=0, deflate=True)
        doc.close()
        return '第 %d--%d 页' % (TOC_FIRST + 1, last + 1)
    doc.set_toc(toc)
    doc.save(dst, garbage=0, deflate=True)
    doc.close()
    return None


CASES = [
    ('title', 1, 'FAIL', 'O3a', '找不到同名项'),
    ('order', 1, 'FAIL', 'O2', '倒退'),
    ('empty', 1, 'FAIL', 'O1', '没有任何书签'),
    ('toc-tail', 0, 'WARN', 'O5', None),      # needle 由 tamper() 返回（含目录跨度）
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
        extra = tamper(kind, dst)
        if needle is None:
            needle = extra or ''
        code, data = run_checker(dst)
        ok = (code == want_exit) and has_needle(data, level, needle) \
            and (cat is None or any(r.get('cat') == cat for r in data if r.get('level') == level))
        cats = sorted({r.get('cat') for r in data if r.get('level') == 'FAIL'})
        print('[%s] %-8s exit=%d want=%d needle=%s FAIL_cats=%s'
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
