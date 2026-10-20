# -*- coding: utf-8 -*-
"""把 PDF 书签页与目录链接目标改到印刷目录写的那一页（两遍法，避免边改边读）。

为什么需要：hyperref 的书签锚点在 \\section 排版时取值，标题被挤到下一页时
书签就早一页；印刷目录里的页码来自 .aux，与实际排版一致。于是「点目录会落到
前一页」。这里按印刷目录的页码，把书签页与目录链接目标一起改正。

用法：
    python tools/fix_outline_pages.py pdf/xxx.pdf [--out 新文件]
不给 --out 时原地改写（先写临时文件再替换）。
"""
import argparse
import os
import re
import sys

import pymupdf

sys.stdout.reconfigure(encoding='utf-8')

TOC_FIRST = 1          # 目录从第 2 页（索引 1）开始
MIN_TOC_LINKS = 20


def norm(s):
    s = re.sub(r'^(?:[A-Za-z]{1,2}[.．])?\d+(?:[.．]\d+)*[.．]?(?![-\u2010-\u2015])', '', s or '')
    s = re.sub(r'[\s\u200b]+', '', s)
    s = (s.replace('“', '“').replace('”', '“').replace('″', '“')
         .replace('‘', "'").replace('’', "'"))
    return re.sub(r'[^0-9A-Za-z\u4e00-\u9fff"\'()（）,，、:：;；.。!！?？\-–—/]', '', s)


def printed_page(page, rect, words):
    """题注行右侧右对齐的纯数字 = 印刷目录页码。"""
    yc = (rect.y0 + rect.y1) / 2.0
    best = None
    for w in words:
        if not re.fullmatch(r'\d{1,3}', w[4]):
            continue
        if abs((w[1] + w[3]) / 2.0 - yc) > 6.0:
            continue
        if w[0] < rect.x1 - 1.0:
            continue
        if best is None or w[2] > best[2]:
            best = w
    return int(best[4]) if best else None


def count_toc_pages(doc):
    """目录页数：先数连续「每页 ≥20 个链接」的前缀，再补短尾页。

    短尾页（目录最后只排几行）链接数可能低于阈值，但它仍是目录：判据是链接
    矩形内的文字能与某个书签同名、且同一行右侧有页码数字。
    """
    n = 0
    for pno in range(TOC_FIRST, doc.page_count):
        if len(doc[pno].get_links()) >= MIN_TOC_LINKS:
            n += 1
        else:
            break
    titles = {norm(t) for _l, t, _p in doc.get_toc()}
    for pno in range(TOC_FIRST + n, min(TOC_FIRST + n + 3, doc.page_count)):
        page = doc[pno]
        words = page.get_text('words')
        hit = False
        for lk in page.get_links():
            rect = lk.get('from')
            if rect is not None and norm(page.get_textbox(rect)) in titles \
                    and printed_page(page, rect, words):
                hit = True
                break
        if hit:
            n += 1
    return max(n, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf')
    ap.add_argument('--out')
    args = ap.parse_args()

    doc = pymupdf.open(args.pdf)
    toc = doc.get_toc()
    normed = [norm(t) for _l, t, _p in toc]
    span = list(range(min(count_toc_pages(doc) + TOC_FIRST, doc.page_count)))

    # 第一遍：只读，收集 (页号, 链接矩形, 标题, 印刷页, 现有目标页)
    entries = []
    for pno in span:
        page = doc[pno]
        links = page.get_links()
        if not links:
            continue
        words = page.get_text('words')
        for lk in links:
            rect = lk.get('from')
            if rect is None:
                continue
            entries.append([pno, rect, page.get_textbox(rect).strip(),
                            printed_page(page, rect, words), lk.get('page')])

    # 配对书签（同名条目按文档顺序）
    cursor = 0
    plan = []          # 每页要重建的链接：(pno, rect, 目标页 1-based)
    sec_page = {}      # 小节号 -> 正确页（供正文里的 \S\ref 链接复用）
    fixed_toc = 0
    for pno, rect, title, printed, target in entries:
        key = norm(title)
        idx = next((i for i in range(cursor, len(normed)) if normed[i] == key), None)
        if idx is None:
            idx = next((i for i in range(len(normed)) if normed[i] == key), None)
            if idx is None:
                continue
        cursor = idx + 1
        want = printed or (toc[idx][2])
        if toc[idx][2] != want:
            toc[idx][2] = want
            fixed_toc += 1
        plan.append((pno, rect, want))
        num = re.match(r'((?:[A-Za-z]{1,2}\.)?\d+(?:\.\d+)*)', title.strip())
        if num:
            sec_page[num.group(1).rstrip('.')] = want

    # 第二遍：清空目录页旧链接，按印刷位置（先上后下、先左后右）重建。
    # 必须整页重建：只改落点不对的那几条会把注释顺序打乱，目录行序跟着乱，
    # 检查脚本按顺序配对同名条目时就会误判。
    by_page = {}
    for pno, rect, want in plan:
        by_page.setdefault(pno, []).append((rect, want))
    fixed_links = 0
    for pno, items in by_page.items():
        page = doc[pno]
        for lk in page.get_links():
            page.delete_link(lk)
        items.sort(key=lambda it: (round(it[0].y0, 1), round(it[0].x0, 1)))
        for rect, want in items:
            page.insert_link({'kind': pymupdf.LINK_GOTO, 'from': rect,
                              'page': want - 1, 'to': pymupdf.Point(0, 0), 'zoom': 0})
            fixed_links += 1

    doc.set_toc(toc)
    out = args.out or (args.pdf + '.fixed')
    doc.save(out, garbage=3, deflate=True)
    doc.close()

    if not args.out:
        os.replace(out, args.pdf)
        out = args.pdf
    print('entries %d / 书签 %d；改书签 %d 条、重建目录链接 %d 条 -> %s'
          % (len(plan), len(toc), fixed_toc, fixed_links, out))


if __name__ == '__main__':
    main()
