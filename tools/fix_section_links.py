# -*- coding: utf-8 -*-
r"""把附录 F 里「\S\ref」的链接落点校正到 .aux 记录的小节首页。

hyperref 的锚点在 \section 排版时取值；标题被挤到下一页时，链接会落到前一页底部。
这里按 .aux 的 label 页码重设这些链接的目标页，页内链接按 (y, x) 重建以保持顺序。

用法：python -X utf8 tools/fix_section_links.py <in.pdf> [out.pdf]
"""
import io
import os
import re
import sys

import pymupdf

sys.stdout.reconfigure(encoding='utf-8')

MAIN_TEX = 'src/plasma_physics_textbook_v2.tex'
AUX = 'src/plasma_physics_textbook_v2.aux'
TOC_MARK = '附录'


def norm(s):
    return re.sub(r'\s+', ' ', s).strip()


def main():
    src_pdf = sys.argv[1] if len(sys.argv) > 1 else 'pdf/plasma_physics_textbook_v2.pdf'
    out_pdf = sys.argv[2] if len(sys.argv) > 2 else src_pdf

    aux = {}
    text = io.open(AUX, encoding='utf-8', errors='replace').read()
    for m in re.finditer(r'\\newlabel\{([^}]*)\}\{\{([^}]*)\}\{([^}]*)\}', text):
        aux[m.group(1)] = int(m.group(3)) if m.group(3).isdigit() else None

    tex = io.open(MAIN_TEX, encoding='utf-8', errors='replace').read()
    m = re.search(r'\\appchapter\{附录\s*F.*?(?=\\appchapter\{附录\s*G)', tex, re.S)
    if not m:
        raise SystemExit('找不到附录 F 源区段')
    labels = re.findall(r'\\S\\ref\{([^}]+)\}', m.group(0))

    doc = pymupdf.open(src_pdf)
    toc = doc.get_toc()

    def first_page(ch):
        for lvl, t, p in toc:
            if lvl == 1 and re.sub(r'\s+', '', norm(t)).startswith(TOC_MARK + ch):
                return p
        return None

    f_start, g_start = first_page('F'), first_page('G')
    if f_start is None or g_start is None:
        raise SystemExit('书签里找不到附录 F/G')
    f_end = g_start - 1

    pages = list(range(f_start - 1, min(f_end, doc.page_count)))
    entries = []          # (page_index, y, x, rect, kind, target)
    sec_links = []        # 与 labels 顺序对应的条目
    for pno in pages:
        page = doc[pno]
        for lk in page.get_links():
            rect = lk.get('from')
            if rect is None:
                continue
            txt = page.get_textbox(rect).strip()
            is_sec = bool(re.match(r'^[§\s]*\d+(?:\.\d+)*', txt))
            entries.append([pno, round(rect.y0, 1), round(rect.x0, 1), dict(lk)])
            if is_sec:
                sec_links.append(len(entries) - 1)

    if len(sec_links) != len(labels):
        raise SystemExit('附录 F 的 \\S\\ref 数量（%d）与链接数（%d）不一致'
                         % (len(labels), len(sec_links)))

    fixed = 0
    for slot, lab in zip(sec_links, labels):
        want = aux.get(lab)
        if want is None:
            continue
        cur = entries[slot][3].get('page')
        if cur != want - 1:
            entries[slot][3]['page'] = want - 1
            entries[slot][3]['kind'] = pymupdf.LINK_GOTO
            entries[slot][3].pop('nameddest', None)
            fixed += 1

    # 整页重建：先按 (y, x) 排序，再清空、重插，保持注解顺序
    for pno in pages:
        page = doc[pno]
        page_entries = sorted([e for e in entries if e[0] == pno], key=lambda e: (e[1], e[2]))
        for lk in page.get_links():
            page.delete_link(lk)
        for _p, _y, _x, lk in page_entries:
            page.insert_link(lk)

    # 整页重建后整份写出（不能用 incremental：上游工具已重写过文件，加密信息会冲突）
    if os.path.realpath(out_pdf) == os.path.realpath(src_pdf):
        tmp_out = src_pdf + '.tmp'
        doc.save(tmp_out)
        doc.close()
        os.replace(tmp_out, src_pdf)
    else:
        doc.save(out_pdf)
        doc.close()
    print('附录 F 页 %d–%d，链接 %d 个（其中 \\S\\ref 类 %d 个），改目标页 %d 条 -> %s'
          % (f_start, f_end, len(entries), len(sec_links), fixed, out_pdf))


if __name__ == '__main__':
    main()
