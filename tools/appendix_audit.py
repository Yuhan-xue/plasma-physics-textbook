# -*- coding: utf-8 -*-
"""附录体检：每个附录的页数、体量、图数、节数，以及与其他附录的主题重叠线索。

数据来源：PDF 书签（附录起始页）+ src/*.tex 源码统计。
用法：python -X utf8 tools/appendix_audit.py [--json]
"""
import io
import json
import os
import re
import sys

import pymupdf

sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
PDF = os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf')

main = io.open(os.path.join(SRC, 'plasma_physics_textbook_v2.tex'), encoding='utf-8').read()

# 主文件里的附录顺序 → 源文件（A–H 的正文直接写在主文件里，没有 \input）
entries = []          # (字母, 标题, 源文件名或 None)
inline_spans = []     # 内联附录的字节区间
marks = [(m.start(), m.group(1), m.group(2), m.group(4))
         for m in re.finditer(
             r'\\appchapter\{附录\s*([A-Z]{1,2})[：:]([^}]*)\}(\s*\\input\{([^}]+)\})?', main)]
for i, (pos, letter, title, fname) in enumerate(marks):
    end = marks[i + 1][0] if i + 1 < len(marks) else len(main)
    entries.append((letter, title, (fname or '').strip() or None))
    inline_spans.append((pos, end))

# PDF 书签里附录的起始页
doc = pymupdf.open(PDF)
starts = {}
for lvl, title, page in doc.get_toc():
    m = re.match(r'附录\s*([A-Z]{1,2})[：:]', title.strip())
    if m and m.group(1) not in starts:
        starts[m.group(1)] = page
total = doc.page_count

rows = []
for i, (letter, title, fname) in enumerate(entries):
    p0 = starts.get(letter)
    nxt = None
    for j in range(i + 1, len(entries)):
        if starts.get(entries[j][0]):
            nxt = starts[entries[j][0]]
            break
    p1 = (nxt - 1) if nxt else total
    pages = (p1 - p0 + 1) if p0 else None
    if fname:
        path = os.path.join(SRC, fname + '.tex')
        t = io.open(path, encoding='utf-8').read() if os.path.exists(path) else ''
    else:
        t = main[inline_spans[i][0]:inline_spans[i][1]]
        fname = '（主文件内）'
    rows.append({
        'letter': letter,
        'title': title,
        'file': fname,
        'pages': pages,
        'bytes': len(t.encode('utf-8')),
        'sections': len(re.findall(r'\\section\{', t)),
        'figures': len(re.findall(r'\\includegraphics', t)),
        'tables': len(re.findall(r'\\begin\{(table|longtable)\}', t)),
        'boxes': len(re.findall(r'\\begin\{(tip|insight|warning)\}', t)),
        'exercises': len(re.findall(r'\\begin\{example\}', t)),
    })

if len(rows) < 40 or not any(r['pages'] for r in rows):
    raise SystemExit('解析失败：只认出 %d 个附录（检查 \\appchapter 的写法是否变了）' % len(rows))
rows_by_pages = sorted([r for r in rows if r['pages']], key=lambda r: -r['pages'])
print('附录总数 %d；有书签页码的 %d 个' % (len(rows), len(rows_by_pages)))
print('%-4s %-30s %5s %8s %4s %4s %4s %4s' %
      ('字母', '标题', '页数', '字节', '节', '图', '表', '盒'))
for r in rows_by_pages:
    print('%-4s %-30s %5d %8d %4d %4d %4d %4d' %
          (r['letter'], r['title'][:28], r['pages'], r['bytes'],
           r['sections'], r['figures'], r['tables'], r['boxes']))

tot_pages = sum(r['pages'] for r in rows_by_pages)
print('附录合计页数 %d（占全书 %d 页的 %.0f%%）' % (tot_pages, total, 100.0 * tot_pages / total))
print('前 10 大附录占附录页数的 %.0f%%' %
      (100.0 * sum(r['pages'] for r in rows_by_pages[:10]) / tot_pages))

if '--json' in sys.argv:
    out = os.path.join(ROOT, 'tools', 'appendix_audit.json')
    io.open(out, 'w', encoding='utf-8').write(json.dumps(rows, ensure_ascii=False, indent=1))
    print('写出', out)
