# -*- coding: utf-8 -*-
"""生成 docs/figures.md：面向读者的图册（按章/附录列出全部插图）。

数据来源：figures/CREDITS.md（一行一张图：文件名 | 中文名 | 作者 | 许可 | 来源页 | 用在）。
缩略图：docs/figures-thumbs/<同名>.jpg（宽 460 以内），由本脚本生成。
用法：python -X utf8 tools/make_figures_gallery.py [--write]
"""
import io
import os
import re
import sys

from PIL import Image

sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CREDITS = os.path.join(ROOT, 'figures', 'CREDITS.md')
THUMBS = os.path.join(ROOT, 'docs', 'figures-thumbs')
OUT = os.path.join(ROOT, 'docs', 'figures.md')
WRITE = '--write' in sys.argv

ROW = re.compile(r'^\| `([^`]+)` \| (.*?) \| (.*?) \| (.*?) \| (\S+) \| (.*?) \|\s*$')


def read_rows():
    rows = []
    with io.open(CREDITS, encoding='utf-8') as f:
        for line in f:
            m = ROW.match(line.rstrip('\n'))
            if m:
                fname, name, author, lic, url, where = m.groups()
                rows.append({'file': fname, 'name': name, 'author': author,
                             'lic': lic, 'url': url, 'where': where.strip()})
    return rows


def order_key(where):
    """第 N 章 → (0, N)；附录 X → (1, 字母序号)；配套手册 J/K/L → (2, ...)"""
    m = re.search(r'第\s*(\d+)\s*章', where)
    if m:
        return (0, int(m.group(1)), '')
    m = re.search(r'附录\s*([A-Z]{1,2})', where)
    if m:
        s = m.group(1)
        n = (ord(s[0]) - ord('A')) * 26 + (len(s) - 1)
        return (1, n, s)
    return (2, 99, where)


def group_title(where):
    m = re.search(r'第\s*(\d+)\s*章', where)
    if m:
        return '第 %s 章' % m.group(1)
    m = re.search(r'附录\s*([A-Z]{1,2})', where)
    if m:
        s = m.group(1)
        if s in ('J', 'K', 'L'):
            return '附录 %s（配套手册）' % s
        return '附录 %s' % s
    return where or '其他'


def make_thumb(row):
    src = os.path.join(ROOT, 'figures', row['file'])
    stem = os.path.splitext(row['file'])[0]
    dst = os.path.join(THUMBS, stem + '.jpg')
    if not os.path.exists(src):
        return None
    if os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        return dst
    if not WRITE:
        return dst
    os.makedirs(THUMBS, exist_ok=True)
    im = Image.open(src)
    if im.mode not in ('RGB', 'L'):
        im = im.convert('RGB')
    im.thumbnail((460, 460), Image.LANCZOS)
    im.save(dst, 'JPEG', quality=80, optimize=True)
    return dst


def main():
    rows = read_rows()
    rows.sort(key=lambda r: order_key(r['where']))
    print('读到 %d 张图' % len(rows))
    groups = {}
    for r in rows:
        groups.setdefault(group_title(r['where']), []).append(r)

    lines = []
    lines.append('# 图册：全书 95 张插图一览')
    lines.append('')
    lines.append('> 这本书用图讲物理：全书 95 张插图，全部是实拍照片或公开示意图，'
                 '逐张登记了作者、许可与来源页。')
    lines.append('> 下面按「章 / 附录」把它们排开，点击任意一张可看原图；'
                 '完整清单与来源链接见 [`figures/CREDITS.md`](../figures/CREDITS.md)。')
    lines.append('')
    lines.append('合计 **%d 张**。若想按来源查看（NASA、ITER、LLNL、EUROfusion…），'
                 '见文末统计。' % len(rows))
    lines.append('')
    for g, items in groups.items():
        if len(items) == 1 and g == '其他':
            continue
        lines.append('## %s（%d 张）' % (g, len(items)))
        lines.append('')
        lines.append('| 图 | 说明 | 作者 | 许可与来源 |')
        lines.append('|---|---|---|---|')
        for r in items:
            thumb = make_thumb(r)
            if thumb is None:
                print('缺图：', r['file'])
                continue
            rel = os.path.relpath(thumb, os.path.join(ROOT, 'docs')).replace('\\', '/')
            img = '<img src="%s" width="230">' % rel
            lic = r['lic']
            link = r['url']
            lines.append('| %s | %s | %s | %s（[来源页](%s)） |'
                         % (img, r['name'], r['author'], lic, link))
        lines.append('')

    # 来源统计（口径统一在 tools/figure_sources.py）
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from figure_sources import counts as source_counts
    c = source_counts(rows)
    lines.append('## 来源统计')
    lines.append('')
    lines.append('- 政府机构、国家实验室与大学/研究所：**%d 张**'
                 '（NASA、ESA/SDO、ITER、LLNL、INL、EUROfusion、USAF、KSTAR 等）' % c['inst'])
    lines.append('- 个人摄影师、企业与开放社区：**%d 张**' % c['person'])
    lines.append('')
    lines.append('许可只用两种：**公有领域**与 **CC BY / CC BY-SA**（可商用、可改、需署名）。')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('生成：`python -X utf8 tools/make_figures_gallery.py --write`'
                 '（数据取自 `figures/CREDITS.md`，缩略图在 `docs/figures-thumbs/`）。')

    text = '\n'.join(lines) + '\n'
    if WRITE:
        io.open(OUT, 'w', encoding='utf-8', newline='\n').write(text)
        print('写出', OUT, len(text), '字符')
    else:
        print('（干跑）将写出 docs/figures.md，%d 字符' % len(text))
        print('\n'.join(lines[:24]))


if __name__ == '__main__':
    main()
