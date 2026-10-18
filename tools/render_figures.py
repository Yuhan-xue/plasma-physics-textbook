#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把书里的自绘图渲染成 PNG，连图注一起裁出来——给"人眼复核"用的那一半。
====================================================================

为什么需要它：`check/figure_layout.py` 只告诉你"几何判据没过"，但不告诉你
"长成什么样"。重绘的人需要一个 200 DPI、边界正确、带图注的小图，才能判断是
真坏了还是判据太挑。所以本脚本复用同一套定位与图区算法（import check/figure_layout），
保证"量的"和"看的"是同一块地方，不会出现检查说 A、截图截到 B。

裁切范围：左/右 = 版心（图若越出版心就跟着放宽，好让人看见越界），
上 = 图区上沿再留 3 pt，下 = 题注末行再留 3 pt——这样"图有没有压到题注"
一眼可判。

用法
----
    python -X utf8 tools/render_figures.py --all                    # 全书 36 张
    python -X utf8 tools/render_figures.py --label fig:glow          # 单张（可重复）
    python -X utf8 tools/render_figures.py --all --dpi 300 --out tmp/figcheck300

参数：
    --all           渲染全部自画图（与 --label 二选一，可同时给）
    --label LABEL   只渲染该 label，可重复
    --dpi N         渲染分辨率，默认 200
    --out DIR       PNG 目录，默认 tmp/figcheck
    --pdf PATH      编译产物（默认 pdf/plasma_physics_textbook_v2.pdf）
    --src DIR       书稿目录（默认 src）

每张图打印一行摘要：label | 图号 | 页 | 宽×高（pt 与 cm）| 占版心 | 是否超 86% 版心
文件名：`<figkey>.png`（figkey = label 去掉 `fig:`；无 label 的图用 `文件-L行号`）。

退出码：0 —— 渲染是给人看的辅助手段，不承担门禁职责（几何判据在 check 里）。
"""

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'check'))
import figure_layout as fl                                   # noqa: E402

DEFAULT_OUT = os.path.join(ROOT, 'tmp', 'figcheck')


def figkey(fig):
    """文件名：label 去前缀；没有 label 就用「文件-L行号」。"""
    lab = fig.get('label')
    if lab:
        return lab[4:] if lab.startswith('fig:') else lab
    stem = os.path.splitext(fig['file'])[0]
    return '%s-L%d' % (stem, fig['line'])


def clip_rect(page, box, ent):
    """图区 + 题注的渲染窗口（左右按版心，越界时跟着放宽）。"""
    x0 = min(fl.TEXT_L, box[0]) - 3.0
    x1 = max(fl.TEXT_R, box[2]) + 3.0
    y0 = max(0.0, box[1] - 3.0)
    y1 = (ent['y1'] if ent else box[3]) + 3.0
    x0 = max(0.0, x0)
    x1 = min(page.rect.width, x1)
    y1 = min(page.rect.height, y1)
    return x0, y0, x1, y1


def main():
    ap = argparse.ArgumentParser(description='渲染自画图为 PNG（含图注）')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--label', action='append', default=None, metavar='LABEL')
    ap.add_argument('--dpi', type=int, default=200)
    ap.add_argument('--out', default=DEFAULT_OUT)
    ap.add_argument('--pdf', default=fl.DEFAULT_PDF)
    ap.add_argument('--src', default=fl.DEFAULT_SRC)
    args = ap.parse_args()

    if not args.all and not args.label:
        ap.error('要么 --all，要么 --label fig:xxx')

    try:
        import pymupdf as fitz
    except ImportError:
        print('未安装 PyMuPDF（pip install PyMuPDF）')
        return 0

    if not os.path.exists(args.pdf):
        print('找不到编译产物：%s' % args.pdf)
        return 0
    if not os.path.isdir(args.out):
        os.makedirs(args.out)

    figs = [f for f in fl.discover_figures(args.src) if not f.get('skip')]
    if args.label:
        want = [s for s in args.label]
        figs = [f for f in figs
                if any(s in (f.get('label') or '') or s == figkey(f) for s in want)]

    doc = fitz.open(args.pdf)
    aux = fl.parse_aux_labels(os.path.join(args.src, fl.AUX_NAME))
    by_num, every = fl.build_caption_index(doc)

    n_ok = n_bad = 0
    big = []
    print('%-32s %-6s %-4s %-16s %-7s %-6s %s'
          % ('label', '图号', '页', '宽×高（cm）', '占版心', '结论', 'PNG'))
    for f in figs:
        page_no, ent, method = fl.locate_figure(f, by_num, every, aux)
        key = figkey(f)
        if page_no is None:
            n_bad += 1
            print('%-32s %-6s %-4s %-16s %-7s %-6s %s'
                  % (f.get('label') or key, '—', '—', '—', '—', '未定位', method))
            continue
        page = doc[page_no - 1]
        box, n_items = fl.region_for_caption(page, ent) if ent is not None else (None, 0)
        if box is None:
            n_bad += 1
            print('%-32s %-6s %-4s %-16s %-7s %-6s %s'
                  % (f.get('label') or key, ent['num'], page_no, '—', '—', '空图区', '图区没内容'))
            continue
        w_pt = box[2] - box[0]
        h_pt = box[3] - box[1]
        w_frac = w_pt / fl.TEXT_W
        r = clip_rect(page, box, ent)
        pix = page.get_pixmap(clip=fitz.Rect(*r), dpi=args.dpi)
        path = os.path.join(args.out, key + '.png')
        pix.save(path)
        n_ok += 1
        flag = '偏大' if w_frac > fl.W_BIG else '正常'
        if w_frac > fl.W_BIG or h_pt / fl.PT_CM > fl.H_HUGE:
            big.append((f.get('label') or key, ent['num'], page_no,
                        w_pt / fl.PT_CM, h_pt / fl.PT_CM, w_frac))
        print('%-32s %-6s %-4d %-16s %-7s %-6s %s'
              % (f.get('label') or key, ent['num'], page_no,
                 '%.2f×%.2f' % (w_pt / fl.PT_CM, h_pt / fl.PT_CM),
                 '%.0f%%' % (w_frac * 100), flag,
                 os.path.relpath(path, ROOT)))

    print()
    print('渲染 %d 张，未完成 %d 张；PNG 在 %s（%d DPI）'
          % (n_ok, n_bad, os.path.relpath(args.out, ROOT), args.dpi))
    if big:
        print('偏大/过高的图（> %.0f%% 版心 或 > %.0f cm）：' % (fl.W_BIG * 100, fl.H_HUGE))
        for lab, num, pno, wcm, hcm, wf in big:
            print('  %-32s 图 %-6s p%-4d 宽 %5.2f cm（%.0f%%） 高 %5.2f cm'
                  % (lab, num, pno, wcm, wf * 100, hcm))
    else:
        print('没有偏大的图（全部 ≤ %.0f%% 版心、≤ %.0f cm）' % (fl.W_BIG * 100, fl.H_HUGE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
