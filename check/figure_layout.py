#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》图形版式自动守卫（基于编译后 PDF 的几何）
==================================================================

背景（PP-U01）：图 9.1／9.2 曾经出现节点互相压盖、标签被框边裁切的问题，
而旧的检查只看句长与颜色，看不出版式崩溃。本模块直接量 PDF 的几何：

  P. 版式几何 —— 对指定图形的页面：
      1) 填充节点框两两相交面积 > 1 pt²            → FAIL（节点交叠）
      2) 文字与某节点框相交 >5%，但落在框内 <95%    → FAIL（标签被框压住/裁切）
      3) 一块文字同时落进两个节点框（面积和 >60%）  → FAIL（文字跨框）
      4) 节点框越出版心（A4，left 1.75 cm / right 1.95 cm）→ FAIL（超宽）
      5) 页面上找不到目标图形                      → FAIL（图形消失，不能静默通过）

  Q. 示意图几何（DIAGRAMS：没有节点框的图，如瑞利--泰勒）：
      6) 两块文字互相压住（交叠面积/较小者 > 50%）  → FAIL（文字叠字）
      7) 彩色描边路径（箭头、连线）穿过文字包围盒   → FAIL（线条穿字）

  规则 6/7 的由来：v1.7.0 目视检查第 114 页时发现图 9.3 的红色重力箭头与
  $\vect{g}$ 标签正压在「重流体 ($\rho_2$)」区域标签上（箭杆自 x=290.6 pt 竖直
  穿过标签包围盒 [254.9, 583.9, 326.2, 600.5]），而原有的节点框规则查不到
  这类「无节点框示意图」的叠字。规则 6 只比 word 级包围盒，不做数学上下标
  拆分（否则全书会误报：PyMuPDF 会把 $v_t^2$ 拆成 `v2` 与 `t` 两个 word）。

依赖 PyMuPDF；未安装时报 WARN 并跳过几何项（不会静默通过）。

用法：
    python check/figure_layout.py                      # 用 pdf/plasma_physics_textbook_v2.pdf
    python check/figure_layout.py --pdf <path> --json

退出码：0 = 无 FAIL；1 = 存在 FAIL
"""

import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PDF = os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf')

# 版心：A4 宽 21 cm，geometry 设 left=1.75 cm、right=1.95 cm（v1.9.0 版面密度规格）
PT_CM = 72.0 / 2.54
PAGE_W = 21.0 * PT_CM
TEXT_L = 1.75 * PT_CM
TEXT_R = PAGE_W - 1.95 * PT_CM

# 待守护的图形：题注中的特征词（PDF 文本层，必须足够独特）→ 说明
FIGURES = [
    ('正反馈回路', '图 9.1 漂移波正反馈回路'),
    ('分类树状图', '图 9.2 等离子体不稳定性分类树'),
]

# 没有「节点框」的示意图：改查文字叠字与线条穿字（规则 6/7）
DIAGRAMS = [
    ('尖峰处重流体下沉', '图 9.3 瑞利-泰勒不稳定性物理图像'),
]

ALL_FIGURES = FIGURES + DIAGRAMS

RESULTS = []


def rec(level, cat, msg):
    RESULTS.append((level, cat, msg))


def inter_area(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if (w > 0 and h > 0) else 0.0


def area(a):
    return max(a[2] - a[0], 0.0) * max(a[3] - a[1], 0.0)


def node_boxes(page, ytop, ybot):
    """页面上一块区域内的「节点框」：有填充、尺寸像节点（不是整页底纹）。"""
    out = []
    for d in page.get_drawings():
        f = d.get('fill')
        if f is None or d.get('type') not in ('f', 'fs'):
            continue
        if all(c > 0.995 for c in f[:3]):          # 纯白填充不算
            continue
        r = d['rect']
        box = (r.x0, r.y0, r.x1, r.y1)
        a = area(box)
        if a < 20 or a > 6000:
            continue
        if (box[2] - box[0]) < 10 or (box[3] - box[1]) < 8:
            continue
        if box[3] < ytop or box[1] > ybot:         # 只取图形所在纵向区间
            continue
        out.append(box)
    return out


def analyze(page, key, tag, ytop, ybot):
    h = page.rect.height
    t = page.get_text()
    if key not in t:
        rec('FAIL', 'P', '%s：在本页文本层找不到关键字「%s」（图形可能已被删除或改名）' % (tag, key))
        return False

    boxes = node_boxes(page, ytop, ybot)
    if len(boxes) < 3:
        # 题注关键字可能只是出现在插图目录/正文引用里，那种页面没有图形几何
        rec('INFO', 'P', '%s：第 %d 页含关键字但只识别到 %d 个节点框，跳过几何判定'
            % (tag, page.number + 1, len(boxes)))
        return False

    # 1) 节点交叠
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a = inter_area(boxes[i], boxes[j])
            if a > 1.0:
                bad.append((a, i, j))
    if bad:
        bad.sort(reverse=True)
        rec('FAIL', 'P', '%s：%d 对节点框交叠（最严重 %.1f pt²）'
            % (tag, len(bad), bad[0][0]))
    else:
        rec('PASS', 'P', '%s：%d 个节点框两两不交叠（阈值 1 pt²）' % (tag, len(boxes)))

    # 2)/(3) 文字与框的关系（只看图形区间内的文字，避免误伤正文）
    clipped, crossing = [], []
    for w in page.get_text('words'):
        wb = (w[0], w[1], w[2], w[3])
        if wb[3] < ytop or wb[1] > ybot:
            continue
        wa = area(wb)
        if wa <= 0:
            continue
        hits = []
        for b in boxes:
            ia = inter_area(wb, b)
            if ia / wa > 0.05:
                hits.append((ia / wa, b))
        if not hits:
            continue
        ratio, bb = max(hits, key=lambda t: t[0])
        # 字形包围盒的 ascent 可能比节点填充框高出不到 1 pt（实测 0.53 pt），
        # 那是量测噪声而不是裁切；只有「落在框内不足 90% 且真的探出框外 >1 pt」才算裁切。
        over = max(bb[0] - wb[0], wb[2] - bb[2], bb[1] - wb[1], wb[3] - bb[3], 0.0)
        if ratio < 0.90 and over > 1.0:
            clipped.append((w[4], ratio, over))
        if sum(r for r, _ in hits) > 0.60 and len(hits) > 1:
            crossing.append((w[4], sum(r for r, _ in hits)))
    if clipped:
        rec('FAIL', 'P', '%s：%d 处文字被节点框边裁切（例：%s 落在框内 %.0f%%，探出 %.1f pt）'
            % (tag, len(clipped), clipped[0][0], clipped[0][1] * 100, clipped[0][2]))
    else:
        rec('PASS', 'P', '%s：没有文字被框边裁切（文字要么整块在框内、要么整块在框外）' % tag)
    if crossing:
        rec('FAIL', 'P', '%s：%d 处文字同时压住两个节点框（例：%s 覆盖 %.0f%%）'
            % (tag, len(crossing), crossing[0][0], crossing[0][1] * 100))
    else:
        rec('PASS', 'P', '%s：没有文字跨在两个节点框上' % tag)

    # 4) 版心
    over = [b for b in boxes if b[0] < TEXT_L - 1.0 or b[2] > TEXT_R + 1.0]
    if over:
        worst = max(over, key=lambda b: max(TEXT_L - b[0], b[2] - TEXT_R))
        rec('FAIL', 'P', '%s：%d 个节点框越出版心 [%.1f, %.1f] pt（最远 %.1f pt）'
            % (tag, len(over), TEXT_L, TEXT_R, max(TEXT_L - worst[0], worst[2] - TEXT_R)))
    else:
        rec('PASS', 'P', '%s：全部节点框在版心内（%.1f–%.1f pt，即 1.75/1.95 cm 页边距）'
            % (tag, TEXT_L, TEXT_R))

    rec('INFO', 'P', '%s：本页第 %d 页，识别节点框 %d 个，图形纵向区间 %.0f–%.0f pt'
        % (tag, page.number + 1, len(boxes), ytop, ybot))
    return True


def analyze_diagram(page, key, tag, ytop, ybot):
    """示意图（无节点框）：查「文字叠字」与「彩色线条穿字」。"""
    if key not in page.get_text():
        rec('FAIL', 'P', '%s：在本页文本层找不到关键字「%s」（图形可能已被删除或改名）' % (tag, key))
        return False

    words = [w for w in page.get_text('words') if w[1] >= ytop and w[3] <= ybot]
    if not words:
        rec('FAIL', 'P', '%s：纵向区间 %.0f–%.0f pt 内没有任何文字（图形可能已被删除）'
            % (tag, ytop, ybot))
        return False

    # 6) 文字互相压住
    bad = []
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            a, b = words[i], words[j]
            ia = inter_area(a[:4], b[:4])
            if ia <= 0:
                continue
            amin = min(area(a[:4]), area(b[:4]))
            if amin > 0 and ia / amin > 0.50:
                bad.append((ia / amin, a[4], b[4]))
    if bad:
        bad.sort(reverse=True)
        rec('FAIL', 'P', '%s：%d 处文字互相压住（最严重：%s 覆盖 %s 的 %.0f%%）'
            % (tag, len(bad), bad[0][1], bad[0][2], bad[0][0] * 100))
    else:
        rec('PASS', 'P', '%s：%d 块文字互不压盖（阈值：交叠不超过较小者的 50%%）'
            % (tag, len(words)))

    # 7) 彩色描边路径（箭头/连线）穿过文字
    pierced = []
    for d in page.get_drawings():
        col = d.get('color')
        if not col or d.get('type') not in ('s', 'fs'):
            continue
        rgb = col[:3]
        if max(rgb) - min(rgb) < 0.15:          # 黑/灰线条（正文规则线、黑色曲线）不算
            continue
        r = d['rect']
        if max(r.width, r.height) < 10.0:       # 太短的装饰线段不算
            continue
        wdt = max(d.get('width') or 0.0, 0.6)
        rr = (r.x0 - wdt, r.y0 - wdt, r.x1 + wdt, r.y1 + wdt)
        if rr[3] < ytop or rr[1] > ybot:
            continue
        for w in words:
            if inter_area(rr, w[:4]) > 3.0:
                pierced.append((w[4], rr))

    if pierced:
        rec('FAIL', 'P', '%s：%d 处彩色线条/箭头压住文字（例：%s 被 x=%.1f–%.1f pt 的线条穿过）'
            % (tag, len(pierced), pierced[0][0], pierced[0][1][0], pierced[0][1][2]))
    else:
        rec('PASS', 'P', '%s：没有彩色线条/箭头压住文字' % tag)

    rec('INFO', 'P', '%s：本页第 %d 页，区间 %.0f–%.0f pt，文字 %d 块（示意图规则 6/7）'
        % (tag, page.number + 1, ytop, ybot, len(words)))
    return True


def figure_bands(page, entries=None):
    """本页每个目标图形对应的纵向区间：上一个题注锚点之后 → 本题注锚点下方 40 pt。

    题目中的图是「图在上、题注在下」，所以区间取题注上方到上一个题注之间；
    若只按固定高度（如题注上方 300 pt）截取，高图的顶部节点会落在区间外，
    检查就漏掉了它们（v1.7.0 自查发现：图 9.2 有 3 个顶部节点未被覆盖）。
    """
    if entries is None:
        entries = ALL_FIGURES
    words = page.get_text('words')
    anchors = []
    for key, tag in entries:
        hit = [w for w in words if key in w[4]]
        if not hit:
            continue
        y = min(w[1] for w in hit)
        anchors.append((y, key, tag))
    anchors.sort()
    bands = []
    prev = 0.0
    for y, key, tag in anchors:
        bands.append((key, tag, max(0.0, prev + 10.0), min(page.rect.height, y + 40.0)))
        prev = y
    return bands



def _rect(bbox):
    """把 bbox（4 元组）转成与绘图 rect 同型的矩形对象。"""
    import pymupdf as _f
    return _f.Rect(bbox)


# ---------------------------------------------------------------- 图内文字

MIN_FIG_FONT = 6.5          # pt：图内标注的最小可读字号
OVERLAP_FRAC = 0.25         # 交叠面积 / 较小 span 面积 的告警阈值


def figure_regions(page):
    """本页的「图形区」矩形列表：排除页眉页脚带与整栏宽的盒子框。"""
    W = page.rect.width
    H = page.rect.height
    rects = []
    for d in page.get_drawings():
        r = d['rect']
        if r.y0 < 70.0 or r.y1 > H - 70.0:
            continue
        if r.width > 0.85 * W and r.height > 60.0:      # tcolorbox 之类的整栏框
            continue
        if r.width < 6.0 or r.height < 6.0:
            continue
        rects.append(r)
    return rects


def caption_bands(page):
    """本页所有「图 X.Y」题注行给出的 figure 纵向区间。

    不能直接用 `figure_bands()`：它只覆盖 ALL_FIGURES 登记过的图（9.1--9.3），
    其余 17 张图会被整页跳过——真重叠也就漏掉了（负对照实测：限定到 ALL_FIGURES 后
    旧 PDF 里图 2.1 的标签互压不再报出）。这里改为认页上**任何**一条图题注行，
    区间取「上一条题注下方 10 pt → 本题注下方 40 pt」，与 figure_bands 同口径。
    """
    anchors = []
    for blk in page.get_text('dict')['blocks']:
        for ln in blk.get('lines', []):
            txt = ''.join(sp['text'] for sp in ln['spans']).strip()
            if re.match(r'^图\s*\d+\.\d+', txt):
                anchors.append(ln['bbox'][1])
    anchors.sort()
    bands, prev = [], 0.0
    for y in anchors:
        bands.append((max(0.0, prev + 10.0), min(page.rect.height, y + 40.0)))
        prev = y
    return bands


def check_figure_text_overlap(doc):
    """扫描全书：figure float 内的文字块两两重叠、以及字号过小的标注。

    检查范围限定在**真正的 figure** 里：纵向区间由题注锚点给出（见 `caption_bands()`）。
    此前按「页上有 drawing 就算图区」扫描，会把显示公式也算进去——分式的横线是 drawing
    元素，而分子/分母的 span 矩形上下叠置、天然交叠，于是 p87 的 Vlasov 方程
    $\\partial f/\\partial x$ 被报成「图内文字重叠」。
    """
    over, small = [], []
    n_span = 0
    for pno, page in enumerate(doc, 1):
        bands = caption_bands(page)
        if not bands:
            continue
        regions = figure_regions(page)
        spans = []
        d = page.get_text('dict')
        for bi, blk in enumerate(d['blocks']):
            for ln in blk.get('lines', []):
                for sp in ln['spans']:
                    txt = sp['text'].strip()
                    r = _rect(sp['bbox'])
                    if len(txt) < 2:
                        continue
                    if r.y0 < 70.0 or r.y1 > page.rect.height - 70.0:
                        continue
                    if not any(yt <= r.y0 and r.y1 <= yb for yt, yb in bands):
                        continue
                    if not any(r.intersects(rg) for rg in regions):
                        continue
                    spans.append((r, txt, sp['size'], bi))
                    if sp['size'] < MIN_FIG_FONT:
                        small.append('p%d：%.1f pt「%s」' % (pno, sp['size'], txt[:16]))
        n_span += len(spans)
        for i in range(len(spans)):
            for j in range(i + 1, len(spans)):
                ra, ta, _, ba = spans[i]
                rb, tb, _, bb = spans[j]
                if ba == bb:
                    continue                     # 同一文本块内不判（数学 span 常轻微重叠）
                # 公式被拆成相邻 span：同基线水平紧邻，或分式的分子/分母（上下叠置、
                # 水平居中对齐）。两者都是纯符号片段时才豁免——中文标签的真实碰撞
                # 不会同时满足「纯符号」与「上下居中对齐」。
                _syms = set('abcdefghijklmnopqrstuvwxyz'
                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
                            '∂∇×·<>/=+-^_{}()[]|,.;:')
                _cya = (ra.y0 + ra.y1) / 2.0
                _cyb = (rb.y0 + rb.y1) / 2.0
                _cxa = (ra.x0 + ra.x1) / 2.0
                _cxb = (rb.x0 + rb.x1) / 2.0
                _dy = abs(_cya - _cyb)
                if set(ta) <= _syms and set(tb) <= _syms and (
                        _dy < 2.0 or (_dy < 22.0 and abs(_cxa - _cxb) < 10.0)):
                    continue      # 同一个公式的碎片，不是碰撞
                a = inter_area(ra, rb)
                if a <= 0:
                    continue
                m = min(area(ra), area(rb))
                if m > 0 and a / m > OVERLAP_FRAC:
                    over.append('p%d：交叠 %d%%「%s」×「%s」'
                                % (pno, int(100 * a / m), ta[:14], tb[:14]))
    if small:
        for s in small[:6]:
            rec('FAIL', 'P', '图内标注字号过小（< %.1f pt）：%s' % (MIN_FIG_FONT, s))
        if len(small) > 6:
            rec('FAIL', 'P', '图内标注字号过小：另有 %d 处' % (len(small) - 6))
    if over:
        for s in over[:6]:
            rec('FAIL', 'P', '图内文字重叠：%s' % s)
        if len(over) > 6:
            rec('FAIL', 'P', '图内文字重叠：另有 %d 处' % (len(over) - 6))
    if not small and not over:
        rec('PASS', 'P', '插图内文字：%d 个 span 无重叠、最小字号 ≥ %.1f pt'
            % (n_span, MIN_FIG_FONT))


# ---------------------------------------------------------------- 插图配色规范（静态）

# 与 tmp/palette_fix_all.py 使用**同一个正则**：三种写法（fill=blue!20 / [blue!20] / \color{blue}）
# 一次覆盖；`!black` 这类混色目标由前置否定 (?<![A-Za-z!]) 放行；
# 前置再加 `\\`：`\blacksquare` 这类命令名里出现的 black 不是颜色（v1.26.39 修）。
PALETTE_COLORS = ('figBlue', 'figVerm', 'figOrange', 'figGreen', 'figPurple', 'figSky', 'figGray', 'figGray', 'figInk', 'figOrange', 'figPurple', 'figGray', 'figPurple', 'figGreen', 'figGreen', 'figPurple')
PALETTE_PAT = re.compile('(?<![A-Za-z!\\\\])(magenta|orange|purple|yellow|violet|green|black|brown|olive|blue|cyan|gray|grey|teal|pink|red)(!\\d+(?:![a-z]+)?)?')
FIG_ENV = re.compile(r'\\begin\{figure\}(?:.|\n)*?\\end\{figure\}')


def check_figure_palette():
    """静态扫描：插图里不得再出现写死颜色（一律用统一色板 figXxx）。"""
    import glob as _glob
    src = os.path.join(ROOT, 'src')
    bad = []
    for path in sorted(_glob.glob(os.path.join(src, '*.tex'))):
        raw = open(path, encoding='utf-8').read()
        raw = re.sub(r'(?<!\\)%.*', '', raw)
        for m in FIG_ENV.finditer(raw):
            block = m.group(0)
            for wm in PALETTE_PAT.finditer(block):
                # 只认「颜色用法」：前面紧邻 [ = { , 之一（fill=blue!20 / [blue!20] / \color{blue} / draw=red）。
                # 否则是文件名或普通单词（如 ../figures/red_sprites_iss.jpg，v1.26.56 修）。
                prev = block[max(0, wm.start() - 1):wm.start()]
                if prev not in ('[', '=', '{', ','):
                    continue
                line = raw[:m.start() + wm.start()].count('\n') + 1
                bad.append('%s:%d 用了写死颜色 %s' % (os.path.basename(path), line, wm.group(0)))
    if bad:
        for b in bad[:6]:
            rec('FAIL', 'P', '插图配色规范：%s（应改用统一色板 figBlue/figVerm/…）' % b)
        if len(bad) > 6:
            rec('FAIL', 'P', '插图配色规范：另有 %d 处' % (len(bad) - 6))
    else:
        rec('PASS', 'P', '插图配色规范：无写死颜色（全部走统一色板）')


def check_photo_sources():
    """静态扫描：正文引用的实拍/机构图必须①文件存在②图环境里有 \\figsource 署名
    ③在 figures/CREDITS.md 里登记④单张不超过 1.5 MB。"""
    import glob as _glob
    src = os.path.join(ROOT, 'src')
    figdir = os.path.join(ROOT, 'figures')
    credits_path = os.path.join(figdir, 'CREDITS.md')
    credits = open(credits_path, encoding='utf-8').read() if os.path.exists(credits_path) else ''
    inc = re.compile(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}')
    missing, no_src, unlisted, heavy, total = [], [], [], [], 0
    for path in sorted(_glob.glob(os.path.join(src, '*.tex'))):
        raw = open(path, encoding='utf-8').read()
        raw_nc = re.sub(r'(?<!\\)%.*', '', raw)
        for m in inc.finditer(raw_nc):
            rel = m.group(1).strip()
            total += 1
            p = os.path.normpath(os.path.join(src, rel))
            base = os.path.basename(rel)
            line = raw_nc[:m.start()].count('\n') + 1
            if not os.path.exists(p):
                missing.append('%s:%d -> %s' % (os.path.basename(path), line, rel))
                continue
            sz = os.path.getsize(p)
            if sz > 1500 * 1024:
                heavy.append('%s（%.2f MB）' % (base, sz / 1048576.0))
            if base not in credits:
                unlisted.append(base)
        for fm in FIG_ENV.finditer(raw_nc):
            blk = fm.group(0)
            if '\\includegraphics' in blk and '\\figsource{' not in blk:
                line = raw_nc[:fm.start()].count('\n') + 1
                no_src.append('%s:%d' % (os.path.basename(path), line))
    if total == 0:
        rec('FAIL', 'P', '实拍图：正文里一处 \\includegraphics 都没有（配图被删？）')
        return
    if missing:
        for b in missing[:5]:
            rec('FAIL', 'P', '实拍图文件缺失：%s' % b)
    if no_src:
        for b in no_src[:5]:
            rec('FAIL', 'P', '图环境缺 \\figsource 署名：%s' % b)
    if unlisted:
        for b in sorted(set(unlisted))[:5]:
            rec('FAIL', 'P', '实拍图未在 figures/CREDITS.md 登记：%s' % b)
    if heavy:
        for b in heavy[:3]:
            rec('WARN', 'P', '实拍图偏大（建议 ≤1.5 MB）：%s' % b)
    if not (missing or no_src or unlisted):
        rec('PASS', 'P', '实拍图：%d 处引用均存在、均有署名、均已在 CREDITS.md 登记' % total)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdf', default=DEFAULT_PDF)
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.pdf):
        rec('FAIL', 'P', '找不到编译产物：%s（先运行 check/full_check.py 编译）' % args.pdf)
    else:
        try:
            import pymupdf as fitz          # 新版包名（PyMuPDF ≥1.24）
        except ImportError:
            try:
                import fitz                 # 旧版包名
            except ImportError:
                rec('WARN', 'P', '未安装 PyMuPDF，跳过图形几何检查（install: pip install PyMuPDF）')
                fitz = None
        if fitz is not None:
            check_figure_palette()
            check_photo_sources()
            doc = fitz.open(args.pdf)
            seen = set()
            box_keys = set(k for k, _ in FIGURES)
            for page in doc:
                for key, tag, ytop, ybot in figure_bands(page):
                    if key in box_keys:
                        ok = analyze(page, key, tag, ytop, ybot)
                    else:
                        ok = analyze_diagram(page, key, tag, ytop, ybot)
                    if ok:
                        seen.add(tag)
            check_figure_text_overlap(doc)
            doc.close()
            for _, tag in ALL_FIGURES:
                if tag not in seen:
                    rec('FAIL', 'P', '%s：在整本 PDF 中都没找到（图形消失或题注被改）' % tag)

    if args.json:
        print(json.dumps([{'level': a, 'cat': b, 'msg': c} for a, b, c in RESULTS],
                         ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》图形版式自动守卫')
    print('=' * 78)
    for lvl in ['FAIL', 'WARN', 'INFO', 'PASS']:
        rows = [r for r in RESULTS if r[0] == lvl]
        if not rows:
            continue
        print('\n--- %s (%d) ---' % (lvl, len(rows)))
        for _, cat, msg in rows:
            print('  [%s] %s' % (cat, msg))
    nf = sum(1 for r in RESULTS if r[0] == 'FAIL')
    nw = sum(1 for r in RESULTS if r[0] == 'WARN')
    print('\n' + '=' * 78)
    print('汇总：FAIL %d ／ WARN %d ／ 合计 %d 项' % (nf, nw, len(RESULTS)))
    print('=' * 78)
    return 1 if nf else 0


if __name__ == '__main__':
    sys.exit(main())
