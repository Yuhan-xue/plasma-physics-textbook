#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》图形版式自动守卫（量编译后 PDF 的几何）
================================================================

为什么要这个文件：图坏在哪可以一眼看出，但改一次图就可能重新长回来；把"看得出
坏了"写成可重跑的断言，重绘的人才有回执。所有判据都从 PDF 取值，不信任源码意图。

规则组（编号沿用 v1.7.0，阈值一字未改）
--------------------------------------
  P. 带节点框的图（流程图、分类树一类）：
      1) 填充节点框两两相交面积 > 1 pt²                  → FAIL（节点交叠）
      2) 文字与某框相交 >5%，但落在框内 <90% 且探出 >1 pt  → FAIL（标签被框边裁切）
      3) 一块文字同时落进两个框（面积和 >60%）            → FAIL（文字跨框）
      4) 节点框越出版心（A4，left 1.75 cm / right 1.95 cm）→ FAIL（超宽）
  Q. 没有节点框的示意图（瑞利--泰勒、曲线图一类）：
      6) 两块文字互相压住（交叠面积/较小者 > 50%）        → FAIL（文字叠字）
      7) 彩色描边路径（箭头、连线）穿过文字包围盒          → FAIL（线条穿字）
  另加两条通用项：图的实测外框越出版心（FAIL）；图宽 > 86% 版心、宽 > 95% 版心、
  高 > 12 cm（WARN，只提醒"这张图明显过大"，不算硬错）。

  规则 2 的 90%/1 pt 与规则 6 只比 word 级包围盒，都是 v1.7.0 定的：字形包围盒的
  ascent 会比节点填充框高出不到 1 pt，那是量测噪声不是裁切；而 PyMuPDF 会把
  $v_t^2$ 拆成 `v2` 与 `t` 两个 word，做数学上下标拆分会让全书误报。

定位方式（v1.27.0 起不再手写 3 张图的词表）
------------------------------------------
1. 扫 `src/*.tex`：每个 `\\begin{tikzpicture}` 是一张候选图，其后紧邻的
   `\\caption{...}` 与 `\\label{fig:...}` 是它的身份。封面 overlay（`remember picture`
   铺满整页做背景）不是插图，跳过；独立文档（自带 `\\documentclass`，如配套手册）
   不在本书 PDF 里，也跳过但在输出里列出。
2. 页码优先取 `src/plasma_physics_textbook_v2.aux` 的 `\\newlabel{fig:xxx}{{图号}{页}}`；
   aux 缺失就用题注里「最长的一段连续汉字」去 PDF 里搜题注行。为什么不是机械取前
   8 个汉字：PDF 文本层会把 `$I$-$V$`、`$n\\tau$` 这类数学从汉字串中间切断，拼出来
   的串在 PDF 里并不存在（实测朗缪尔探针那张图会抓不到）。题目里"40 行窗口"也不够：
   附录导航图的题注在 tikzpicture 结束后的第 62 行、DC 辉光图在第 78 行，故按
   "tikz 结束 + 1500 字符"取，遇到 `\\end{figure}`/段落命令提前收手。
3. 图区 = 题注上方那块「连续内容」：从题注行往上逐块吃，遇到正文段落行、上一个题注、
   tcolorbox 大边框、页眉页脚带，或与上一块的空白 > 40 pt 就停。按固定高度截
   （如"题注上方 300 pt"）会把学习路径图上方那张表格一起吞进图区，节点框规则就会
   把表格判成图。
4. 两种定位都失败 → WARN 并列出源码位置。**绝不静默通过**。

用法
----
    python -X utf8 check/figure_layout.py                      # 全书自画图体检表
    python -X utf8 check/figure_layout.py --only fig:glow      # 只查一张（可重复）
    python -X utf8 check/figure_layout.py --json               # 机读（含逐图实测尺寸）
    python -X utf8 check/figure_layout.py --pdf tmp/x.pdf --src src
    python -X utf8 check/figure_layout.py --no-legacy          # 关掉旧词表兼容检查

参数：
    --pdf PATH     编译产物（默认 pdf/plasma_physics_textbook_v2.pdf）
    --src DIR      书稿目录（默认 src），同目录下找 .aux
    --only SUBSTR  只检查 label 含该子串的图，可重复；给了它默认不再跑旧词表
    --json         输出 JSON（figures + 旧 results + summary）
    --legacy       强制跑旧词表检查（默认就跑，除了配合 --only）
    --no-legacy    关掉旧词表检查

退出码：0 = 无 FAIL；1 = 存在 FAIL（含 PDF 缺失）
"""

import argparse
import glob
import json
import os
import re
import sys
import unicodedata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PDF = os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf')
DEFAULT_SRC = os.path.join(ROOT, 'src')
AUX_NAME = 'plasma_physics_textbook_v2.aux'
MAIN_TEX = AUX_NAME[:-4] + '.tex'

# 版心：A4 宽 21 cm，geometry 设 left=1.75 cm、right=1.95 cm（v1.9.0 版面密度规格）
PT_CM = 72.0 / 2.54
PAGE_W = 21.0 * PT_CM
TEXT_L = 1.75 * PT_CM
TEXT_R = PAGE_W - 1.95 * PT_CM
TEXT_W = TEXT_R - TEXT_L

# 图区判定的常数：上下 52 pt 之外是页眉页脚；块与块之间超过 40 pt 的空白视为「图结束了」。
# 这几个数按本书实测的正文行距（18.8 pt）与浮动体间距定。
BAND_TOP = 52.0
BAND_BOTTOM = 52.0
GAP_MAX = 40.0

# 「明显过大」的提醒线（不算硬错，只进 WARN）
W_BIG = 0.86        # 宽占版心 86% 以上 → 偏大
W_HUGE = 0.95       # 宽占版心 95% 以上 → 明显过大
H_HUGE = 12.0       # 高 12 cm 以上 → 明显过高

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


# ---------------------------------------------------------------- 规则谓词
#
# 这里只放「判定」不产生文案：旧词表路径与新自动发现路径共用同一套阈值，
# 免得同一个物理语义在两处各写一遍、日后走样。

def is_nested(a, b):
    """一个框基本包住另一个框：那是「框里套框」的同一处图形，不是两个框打架。

    本书不少节点是「带底色的外框 + 白底小标签」两层叠着画的，若按「两框相交
    > 1 pt² 即交叠」判，这类图会整片误报（实测图 5.4／5.5 的嵌套框各占 731.9 pt²）。
    """
    smaller = min(area(a), area(b))
    return smaller > 0 and inter_area(a, b) / smaller > 0.90


def overlap_pairs(boxes):
    """规则 1：节点框两两相交面积 > 1 pt² 的配对（按严重度降序）。"""
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if is_nested(boxes[i], boxes[j]):
                continue
            a = inter_area(boxes[i], boxes[j])
            if a > 1.0:
                bad.append((a, i, j))
    bad.sort(reverse=True)
    return bad


def _hit_ratios(boxes, wb):
    """一个字/词落在各节点框内的面积占比（只留 >5% 的）。"""
    wa = area(wb)
    if wa <= 0:
        return []
    hits = []
    for b in boxes:
        ia = inter_area(wb, b)
        if ia / wa > 0.05:
            hits.append((ia / wa, b))
    return hits


def clip_hits(boxes, words, ytop, ybot):
    """规则 2：文字压到框边（落在框内不足 90%，且真的探出框外 >1 pt）。"""
    out = []
    for w in words:
        wb = (w[0], w[1], w[2], w[3])
        if wb[3] < ytop or wb[1] > ybot:
            continue
        hits = _hit_ratios(boxes, wb)
        if not hits:
            continue
        ratio, bb = max(hits, key=lambda t: t[0])
        over = max(bb[0] - wb[0], wb[2] - bb[2], bb[1] - wb[1], wb[3] - bb[3], 0.0)
        if ratio < 0.90 and over > 1.0:
            out.append((w[4], ratio, over))
    return out


def cross_hits(boxes, words, ytop, ybot):
    """规则 3：一块文字同时落进两个框（面积和 > 60%）。

    先把「框套框」的两层并成一处，再数到底跨了几个框：否则一个坐在两层框里的
    普通标签会被算成「覆盖 200%」。
    """
    out = []
    for w in words:
        wb = (w[0], w[1], w[2], w[3])
        if wb[3] < ytop or wb[1] > ybot:
            continue
        hits = _hit_ratios(boxes, wb)
        merged = [(r, b) for r, b in hits
                  if not any(b2 is not b and area(b2) > area(b) and is_nested(b, b2)
                             for _, b2 in hits)]
        if sum(r for r, _ in merged) > 0.60 and len(merged) > 1:
            out.append((w[4], sum(r for r, _ in merged)))
    return out


def word_overlaps(words):
    """规则 6：两块文字互相压住（交叠面积/较小者 > 50%）。"""
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
    bad.sort(reverse=True)
    return bad


def _seg_hits_rect(p0, p1, r):
    """线段与矩形是否相交（Liang--Barsky 裁剪）。"""
    x0, y0, x1, y1 = r
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, p0[0] - x0), (dx, x1 - p0[0]), (-dy, p0[1] - y0), (dy, y1 - p0[1])):
        if abs(p) < 1e-12:
            if q < 0:
                return False
            continue
        t = q / p
        if p < 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return False
    return True


def _stroke_segments(d, samples=8):
    """把一段 drawing 拆成折线（曲线按 8 段线性化），用于「线真的穿过字」判定。

    为什么不能只用 d['rect']：一条平滑曲线的包围盒等于整个绘图区，落在绘图区里的
    任何标注都会被误判成「被曲线穿过」。按线段判才分得清"从旁边过"和"从字上过"。
    """
    segs = []
    for it in d.get('items') or []:
        op = it[0]
        if op == 'l':
            segs.append(((it[1].x, it[1].y), (it[2].x, it[2].y)))
        elif op == 're':
            r = it[1]
            p = [(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1)]
            segs += [(p[k], p[(k + 1) % 4]) for k in range(4)]
        elif op == 'qu':
            q = it[1]
            p = [(q.ul.x, q.ul.y), (q.ur.x, q.ur.y), (q.lr.x, q.lr.y), (q.ll.x, q.ll.y)]
            segs += [(p[k], p[(k + 1) % 4]) for k in range(4)]
        elif op == 'c':
            pts = [(it[k].x, it[k].y) for k in (1, 2, 3, 4)]
            prev = pts[0]
            for k in range(1, samples + 1):
                t = float(k) / samples
                mt = 1.0 - t
                x = (mt ** 3 * pts[0][0] + 3 * mt * mt * t * pts[1][0]
                     + 3 * mt * t * t * pts[2][0] + t ** 3 * pts[3][0])
                y = (mt ** 3 * pts[0][1] + 3 * mt * mt * t * pts[1][1]
                     + 3 * mt * t * t * pts[2][1] + t ** 3 * pts[3][1])
                segs.append((prev, (x, y)))
                prev = (x, y)
    return segs


def _occluders(page):
    """「白底标注底板」：近白填充矩形，按绘制顺序返回 [(序号, rect)]。

    为什么需要它：TikZ 里 `node[fill=white]` 会给标签垫一块白底，让曲线从标签背后
    穿过去时被挡住——这是**有意为之**的排版，视觉上完全正常。只看几何会把这类图整片
    误报成「线条穿字」。绘制顺序（content stream 顺序）告诉我们谁盖住谁：底板排在线
    之后，说明它确实压在线上。
    """
    out = []
    for k, d in enumerate(page.get_drawings()):
        f = d.get('fill')
        if f is None or min(f[:3]) < 0.93:     # 只要近白
            continue
        r = d['rect']
        if r.width < 4.0 or r.height < 4.0:
            continue
        out.append((k, (r.x0, r.y0, r.x1, r.y1)))
    return out


def pierce_hits(page, words, ytop, ybot, precise=False):
    """规则 7：彩色描边路径（箭头、连线）穿过文字包围盒。

    precise=True 时按真实线段判定，并要求标签没有白底托着（自动发现路径用）；
    False 时按包围盒判定（旧词表路径保持 v1.7.0 的行为，供
    check/fixtures/gate_figure_layout.py 回归）。
    """
    pierced = []
    draws = page.get_drawings()
    occ = _occluders(page) if precise else []
    for k, d in enumerate(draws):
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
        segs = _stroke_segments(d) if precise else None
        for w in words:
            wb = (w[0], w[1], w[2], w[3])
            if segs:
                wr = (wb[0] - wdt, wb[1] - wdt, wb[2] + wdt, wb[3] + wdt)
                hit = any(_seg_hits_rect(p, q, wr) for p, q in segs)
            else:
                hit = inter_area(rr, wb) > 3.0
            if not hit:
                continue
            wa = area(wb)
            if occ and wa > 0:
                # 线上盖着一块白底板 → 标签本来就把线挡住了，不是缺陷
                if any(kk > k and inter_area(ob, wb) / wa >= 0.85 for kk, ob in occ):
                    continue
            pierced.append((w[4], rr))
    return pierced


# ---------------------------------------------------------------- 题注 / 图区

CJK_RE = re.compile(r'[\u4e00-\u9fff]')
CAP_ANCHOR_RE = re.compile(r'^\s*图\s*([0-9A-Z]+(?:\.[0-9]+)?)(?![\d.])')


def _brace_block(s, pos):
    """s[pos] 是 '{'：返回配对括号内的内容与结束位置。"""
    if pos >= len(s) or s[pos] != '{':
        return None, pos
    depth, i, out = 0, pos, []
    while i < len(s):
        ch = s[i]
        if ch == '{':
            depth += 1
            if depth == 1:
                i += 1
                continue
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return ''.join(out), i + 1
        out.append(ch)
        i += 1
    return None, pos


def _top_groups(s):
    """把 '{{a}{b}{c}}' 的内层切成 ['a','b','c']。"""
    groups, i = [], 0
    while i < len(s):
        if s[i] == '{':
            g, j = _brace_block(s, i)
            if g is None:
                break
            groups.append(g)
            i = j
        else:
            i += 1
    return groups


def parse_aux_labels(path):
    """`\\newlabel{fig:xxx}{{图号}{页码}{题注}{锚点}{}}` → {label: [图号, 页码, 题注, ...]}

    不能按 `{{...}{...}}` 硬切：本书第三段是完整题注（带 `\\textbf{}`、可跨行），
    必须走配对括号解析。
    """
    if not os.path.exists(path):
        return {}
    raw = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for m in re.finditer(r'\\newlabel\{(fig:[^}]*)\}\{', raw):
        inner, _ = _brace_block(raw, m.end() - 1)
        if inner is None:
            continue
        out[m.group(1)] = _top_groups(inner)
    return out


def _strip_tex(s):
    """题注源码压成可搜的纯文字：去数学、去命令、去括号。"""
    s = re.sub(r'\$[^$]*\$', ' ', s or '')
    s = re.sub(r'\\(?:textbf|emph|textit|text|mbox)\{', '{', s)
    s = re.sub(r'\\[a-zA-Z]+', ' ', s)
    s = s.replace('{', '').replace('}', '').replace('~', ' ').replace('--', '-')
    return s


def search_key(caption):
    """题注里最长的一段连续汉字（≤8 字），用作 PDF 搜索键。"""
    t = _strip_tex(caption)
    best, cur = '', ''
    for ch in t:
        if CJK_RE.match(ch):
            cur += ch
            best = cur if len(cur) > len(best) else best
        else:
            cur = ''
    if len(best) >= 4:
        return best[:8]
    m = re.search(r'[A-Za-z]{5,}', t)
    return m.group(0) if m else (best or t.strip()[:12])


def _norm(s):
    return re.sub(r'\s+', '', s or '')


def _caption_in(window):
    """窗口里第一个 `\\caption` / `\\caption[短]{长}` 的长题注。"""
    i = window.find('\\caption')
    if i < 0:
        return None
    j = i + len('\\caption')
    while j < len(window) and window[j] in ' \t\r\n':
        j += 1
    if j < len(window) and window[j] == '[':
        k = window.find(']', j)
        if k > 0:
            j = k + 1
        while j < len(window) and window[j] in ' \t\r\n':
            j += 1
    if j >= len(window) or window[j] != '{':
        return None
    text, _ = _brace_block(window, j)
    return text


def discover_figures(src_dir=DEFAULT_SRC):
    """扫 `src/*.tex`，返回全书自画图的身份表（每张一个 dict）。

    窗口取「`\\end{tikzpicture}` 之后 1500 字符」，并在 `\\end{figure}` / `\\begin{figure}`
    / `\\begin{table}` / `\\section` 处提前收手：既够到 tikz 很长、题注在 60 行之后
    的图（附录导航图、DC 辉光图），又不会把下一张图的题注错认成自己的。
    """
    figs = []
    for path in sorted(glob.glob(os.path.join(src_dir, '*.tex'))):
        raw = open(path, encoding='utf-8', errors='replace').read()
        if '\\begin{tikzpicture}' not in raw:
            continue
        base = os.path.basename(path)
        standalone = ('\\documentclass' in raw) and (base != MAIN_TEX)
        for m in re.finditer(r'\\begin\{tikzpicture\}(\[[^\]]*\])?', raw):
            opts = (m.group(1) or '').strip()
            end = raw.find('\\end{tikzpicture}', m.end())
            body_end = end if end >= 0 else m.end()
            window = raw[body_end:body_end + 1500]
            cut = len(window)
            for stop in ('\\end{figure}', '\\begin{figure}', '\\begin{table}',
                         '\\end{table}', '\\begin{tikzpicture}', '\\section{',
                         '\\subsection{'):
                k = window.find(stop)
                if 0 <= k < cut:
                    cut = k
            window = window[:cut]
            fig = {
                'label': None,
                'caption': _caption_in(window),
                'file': base,
                'line': raw[:m.start()].count('\n') + 1,
                'opts': opts,
                'skip': None,
            }
            gm = re.search(r'\\label\{(fig:[^}]*)\}', window)
            if gm:
                fig['label'] = gm.group(1)
            if 'overlay' in opts and ('remember picture' in opts or 'current page' in opts):
                fig['skip'] = '封面 overlay（remember picture 铺满整页，不是插图）'
            elif standalone:
                fig['skip'] = '独立文档（自带 \\documentclass），不在本书 PDF 内'
            figs.append(fig)
    return figs


def page_lines(page):
    """页面文本行：[(x0, y0, x1, y1, 字号, 文本)]，已剔除页眉页脚带。"""
    H = page.rect.height
    out = []
    for blk in page.get_text('dict')['blocks']:
        for ln in blk.get('lines', []):
            txt = ''.join(sp['text'] for sp in ln['spans']).strip()
            if not txt:
                continue
            bb = ln['bbox']
            if bb[1] < BAND_TOP or bb[3] > H - BAND_BOTTOM:
                continue
            out.append((bb[0], bb[1], bb[2], bb[3],
                        max((sp['size'] for sp in ln['spans']), default=0.0), txt))
    return out


def caption_entries(page):
    """本页题注：把「图 N.N」锚点行与右侧（同一视觉行）/下方（折行）的续行拼成一条。

    PDF 文本层把「图 9.4」和题注正文劈成两个 line 对象（同一 y），只取锚点行会
    丢掉整条题注的字。
    """
    lines = sorted(page_lines(page))
    out = []
    for (x0, y0, x1, y1, sz, txt) in lines:
        m = CAP_ANCHOR_RE.match(txt)
        if not m:
            continue
        num = m.group(1)
        bb = [x0, y0, x1, y1]
        text = [txt]
        # 同一视觉行右半段
        for (ax0, ay0, ax1, ay1, asz, atxt) in lines:
            if ax0 >= x1 - 1.0 and abs(ay0 - y0) <= 3.0:
                bb[2] = max(bb[2], ax1)
                bb[3] = max(bb[3], ay1)
                text.append(atxt)
        # 折行续行：行距紧、字号一致。不要求左端贴版心——题注最后一行常居中缩进
        # （实测图 3.5 的第三行从 x=117 起），按左端卡会把最后一行漏在题注外，
        # 那一行随后就被上方那张图当成自己的内容吃掉。
        prev_bottom = bb[3]
        for (ax0, ay0, ax1, ay1, asz, atxt) in lines:
            if ay0 < prev_bottom - 1.0:
                continue
            if ax0 < x0 - 3.0:
                continue
            if ay0 - prev_bottom > 10.0:
                continue
            if abs(asz - sz) > 1.6:
                continue
            if CAP_ANCHOR_RE.match(atxt):
                break
            bb[2] = max(bb[2], ax1)
            bb[3] = max(bb[3], ay1)
            text.append(atxt)
            prev_bottom = ay1
        out.append({'page': page.number + 1, 'num': num,
                    'x0': bb[0], 'y0': bb[1], 'x1': bb[2], 'y1': bb[3],
                    'size': sz, 'text': ''.join(text)})
    return out


def build_caption_index(doc):
    """一次扫全书，建「图号 → 题注条目」索引（避免逐图重扫 615 页）。"""
    by_num, every = {}, []
    for page in doc:
        for e in caption_entries(page):
            by_num.setdefault(e['num'], []).append(e)
            every.append(e)
    return by_num, every


def _page_items(page):
    """本页所有可参与「图区」判定的内容块（文字行 + 绘图矩形）。"""
    H = page.rect.height
    items = []
    for (x0, y0, x1, y1, sz, txt) in page_lines(page):
        items.append({'kind': 'text', 'bb': (x0, y0, x1, y1), 'size': sz, 'text': txt})
    page_area = page.rect.width * H
    for d in page.get_drawings():
        r = d['rect']
        bb = (r.x0, r.y0, r.x1, r.y1)
        if bb[3] < BAND_TOP or bb[1] > H - BAND_BOTTOM:
            continue
        if area(bb) > 0.60 * page_area:                 # 整页底纹
            continue
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        if w > 0.95 * page.rect.width and h > 0.5 * H:  # 整页色带
            continue
        items.append({'kind': 'draw', 'bb': bb, 'w': w, 'h': h})
    return items


def _is_caption_line(it):
    return it['kind'] == 'text' and bool(CAP_ANCHOR_RE.match(it['text']))


def _is_body_line(it):
    """正文段落行：从版心左边界起、且够宽（或字号达到正文字号）。"""
    if it['kind'] != 'text':
        return False
    w = it['bb'][2] - it['bb'][0]
    if it['bb'][0] > TEXT_L + 3.0:
        return False
    return w >= 0.88 * TEXT_W or (it['size'] >= 11.6 and w >= 0.45 * TEXT_W)


def _is_box_frame(it):
    """tcolorbox 一类整栏大边框：宽到几乎通栏、又够高。"""
    if it['kind'] != 'draw':
        return False
    w = it['bb'][2] - it['bb'][0]
    return w >= 0.90 * TEXT_W and (it['bb'][3] - it['bb'][1]) >= 60.0


def figure_region(page, cap_top, cap_excl=(), gap_max=None):
    """题注上方那块「连续内容」的包围盒 → ((x0,y0,x1,y1) 或 None, 块数)。

    两步走，而不是「按 bottom 排序逐块往上吃」：
      1) 从题注顶往上找最近的、位于图区之上的块，声明 [它的顶, 图区底] 这段纵向区间
         属于图；
      2) 把**纵向落在这段区间里**的块全部收进来。
    只做第 1 步会漏东西：图里常有一条贯穿全高的竖线/箭头，它一被吃进来就会把
    「上沿」直接推到页面顶部，夹在中间的文字（节点标签、图例）就永远轮不到——
    实测图 3.5 只量到 0.19 cm 宽、学习路径图只量到 45% 版心，都是这个原因。

    `cap_excl` 是本页其它题注的包围盒：撞上它就说明已经越过上一张图，收手
    （否则两张图挤在一页时，上面那张会被下面那张一起吃掉）。

    停止条件：与上一块的空白 > gap_max、正文段落行、其它题注、tcolorbox 大边框、
    页眉页脚带。通栏细线（页内分隔线）既不算图、也不当边界。
    """
    if gap_max is None:
        gap_max = GAP_MAX
    items = _page_items(page)
    ytop = cap_top
    ybot = cap_top
    used = set()
    blocked = False

    def in_excl(it):
        bb = it['bb']
        for (x0, y0, x1, y1) in cap_excl:
            if bb[1] >= y0 - 2.0 and bb[3] <= y1 + 2.0:
                return True
        return False

    while not blocked:
        cand = [(i, it) for i, it in enumerate(items)
                if i not in used and it['bb'][3] <= ytop + 2.5]
        if not cand:
            break
        i, it = max(cand, key=lambda t: t[1]['bb'][3])
        if ytop - it['bb'][3] > gap_max:
            break
        if _is_caption_line(it) or _is_body_line(it) or _is_box_frame(it) or in_excl(it):
            break
        if (it['kind'] == 'draw' and (it['bb'][2] - it['bb'][0]) >= 0.90 * TEXT_W
                and (it['bb'][3] - it['bb'][1]) <= 4.0):
            used.add(i)
            continue                        # 通栏细线（页内分隔线）不算图
        used.add(i)
        ytop = min(ytop, it['bb'][1])
        for j, jt in enumerate(items):
            if j in used:
                continue
            bb = jt['bb']
            if bb[1] < ytop - 1.0 or bb[3] > ybot + 1.0:
                continue                    # 不在刚声明的纵向区间里
            if _is_caption_line(jt) or in_excl(jt):
                # 声明越过了别的题注：上沿收到题注下方，此后停止外扩
                ytop = max(ytop, bb[3] + 1.0)
                blocked = True
                continue
            if _is_body_line(jt) or _is_box_frame(jt):
                continue
            used.add(j)
    if not used:
        return None, 0
    bbs = [items[i]['bb'] for i in used]
    box = (min(b[0] for b in bbs),
           max(ytop, min(b[1] for b in bbs)),
           max(b[2] for b in bbs),
           min(ybot + 1.0, max(b[3] for b in bbs)))
    return box, len(used)


def region_for_caption(page, ent):
    """给一条题注算图区：自动把本页其它题注的范围当成禁区。"""
    if ent is None:
        return None, 0
    others = [(c['x0'], c['y0'], c['x1'], c['y1'])
              for c in caption_entries(page) if abs(c['y0'] - ent['y0']) > 3.0]
    return figure_region(page, ent['y0'], cap_excl=others)


def locate_figure(fig, by_num, every, aux):
    """定位图所在页 → (页码, 题注条目, 定位方式) 或 (None, None, 失败原因)。"""
    lab = fig.get('label')
    num, page = None, None
    g = aux.get(lab) if lab else None
    if g:
        if g[0].strip():
            num = g[0].strip()
        if len(g) > 1 and g[1].strip().isdigit():
            page = int(g[1])
    if num and num in by_num:
        cands = by_num[num]
        if page is not None:
            near = [c for c in cands if abs(c['page'] - page) <= 2]
            if near:
                return near[0]['page'], near[0], 'aux 图号 %s（aux 记第 %d 页）' % (num, page)
        if cands:
            return cands[0]['page'], cands[0], 'aux 图号 %s（页码有漂移）' % num
    key = search_key(fig.get('caption') or '')
    if key:
        hits = [e for e in every if key in _norm(e['text'])]
        if page is not None and hits:
            hits.sort(key=lambda e: abs(e['page'] - page))
        if hits:
            return hits[0]['page'], hits[0], '题注关键词「%s」' % key
    if not lab:
        return None, None, '源码里既没有 \\label{fig:...} 也没有题注 → 无法在 PDF 定位'
    if not num:
        return None, None, 'aux 里没有该 label 的页码，题注关键词也没命中'
    return None, None, 'aux 给出图号 %s，但 PDF 全书找不到这条题注' % num


# ---------------------------------------------------------------- 单图判定

def analyze_figure(page, fig, ent, method, words=None):
    """对一张图跑 P/Q 两组规则，返回逐图记录（含实测尺寸）。

    有 ≥3 个节点框 → 走 P 组；否则说明它是示意图 → 走 Q 组。两组共用的「越出版心」
    与「明显过大」对所有图都查。
    """
    row = {
        'label': fig.get('label'), 'file': fig['file'], 'line': fig['line'],
        'caption': (fig.get('caption') or '').strip(),
        'number': (ent or {}).get('num'), 'page': page.number + 1,
        'locate': method, 'issues': [],
        'w_pt': None, 'h_pt': None, 'w_cm': None, 'h_cm': None, 'w_frac': None,
    }
    cap_top = ent['y0'] if ent else page.rect.height
    box, n_items = (region_for_caption(page, ent) if ent is not None
                    else figure_region(page, cap_top))
    if box is None:
        row['issues'].append(('WARN', 'G', '题注上方找不到任何图形内容（图可能已被删）'))
        return row
    w_pt = box[2] - box[0]
    h_pt = box[3] - box[1]
    row.update({'w_pt': round(w_pt, 1), 'h_pt': round(h_pt, 1),
                'w_cm': round(w_pt / PT_CM, 2), 'h_cm': round(h_pt / PT_CM, 2),
                'w_frac': round(w_pt / TEXT_W, 3),
                'box': [round(v, 1) for v in box]})
    if words is None:
        words = page.get_text('words')
    ytop, ybot = box[1], box[3]
    boxes = node_boxes(page, ytop - 1.0, ybot + 1.0)

    # 通用项：图的实测外框越出版心
    over_l = TEXT_L - box[0]
    over_r = box[2] - TEXT_R
    if max(over_l, over_r) > 2.0:
        row['issues'].append(('FAIL', 'P', '实测外框越出版心 %.1f pt（左 %.1f / 右 %.1f）'
                              % (max(over_l, over_r), max(over_l, 0.0), max(over_r, 0.0))))

    if len(boxes) >= 3:
        bad = overlap_pairs(boxes)
        if bad:
            row['issues'].append(('FAIL', 'P', '%d 对节点框交叠（最严重 %.1f pt²）'
                                  % (len(bad), bad[0][0])))
        clipped = clip_hits(boxes, words, ytop, ybot)
        if clipped:
            row['issues'].append(('FAIL', 'P', '%d 处文字被节点框边裁切（例：%s 落在框内 %.0f%%，探出 %.1f pt）'
                                  % (len(clipped), clipped[0][0], clipped[0][1] * 100, clipped[0][2])))
        crossing = cross_hits(boxes, words, ytop, ybot)
        if crossing:
            row['issues'].append(('FAIL', 'P', '%d 处文字同时压住两个节点框（例：%s 覆盖 %.0f%%）'
                                  % (len(crossing), crossing[0][0], crossing[0][1] * 100)))
        over = [b for b in boxes if b[0] < TEXT_L - 1.0 or b[2] > TEXT_R + 1.0]
        if over:
            worst = max(over, key=lambda b: max(TEXT_L - b[0], b[2] - TEXT_R))
            row['issues'].append(('FAIL', 'P', '%d 个节点框越出版心（最远 %.1f pt）'
                                  % (len(over), max(TEXT_L - worst[0], worst[2] - TEXT_R))))
        row['n_boxes'] = len(boxes)
        row['rules'] = 'P（节点框 %d 个）' % len(boxes)
    else:
        band_words = [w for w in words if w[1] >= ytop - 1.0 and w[3] <= ybot + 1.0]
        if not band_words:
            row['issues'].append(('WARN', 'Q', '图区里没有任何文字（只有图形），Q 组规则无法判定'))
        else:
            ov = word_overlaps(band_words)
            if ov:
                row['issues'].append(('FAIL', 'P', '%d 处文字互相压住（最严重：%s 覆盖 %s 的 %.0f%%）'
                                      % (len(ov), ov[0][1], ov[0][2], ov[0][0] * 100)))
            pr = pierce_hits(page, band_words, ytop, ybot, precise=True)
            if pr:
                row['issues'].append(('FAIL', 'P', '%d 处彩色线条/箭头压住文字（例：%s 被 x=%.1f–%.1f pt 的线条穿过）'
                                      % (len(pr), pr[0][0], pr[0][1][0], pr[0][1][2])))
        row['n_boxes'] = len(boxes)
        row['rules'] = 'Q（无节点框）'

    # 明显过大 / 挤压题注：都只进 WARN
    if row['w_frac'] > W_HUGE:
        row['issues'].append(('WARN', 'S', '图宽占版心 %.0f%%（> %.0f%%，明显过大）'
                              % (row['w_frac'] * 100, W_HUGE * 100)))
    elif row['w_frac'] > W_BIG:
        row['issues'].append(('WARN', 'S', '图宽占版心 %.0f%%（> %.0f%%，偏大）'
                              % (row['w_frac'] * 100, W_BIG * 100)))
    if row['h_cm'] > H_HUGE:
        row['issues'].append(('WARN', 'S', '图高 %.1f cm（> %.0f cm，明显过高）'
                              % (row['h_cm'], H_HUGE)))
    if ent is not None and box[3] > ent['y0'] - 4.0:
        row['issues'].append(('WARN', 'S', '图的下沿距题注顶只有 %.1f pt（< 4 pt，挤压题注）'
                              % (ent['y0'] - box[3])))
    return row


# ---------------------------------------------------------------- 旧词表路径（兼容）

def analyze(page, key, tag, ytop, ybot):
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

    bad = overlap_pairs(boxes)
    if bad:
        rec('FAIL', 'P', '%s：%d 对节点框交叠（最严重 %.1f pt²）'
            % (tag, len(bad), bad[0][0]))
    else:
        rec('PASS', 'P', '%s：%d 个节点框两两不交叠（阈值 1 pt²）' % (tag, len(boxes)))

    words = page.get_text('words')
    clipped = clip_hits(boxes, words, ytop, ybot)
    crossing = cross_hits(boxes, words, ytop, ybot)
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

    bad = word_overlaps(words)
    if bad:
        rec('FAIL', 'P', '%s：%d 处文字互相压住（最严重：%s 覆盖 %s 的 %.0f%%）'
            % (tag, len(bad), bad[0][1], bad[0][2], bad[0][0] * 100))
    else:
        rec('PASS', 'P', '%s：%d 块文字互不压盖（阈值：交叠不超过较小者的 50%%）'
            % (tag, len(words)))

    pierced = pierce_hits(page, words, ytop, ybot, precise=False)
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

    不能直接用 `figure_bands()`：它只覆盖 ALL_FIGURES 登记过的图（第 9 章那 3 张），
    其余图会被整页跳过——真重叠也就漏掉了（负对照实测：限定到 ALL_FIGURES 后
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
PALETTE_PAT = re.compile('(?<![A-Za-z!\\\\])(magenta|orange|purple|yellow|violet|green|black|brown|olive|blue|cyan|gray|grey|teal|pink|red)(!\\d+(?:![a-z]+)?)?')
FIG_ENV = re.compile(r'\\begin\{figure\}(?:.|\n)*?\\end\{figure\}')


def check_figure_palette(src_dir=DEFAULT_SRC):
    """静态扫描：插图里不得再出现写死颜色（一律用统一色板 figXxx）。"""
    bad = []
    for path in sorted(glob.glob(os.path.join(src_dir, '*.tex'))):
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


def check_photo_sources(src_dir=DEFAULT_SRC):
    """静态扫描：正文引用的实拍/机构图必须①文件存在②图环境里有 \\figsource 署名
    ③在 figures/CREDITS.md 里登记④单张不超过 1.5 MB。"""
    figdir = os.path.join(ROOT, 'figures')
    credits_path = os.path.join(figdir, 'CREDITS.md')
    credits = open(credits_path, encoding='utf-8').read() if os.path.exists(credits_path) else ''
    inc = re.compile(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}')
    missing, no_src, unlisted, heavy, total = [], [], [], [], 0
    for path in sorted(glob.glob(os.path.join(src_dir, '*.tex'))):
        raw = open(path, encoding='utf-8').read()
        raw_nc = re.sub(r'(?<!\\)%.*', '', raw)
        for m in inc.finditer(raw_nc):
            rel = m.group(1).strip()
            total += 1
            p = os.path.normpath(os.path.join(src_dir, rel))
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


# ---------------------------------------------------------------- 输出

def _disp(s):
    """字符串的终端显示宽度（CJK 按 2 列）。"""
    return sum(2 if unicodedata.east_asian_width(ch) in ('W', 'F') else 1 for ch in s)


def _pad(s, n):
    return s + ' ' * max(0, n - _disp(s))


def _cut(s, n):
    out, w = '', 0
    for ch in s:
        cw = 2 if unicodedata.east_asian_width(ch) in ('W', 'F') else 1
        if w + cw > n:
            break
        out += ch
        w += cw
    return out


VERDICT_ORDER = {'FAIL': 0, 'WARN': 1, 'PASS': 2}


def _verdict(row):
    levels = [i[0] for i in row['issues']]
    if 'FAIL' in levels:
        return 'FAIL'
    if 'WARN' in levels:
        return 'WARN'
    return 'PASS'


def _summary_of(row):
    issues = [i for i in row['issues'] if i[0] != 'PASS']
    if not issues:
        return '几何规则全过'
    return issues[0][2]


def print_table(rows):
    head = ['图号', 'label', '页', '宽cm', '宽/版心', '高cm', '结论', '问题']
    body = []
    for r in rows:
        body.append([
            (r.get('number') or '—'),
            r.get('label') or '(%s:L%d)' % (r['file'], r['line']),
            str(r.get('page') or '—'),
            '%.2f' % r['w_cm'] if r.get('w_cm') is not None else '—',
            '%.0f%%' % (r['w_frac'] * 100) if r.get('w_frac') is not None else '—',
            '%.2f' % r['h_cm'] if r.get('h_cm') is not None else '—',
            _verdict(r), _cut(_summary_of(r), 46),
        ])
    widths = [min(46, max(_disp(h), max((_disp(b[i]) for b in body), default=0)))
              for i, h in enumerate(head)]
    print('  ' + '  '.join(_pad(h, widths[i]) for i, h in enumerate(head)))
    print('  ' + '  '.join('-' * widths[i] for i in range(len(head))))
    for b in body:
        print('  ' + '  '.join(_pad(_cut(c, widths[i]), widths[i]) for i, c in enumerate(b)))


def main():
    ap = argparse.ArgumentParser(description='图形版式自动守卫（量 PDF 几何）')
    ap.add_argument('--pdf', default=DEFAULT_PDF)
    ap.add_argument('--src', default=DEFAULT_SRC)
    ap.add_argument('--only', action='append', default=None, metavar='SUBSTR')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--legacy', action='store_true', help='强制运行旧词表检查')
    ap.add_argument('--no-legacy', action='store_true', help='关闭旧词表检查')
    args = ap.parse_args()

    run_legacy = not args.no_legacy and (args.legacy or not args.only)

    rows = []
    fitz = None
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

    figs = discover_figures(args.src)
    skipped = [f for f in figs if f.get('skip')]
    active = [f for f in figs if not f.get('skip')]
    if args.only:
        active = [f for f in active
                  if any(s in (f.get('label') or '') for s in args.only)]

    if fitz is not None:
        check_figure_palette(args.src)
        check_photo_sources(args.src)
        doc = fitz.open(args.pdf)
        aux = parse_aux_labels(os.path.join(args.src, AUX_NAME))
        by_num, every = build_caption_index(doc)
        for f in active:
            page_no, ent, method = locate_figure(f, by_num, every, aux)
            if page_no is None:
                rows.append({'label': f.get('label'), 'file': f['file'], 'line': f['line'],
                             'caption': (f.get('caption') or '').strip(),
                             'number': None, 'page': None, 'locate': method,
                             'issues': [('WARN', 'G', '找不到页：%s' % method)],
                             'w_pt': None, 'h_pt': None, 'w_cm': None, 'h_cm': None,
                             'w_frac': None})
                continue
            rows.append(analyze_figure(doc[page_no - 1], f, ent, method))

        for r in rows:
            if r.get('page') is None:
                rec('WARN', 'G', '%s（%s:%d）：%s'
                    % (r.get('label') or '(无 label)', r['file'], r['line'], r['locate']))

        if run_legacy:
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
            for _, tag in ALL_FIGURES:
                if tag not in seen:
                    rec('FAIL', 'P', '%s：在整本 PDF 中都没找到（图形消失或题注被改）' % tag)
        check_figure_text_overlap(doc)
        doc.close()

    for r in rows:
        for lv, rule, msg in r['issues']:
            if lv in ('FAIL', 'WARN'):
                rec(lv, rule, '%s（图 %s，p%s）：%s'
                    % (r.get('label') or '(无 label)', r.get('number') or '?',
                       r.get('page') or '?', msg))

    rows.sort(key=lambda r: (VERDICT_ORDER.get(_verdict(r), 9),
                             r.get('page') or 9999, r.get('label') or ''))
    n_fail = sum(1 for r in rows if _verdict(r) == 'FAIL')
    n_warn = sum(1 for r in rows if _verdict(r) == 'WARN')
    n_pass = sum(1 for r in rows if _verdict(r) == 'PASS')
    big = [r for r in rows if (r.get('w_frac') or 0) > W_BIG or (r.get('h_cm') or 0) > H_HUGE]

    if args.json:
        print(json.dumps({
            'pdf': args.pdf,
            'src': args.src,
            'summary': {'figures': len(rows), 'fail': n_fail, 'warn': n_warn,
                        'pass': n_pass, 'skipped': len(skipped), 'large': len(big)},
            'figures': rows,
            'results': [{'level': a, 'cat': b, 'msg': c} for a, b, c in RESULTS],
        }, ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》图形版式自动守卫')
    print('=' * 78)
    print('[发现] %s/*.tex 里 %d 张 \\begin{tikzpicture}：跳过 %d 张，检查 %d 张'
          % (os.path.basename(os.path.normpath(args.src)), len(figs), len(skipped), len(rows)))
    for f in skipped:
        print('        跳过 %s:%d —— %s' % (f['file'], f['line'], f['skip']))
    unlocated = [r for r in rows if r.get('page') is None]
    print('[定位] 命中 %d ／ 未定位 %d（未定位的图在表里标 WARN，不静默通过）'
          % (len(rows) - len(unlocated), len(unlocated)))
    print('[几何] PASS %d ／ WARN %d ／ FAIL %d' % (n_pass, n_warn, n_fail))

    print('\n--- 体检表（FAIL → WARN → PASS）' + '-' * 40)
    print_table(rows)

    print('\n--- 明显过大的图（宽 > %.0f%% 版心 或 高 > %.0f cm）---' % (W_BIG * 100, H_HUGE))
    if big:
        for r in big:
            bits = []
            if (r.get('w_frac') or 0) > W_BIG:
                bits.append('宽 %.1f%% 版心（%.2f cm）' % (r['w_frac'] * 100, r['w_cm']))
            if (r.get('h_cm') or 0) > H_HUGE:
                bits.append('高 %.1f cm' % r['h_cm'])
            print('  %-30s p%-4s 图 %-6s %s'
                  % (r.get('label') or '(无 label)', r.get('page'), r.get('number') or '?',
                     '；'.join(bits)))
    else:
        print('  无（全部 ≤ %.0f%% 版心、高 ≤ %.0f cm）' % (W_BIG * 100, H_HUGE))

    if run_legacy:
        print('\n--- 旧词表兼容检查（第 9 章 3 张图 + 静态项；--no-legacy 可关）---')
        for lv, cat, msg in RESULTS:
            if lv in ('FAIL', 'WARN'):
                print('  [%s] %s' % (cat, msg))
        print('  PASS %d 项（逐项明细用 --json 看）'
              % sum(1 for r in RESULTS if r[0] == 'PASS'))

    nf = sum(1 for r in RESULTS if r[0] == 'FAIL')
    nw = sum(1 for r in RESULTS if r[0] == 'WARN')
    print('\n' + '=' * 78)
    print('汇总：图 %d 张（FAIL %d ／ WARN %d ／ PASS %d），全仓库检查项 FAIL %d ／ WARN %d'
          % (len(rows), n_fail, n_warn, n_pass, nf, nw))
    print('=' * 78)
    return 1 if nf else 0


if __name__ == '__main__':
    sys.exit(main())
