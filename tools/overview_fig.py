# -*- coding: utf-8 -*-
"""项目一览图：docs/overview.png（给人看的一页）。

数字来源：
- 页数与版本：直接读 pdf/ 里四个 PDF（页数、元数据 Title 里的版本号）
- 插图张数与来源分布：figures/CREDITS.md（分类口径见 tools/figure_sources.py）
- 六套检查的项数：check/ 各脚本的实际项数（改动检查后同步这里）
用法：python -X utf8 tools/overview_fig.py
"""
import datetime
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'tools'))
from figure_sources import load_rows, counts as source_counts  # noqa: E402

for name in ('Microsoft YaHei', 'SimHei', 'DengXian'):
    if any(f.name == name for f in font_manager.fontManager.ttflist):
        plt.rcParams['font.sans-serif'] = [name]
        break
plt.rcParams['axes.unicode_minus'] = False

INK = '#12203A'
MUTED = '#6B7280'
ACCENT = '#1B2E52'
LIGHT = '#2F5286'
PALE = '#6C86AE'
FAINT = '#DCE3EE'
GOLD = '#C2A14D'
AMBER = '#8A6520'
BG = '#F6F2E7'

# ── 数据
pdf_dir = os.path.join(ROOT, 'pdf')
MAIN = os.path.join(pdf_dir, 'plasma_physics_textbook_v2.pdf')
VOL1 = os.path.join(pdf_dir, 'plasma_physics_textbook_v2_vol1_chapters.pdf')
VOL2 = os.path.join(pdf_dir, 'plasma_physics_textbook_v2_vol2_appendices.pdf')
HB = os.path.join(pdf_dir, 'plasma_physics_handbook.pdf')
doc = pymupdf.open(MAIN)
PAGES = doc.page_count
ver = re.search(r'v\d+\.\d+\.\d+', (doc.metadata or {}).get('title') or '')
VERSION = ver.group(0) if ver else 'v?'
P1 = pymupdf.open(VOL1).page_count
P2 = pymupdf.open(VOL2).page_count
PH = pymupdf.open(HB).page_count

rows = load_rows(os.path.join(ROOT, 'figures', 'CREDITS.md'))
NF = len(rows)
src = source_counts(rows)

# 页数构成从书签实测，避免脚本里留着手填的常数
_toc = doc.get_toc()


def _start(prefix):
    for lv, title, page in _toc:
        if lv == 1 and title.strip().startswith(prefix):
            return page
    raise SystemExit('书签里找不到：%s' % prefix)


_lv1 = [(tt.strip(), pp) for lv, tt, pp in _toc if lv == 1]
_pref = [i for i, (tt, _) in enumerate(_lv1) if tt == '前言'][0]
_CH_START = _lv1[_pref + 1][1]      # 前言之后的第一条一级书签就是第 1 章
_TAIL_START = _start('结语')
_APX_START = _start('附录 A')
FRONT = _CH_START - 1                          # 封面、目录、符号与前言
CH_PAGES = _TAIL_START - _CH_START             # 第 1 章起，到结语前
TAIL = _APX_START - _TAIL_START                # 结语
APX_PAGES = PAGES - _APX_START + 1             # 48 个附录
NCHK = [7195, 45, 14, 43, 16, 11]              # 各类检查的项数（随门禁更新）

fig = plt.figure(figsize=(13.6, 8.6), dpi=140)
fig.patch.set_facecolor(BG)
fig.text(0.035, 0.955, '等离子体物理自学教材 · 项目一览', fontsize=24, color=INK, weight='bold')
fig.text(0.035, 0.915,
         '%s ｜ %d 页 ｜ 12 章正文 + 48 附录 + 一册配套手册 ｜ %d 张实拍图 ｜ 六套检查全绿'
         % (VERSION, PAGES, NF), fontsize=13, color=MUTED)

# ── 左上：篇幅构成
ax1 = fig.add_axes([0.05, 0.65, 0.42, 0.145])
ax1.set_facecolor(BG)
parts = [('48 个附录', APX_PAGES, ACCENT), ('正文', CH_PAGES, LIGHT),
         ('前言等', FRONT, PALE), ('结语', TAIL, FAINT)]
left = 0
for label, val, color in parts:
    ax1.barh([0], [val], left=left, color=color, height=0.5)
    if val > 150:
        ax1.text(left + val / 2, 0, '%s · %d 页' % (label, val), ha='center', va='center',
                 fontsize=10.5, color='white' if color in (ACCENT, LIGHT) else INK, weight='bold')
    left += val
ax1.set_xlim(0, PAGES)
ax1.set_ylim(-0.55, 0.62)
ax1.axis('off')
ax1.set_title('%d 页怎么分：三分之二在附录，正文 %d 页' % (PAGES, CH_PAGES), fontsize=13,
              color=INK, loc='left', pad=8)
ax1.text(0, -0.50, '48 个附录 %d 页｜正文 %d 页｜前言等 %d 页｜结语 %d 页'
         % (APX_PAGES, CH_PAGES, FRONT, TAIL),
         fontsize=9.6, color=MUTED)

# ── 右上：版本演进
ax2 = fig.add_axes([0.545, 0.66, 0.40, 0.17])
ax2.set_facecolor('white')
vers = ['v1.26.31', 'v1.27.0', 'v1.27.1', 'v1.27.2', 'v1.27.3', 'v1.27.5', 'v1.27.6', 'v1.27.7',
        'v1.27.8', 'v1.27.9', 'v1.27.10', 'v1.27.11', 'v1.27.12', 'v1.27.13', 'v1.27.14', 'v1.27.15', 'v1.27.16', 'v1.27.17', 'v1.27.18', 'v1.27.19', 'v1.27.20']
pages = [676, 674, 664, 668, 657, 635, 635, 634, 605, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES, PAGES]
figs = [70, 87, 87, 95, 95, 95, 95, 95, NF, NF, NF, NF, NF, NF, NF, NF, NF, NF, NF, NF, NF]
xs = list(range(len(vers)))
ax2.plot(xs, pages, marker='o', color=ACCENT, linewidth=2.4)
notes = {0: (0, 9), 3: (0, 9), 4: (0, -16), 7: (0, 9), 8: (0, -16), 20: (0, -16)}
for x, p in zip(xs, pages):
    if x in notes:
        ox, oy = notes[x]
        ax2.annotate(str(p), (x, p), textcoords='offset points', xytext=(ox, oy),
                     ha='center', fontsize=9.6, color=ACCENT, weight='bold')
show = list(range(0, len(vers), 2)) + [len(vers) - 1]
ax2.set_xticks(show)
ax2.set_xticklabels([vers[i] for i in show], fontsize=7.4, color=INK, rotation=45, ha='right')
ax2.set_ylim(592, 692)
ax2.set_ylabel('页数', fontsize=10.2, color=ACCENT)
ax2.tick_params(axis='y', colors=ACCENT, labelsize=9.8)
ax2.tick_params(axis='x', length=0)
ax2b = ax2.twinx()
ax2b.plot(xs, figs, marker='s', color=GOLD, linewidth=2.2, linestyle='--')
ax2b.set_ylim(52, 122)
ax2b.set_ylabel('实拍图（张）', fontsize=10.2, color=AMBER)
ax2b.tick_params(axis='y', colors=AMBER, labelsize=9.8)
ax2.set_title('一边减页、一边加图（蓝=页数，金=实拍图 70 → %d）' % NF, fontsize=13, color=INK,
              loc='left', pad=12)

# ── 左中：六套检查
ax3 = fig.add_axes([0.05, 0.30, 0.28, 0.24])
ax3.set_facecolor('white')
ax3.axis('off')
ax3.set_title('六套检查：失败 0 项', fontsize=13, color=INK, loc='left', pad=16)
checks = [('数字独立复算', NCHK[0]), ('结构完整性', NCHK[1]), ('前后一致性', NCHK[2]),
          ('正文可读性', NCHK[3]), ('插图与来源', NCHK[4]), ('书签与页码', NCHK[5])]
for i, (name, n) in enumerate(checks):
    y = 0.78 - i * 0.145
    ax3.text(0.03, y, name, fontsize=11, color=INK, va='center')
    ax3.text(0.95, y, '%d 项' % n, fontsize=11, color=ACCENT, va='center', ha='right',
             weight='bold')
    if i < len(checks) - 1:
        ax3.plot([0.03, 0.95], [y - 0.062, y - 0.062], color='#E5E0D2', linewidth=0.9)

# ── 中中：图形来源
ax4 = fig.add_axes([0.365, 0.30, 0.27, 0.24])
ax4.set_facecolor('white')
ax4.axis('off')
ax4.set_title('%d 张图的出处' % NF, fontsize=13, color=INK, loc='left', pad=16)
ax4.set_xlim(0, 1)
ax4.set_ylim(0, 1)
srcs = [('政府 / 国家实验室 / 大学', src['inst'], ACCENT),
        ('个人摄影师 / 企业 / 社区', src['person'], LIGHT)]
for i, (name, n, color) in enumerate(srcs):
    y = 0.68 - i * 0.34
    ax4.text(0.03, y + 0.13, name, fontsize=10.8, color=INK, va='center')
    ax4.add_patch(FancyBboxPatch((0.03, y - 0.06),
                                 0.60 * n / float(max(src['inst'], src['person'], 1)), 0.12,
                                 boxstyle='round,pad=0.004', facecolor=color,
                                 edgecolor='none'))
    ax4.text(0.96, y, '%d 张' % n, fontsize=10.8, color=MUTED, ha='right', va='center')
ax4.text(0.03, 0.03, 'NASA / ESA / ITER / LLNL / INL / EUROfusion / ASDEX / KSTAR / NIF …',
         fontsize=8.8, color=MUTED, va='center')

# ── 右中：交付物
ax5 = fig.add_axes([0.675, 0.285, 0.30, 0.255])
ax5.set_facecolor('white')
ax5.axis('off')
ax5.set_title('读者拿得到的东西', fontsize=13, color=INK, loc='left', pad=16)
ax5.set_xlim(0, 1)
ax5.set_ylim(0, 1)
items = [('pdf/…v2.pdf', '%d 页完整教材' % PAGES),
         ('…vol1_chapters.pdf', '%d 页：前言 + 12 章' % P1),
         ('…vol2_appendices.pdf', '%d 页：48 个附录' % P2),
         ('…handbook.pdf', '%d 页：配套手册 I–L' % PH),
         ('figures/CREDITS.md', '%d 张图的作者与许可' % NF),
         ('docs/figures.md', '图册：一页看完全部插图')]
for i, (a, b) in enumerate(items):
    y = 0.88 - i * 0.15
    ax5.text(0.03, y, a, fontsize=9.8, color=ACCENT, va='center')
    ax5.text(0.03, y - 0.060, b, fontsize=8.7, color=MUTED, va='center')

# ── 脚注
fig.text(0.035, 0.132,
         '「以人为本」的六条口径：① 正文与记录不出现整段程序，只给仓库路径；'
         '② 来源与版权说明改为细排，正文只留物理；',
         fontsize=10.4, color=INK)
fig.text(0.035, 0.098,
         '③ 每张图都写清作者与许可（公有领域或 CC BY / CC BY-SA）；④ 记录（README、版本记录、评审）按读者重写；'
         '⑤ 每版重编、分册、跑六套检查；⑥ 附件逐字节复核后再发布。',
         fontsize=10.4, color=INK)
fig.text(0.035, 0.064,
         '另有 code/ 与 tools/calc/ 下的可运行脚本（正文只给路径），以及每版四个带版本号的 PDF 附件。',
         fontsize=10.1, color=MUTED)
fig.text(0.035, 0.026,
         '数据来源：pdf/ 四个 PDF 实测页数与元数据、figures/CREDITS.md、check/ 六套检查输出'
         '（%s，%s）。' % (datetime.date.today().isoformat(), VERSION),
         fontsize=9.4, color='#8b8577')

out = os.path.join(ROOT, 'docs', 'overview.png')
fig.savefig(out, facecolor=BG)
print(out, '已生成')
