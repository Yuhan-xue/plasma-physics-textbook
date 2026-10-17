# -*- coding: utf-8 -*-
"""项目一览图：docs/overview.png（给人看的一页）。

数据来源：README.md 版本表、check/ 六套检查输出、tools/appendix_audit.json、check/ 守卫项数、pdf/ 实测页数。
用法：python -X utf8 tools/overview_fig.py
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

for name in ('Microsoft YaHei', 'SimHei', 'DengXian'):
    if any(f.name == name for f in font_manager.fontManager.ttflist):
        plt.rcParams['font.sans-serif'] = [name]
        break
plt.rcParams['axes.unicode_minus'] = False

INK = '#1d2433'
MUTED = '#5b6675'
ACCENT = '#2f6fb5'
LIGHT = '#6f9fd8'
PALE = '#a9c4e4'
FAINT = '#cfdcf0'
WARN = '#c1622a'
BG = '#f7f5f0'

fig = plt.figure(figsize=(13.6, 8.6), dpi=140)
fig.patch.set_facecolor(BG)
fig.text(0.035, 0.955, '等离子体物理自学教材 · 项目一览', fontsize=24, color=INK, weight='bold')
fig.text(0.035, 0.915, 'v1.27.5 ｜ 635 页 ｜ 12 章正文 + 49 附录 ｜ 95 张实拍图 ｜ 六套检查全绿',
         fontsize=13, color=MUTED)

# ── 左上：篇幅构成
ax1 = fig.add_axes([0.05, 0.635, 0.42, 0.165])
ax1.set_facecolor(BG)
parts = [('49 个附录', 427, ACCENT), ('12 章正文', 177, LIGHT),
         ('符号与前言', 21, PALE), ('结语与索引', 10, FAINT)]
left = 0
for label, val, color in parts:
    ax1.barh([0], [val], left=left, color=color, height=0.5)
    if val > 100:
        ax1.text(left + val / 2, 0, '%s · %d 页' % (label, val), ha='center', va='center',
                 fontsize=11.5, color='white' if color == ACCENT else INK, weight='bold')
    left += val
ax1.set_xlim(0, 635)
ax1.set_ylim(-0.55, 0.62)
ax1.axis('off')
ax1.set_title('635 页怎么分：三分之二在附录，正文 177 页', fontsize=13, color=INK,
              loc='left', pad=8)
ax1.text(0, -0.50, '49 个附录共 427 页；前 10 大附录只占其中 27%，最大一个 22 页',
         fontsize=10.3, color=MUTED)

# ── 右上：版本演进
ax2 = fig.add_axes([0.545, 0.615, 0.40, 0.21])
ax2.set_facecolor('white')
vers = ['v1.26.31', 'v1.27.0', 'v1.27.1', 'v1.27.2', 'v1.27.3', 'v1.27.5']
pages = [676, 674, 664, 668, 657, 635]
figs = [70, 87, 87, 95, 95, 95]
xs = list(range(len(vers)))
ax2.plot(xs, pages, marker='o', color=ACCENT, linewidth=2.4, label='页数')
offsets = [(0, 9), (0, 9), (0, 9), (0, 9), (-2, -18), (6, 8)]
for x, p, off in zip(xs, pages, offsets):
    ax2.annotate(str(p), (x, p), textcoords='offset points', xytext=off,
                 ha='center', fontsize=10.5, color=ACCENT, weight='bold')
ax2.set_xticks(xs)
ax2.set_xticklabels(vers, fontsize=10.2, color=INK)
ax2.set_ylim(628, 692)
ax2.set_ylabel('页数', fontsize=10.2, color=ACCENT)
ax2.tick_params(axis='y', colors=ACCENT, labelsize=9.8)
ax2.tick_params(axis='x', length=0)
ax2b = ax2.twinx()
ax2b.plot(xs, figs, marker='s', color=WARN, linewidth=2.2, linestyle='--', label='实拍图')
for x, f in zip(xs, figs):
    if False:
        ax2b.annotate(str(f), (x, f), textcoords='offset points', xytext=(0, 7),
                      ha='center', fontsize=10.5, color=WARN, weight='bold')
ax2b.set_ylim(52, 122)
ax2b.set_ylabel('实拍图（张）', fontsize=10.2, color=WARN)
ax2b.tick_params(axis='y', colors=WARN, labelsize=9.8)
ax2.set_title('一边减页、一边加图（蓝=页数，橙=实拍图 70 → 95）', fontsize=13, color=INK,
              loc='left', pad=8)

# ── 左中：六套检查
ax3 = fig.add_axes([0.05, 0.315, 0.28, 0.205])
ax3.set_facecolor('white')
ax3.axis('off')
ax3.set_title('六套检查：失败 0 项', fontsize=13, color=INK, loc='left', pad=22)
checks = [('数字独立复算', 7195), ('结构完整性', 49), ('前后一致性', 14),
          ('正文可读性', 59), ('插图与来源', 16), ('书签与页码', 11)]
for i, (name, n) in enumerate(checks):
    y = 0.78 - i * 0.145
    ax3.text(0.03, y, name, fontsize=11, color=INK, va='center')
    ax3.text(0.95, y, '%d 项' % n, fontsize=11, color=ACCENT, va='center', ha='right',
             weight='bold')
    if i < len(checks) - 1:
        ax3.plot([0.03, 0.95], [y - 0.062, y - 0.062], color='#e6e2da', linewidth=0.9)

# ── 中中：图形来源
ax4 = fig.add_axes([0.365, 0.315, 0.27, 0.205])
ax4.set_facecolor('white')
ax4.axis('off')
ax4.set_title('95 张图的出处', fontsize=13, color=INK, loc='left', pad=22)
ax4.set_xlim(0, 1)
ax4.set_ylim(0, 1)
srcs = [('政府与公共机构', 46, ACCENT),
        ('大学与实验室', 27, LIGHT),
        ('企业、个人、开放社区', 22, PALE)]
for i, (name, n, color) in enumerate(srcs):
    y = 0.72 - i * 0.30
    ax4.text(0.03, y + 0.135, name, fontsize=10.8, color=INK, va='center')
    ax4.add_patch(FancyBboxPatch((0.03, y - 0.075), 0.60 * n / 46.0, 0.13,
                                 boxstyle='round,pad=0.004', facecolor=color,
                                 edgecolor='none'))
    ax4.text(0.96, y, '%d 张' % n, fontsize=10.8, color=MUTED, ha='right', va='center')
ax4.text(0.03, 0.04, 'NASA / ITER / LLNL / INL / EUROfusion / ASDEX / KSTAR / TCV / NIF …',
         fontsize=8.8, color=MUTED, va='center')

# ── 右中：交付物
ax5 = fig.add_axes([0.675, 0.315, 0.28, 0.205])
ax5.set_facecolor('white')
ax5.axis('off')
ax5.set_title('读者拿得到的东西', fontsize=13, color=INK, loc='left', pad=22)
ax5.set_xlim(0, 1)
ax5.set_ylim(0, 1)
items = [('pdf/…v2.pdf', '635 页完整教材'),
         ('…vol1_chapters.pdf', '208 页：前言 + 12 章'),
         ('…vol2_appendices.pdf', '427 页：49 个附录'),
              ('…handbook.pdf', '26 页：配套手册 J/K/L'),
         ('figures/CREDITS.md', '95 张图的作者与许可')]
for i, (a, b) in enumerate(items):
    y = 0.74 - i * 0.235
    ax5.text(0.03, y, a, fontsize=10.6, color=ACCENT, va='center')
    ax5.text(0.03, y - 0.105, b, fontsize=9.6, color=MUTED, va='center')

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
         '数据来源：README.md 版本表、check/ 六套检查输出、tools/appendix_audit.json、'
         'pdf/ 实测页数（2026-10-18，v1.27.5）。',
         fontsize=9.4, color='#8b8577')

out = os.path.join(ROOT, 'docs', 'overview.png')
fig.savefig(out, facecolor=BG)
print(out)
