# -*- coding: utf-8 -*-
"""生成 docs/appendix-audit.md：附录体检表（人类版）。

数据：tools/appendix_audit.json（页数/字节/图/表/盒）+ README 版本表 + PDF 书签。
用法：python -X utf8 tools/make_appendix_audit_doc.py [--write]
"""
import io
import json
import os
import re
import sys

import pymupdf

sys.stdout.reconfigure(encoding='utf-8')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = json.load(io.open(os.path.join(ROOT, 'tools', 'appendix_audit.json'), encoding='utf-8'))
doc = pymupdf.open(os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf'))
toc = [t for t in doc.get_toc() if t[0] == 1]
total = doc.page_count

# 章页数
chapters = []
for i, (lvl, title, page) in enumerate(toc):
    if re.match(r'^(数学与统计|什么是等离子体|单粒子运动|流体描述|磁流体力学|等离子体中的波|动理学理论|'
                r'碰撞与输运|等离子体不稳定性|等离子体约束|低温工业|等离子体基础诊断)', title):
        nxt = None
        for lvl2, t2, p2 in toc[i + 1:]:
            if re.match(r'^(数学与统计|什么是等离子体|单粒子运动|流体描述|磁流体力学|等离子体中的波|动理学理论|'
                        r'碰撞与输运|等离子体不稳定性|等离子体约束|低温工业|等离子体基础诊断|结语)', t2):
                nxt = p2
                break
        chapters.append((title, page, (nxt or total + 1) - page))

LAYOUT = {'I': '研究', 'J': '研究', 'K': '研究', 'L': '研究'}
GROUPS = [('工具型附录 A–L', [r['letter'] for r in rows if len(r['letter']) == 1
                            and r['letter'] <= 'L']),
          ('专题附录 M–Z', [r['letter'] for r in rows if len(r['letter']) == 1
                          and r['letter'] > 'L']),
          ('专题附录 AA–AZ', [r['letter'] for r in rows if len(r['letter']) == 2])]
by_letter = {r['letter']: r for r in rows}

appendix_pages = sum(r['pages'] for r in rows)
top10 = sorted(rows, key=lambda r: -r['pages'])[:10]

lines = []
lines.append('# 附录体检：%d 个附录的体量与可动空间' % len(rows))
lines.append('')
lines.append('> 这份表回答一个问题：**书已经很薄了，还想再薄，能动哪里？**'
             '结论先说：能动的地方不多，而且都在「学习辅助」类附录里，不在物理正文里。')
lines.append('')
lines.append('![项目一览](overview.png)')
lines.append('')
lines.append('## 一、体检结论')
lines.append('')
front_pages = toc[1][2] - 1                     # 封面到前言结束
ch_pages = sum(c[2] for c in chapters)
tail_pages = total - front_pages - ch_pages - appendix_pages
lines.append('- 全 **%d 页**：封面与目录 %d 页，前言 %d 页，正文 12 章 %d 页，'
             '附录 %d 个 **%d 页（%.0f%%）**，结语与索引 %d 页。'
             % (total, toc[0][2] - 1, front_pages - (toc[0][2] - 1), ch_pages,
                len(rows), appendix_pages, 100.0 * appendix_pages / total, tail_pages))
lines.append('- 专题附录（M–AZ）体量非常均匀：**多数 7–12 页**，每篇 7–8 张表、1–3 张图、6–7 个语义盒——'
             '这是模板化的结果，也意味着单篇里没有明显的水分。')
lines.append('- 体量最大的五个：%s。'
             % '、'.join('附录 %s（%d 页）' % (r['letter'], r['pages']) for r in top10[:5]))
lines.append('- 前 10 大附录只占附录总页数的 **27%**，即"砍掉最大的几篇"最多省三四十页，'
             '而且会连带删掉守卫绑定的数字。')
lines.append('')
lines.append('## 二、正文 12 章的页数')
lines.append('')
lines.append('| 章 | 起始页 | 页数 |')
lines.append('|---|---:|---:|')
for title, page, n in chapters:
    lines.append('| %s | %d | %d |' % (title, page, n))
lines.append('')
lines.append('## 三、%d 个附录逐个体检（附录 J／K／L 已移入配套手册）' % len(rows))
lines.append('')
for gname, letters in GROUPS:
    lines.append('### %s' % gname)
    lines.append('')
    lines.append('| 附录 | 标题 | 页数 | 图 | 表 | 语义盒 | 建议 |')
    lines.append('|---|---|---:|---:|---:|---:|---|')
    for letter in letters:
        r = by_letter.get(letter)
        if not r:
            continue
        if letter in 'ABCDEFGH':
            advice = '参考表·建议保持'
        elif letter == 'I':
            advice = '**可以瘦身**：非物理内容占多数'
        elif letter in ('J', 'K', 'L'):
            advice = '可考虑独立成册'
        elif r['pages'] >= 12:
            advice = '偏长·可微调'
        else:
            advice = '保持'
        lines.append('| %s | %s | %d | %d | %d | %d | %s |'
                     % (letter, r['title'], r['pages'], r['figures'], r['tables'],
                        r['boxes'], advice))
    lines.append('')
lines.append('## 四、已经做完的一刀，以及还能动的四处')
lines.append('')
lines.append('**已经做完（v1.27.5）**：附录 J（经典教材逐章对照）、K（跨教材经典例题）、'
             'L（分章自测题库）移出主册，编成 26 页的配套手册'
             '（`pdf/plasma_physics_handbook.pdf`），主册 657 → 635 页；'
             '主册附录字母保留 J／K／L 空位，A–I 之后直接接 M。')
lines.append('')
lines.append('还能动的四处（都要同步改断言，见下一节）：')
lines.append('')
lines.append('1. **附录 I（学习策略手册，%d 页）**——主册最大的单篇，讲的是怎么学、不是物理。'
             '可以把 6 篇核心方法留成 6–8 页，其余（论文三遍法、复习排程、同伴教学等）'
             '移到仓库文档里，需要的人按链接查。'
             % by_letter['I']['pages'])
lines.append('2. **主题相近的专题附录合并**——例如 AC（数值方法）与 AE（动理学）的数值部分、'
             'AQ（干涉诊断）与 AN（偏振诊断）。合并每对大约省 6–8 页，'
             '但会打乱「一附录一主题」的查阅习惯。')
lines.append('3. **细排范围再扩一档**——v1.27.3 已把「来源、分级与版权说明」小节与出处段改成细排；'
             '各附录开头的「与相邻附录的分工」段还没动，改小一号字大约省 2–4 页。')
lines.append('4. **图片再收一档**——竖幅照片已按高度定尺、横幅已收到 0.72 页宽；'
             '再收到 0.68 大约省 3–5 页，但会牺牲一点视觉冲击。')

lines.append('## 五、动页数的代价（必须知道的一条）')
lines.append('')
lines.append('- 书里的每个数字都被 `check/physics_regression.py`（7195 项断言）'
             '**独立重算并要求留在它自己那一节里**；表格行、语义盒标签、'
             '`不复制` 声明、`出处：` 段落同样是绑定对象。')
lines.append('- 所以删内容不是"删一段话"，而是"删一段话 + 改掉绑定它的断言"。'
             '2026-10-17 试过一次整节删除（41 个来源小节），直接换来 234 条断言失败，已回滚。')
lines.append('- 稳妥做法：一次只动一个附录，改完立刻跑 `python -X utf8 check/physics_regression.py`'
             '（读末尾"汇总"行），再重编核对页数。')
lines.append('')
lines.append('## 六、这份表怎么来的')
lines.append('')
lines.append('- 页数：从 PDF 书签取每个附录的起始页，逐篇相减（附录内联在主文件里的 A–H 也一并计入）。')
lines.append('- 体量与结构：统计 `src/appendix_*.tex` 的字节数、`\\section` 数、'
             '`\\includegraphics` 数、`table/longtable` 数与语义盒数。')
lines.append('- 生成脚本：`tools/appendix_audit.py`（`--json` 可导出明细）。')

out = os.path.join(ROOT, 'docs', 'appendix-audit.md')
text = '\n'.join(lines) + '\n'
if '--write' in sys.argv:
    io.open(out, 'w', encoding='utf-8', newline='\n').write(text)
    print('写出', out, len(text), '字符')
else:
    print(text[:1500])
