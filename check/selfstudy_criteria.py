# -*- coding: utf-8 -*-
r"""自学标准判据（供 check/full_check.py 第 R 类与 docs/self-study-audit.md 共用）。

三分对象各用各的标准，判据按本书既有约定写成（不要求把三种写法统一）：
  · 章节：正文 + 本章习题精选合并看。答案块接受「解答／答案／解：／分析／结论」五种写法
    ——本书的案例型例题用「分析＋结论」，计算型例题用「解答」，两种都是完整答案；
  · 工具型附录（A–L）：看结构、可检索与互链；
  · 专题附录（M–AZ）：看练习、答案、出处、适用边界、数据表与语义盒。
"""
from __future__ import print_function

import io
import os
import re

TOOL_APX = set('abcdefghijkl')

ANSWER = re.compile(u'\\\\textbf\\{(解答|答案|解：|解|分析|结论|数值估算|讨论)')
EXERCISE = re.compile(u'\\\\begin\\{example\\}')
ANYSEC = re.compile(u'\\\\section\\{([^}]*)\\}')

CH_CRITERIA = [
    (u'章首引导段', lambda b, s: bool(re.search(u'\\\\epigraph|chapterintro|quotation', b))),
    (u'核心目标盒', lambda b, s: u'核心目标' in b),
    (u'前置依赖', lambda b, s: u'前置依赖' in b),
    (u'例题 ≥3（正文+精选）', lambda b, s: len(EXERCISE.findall(b + s)) >= 3),
    (u'答案块成对', lambda b, s: len(ANSWER.findall(b + s)) >= len(EXERCISE.findall(b + s))),
    (u'先猜再看', lambda b, s: (b + s).count(u'\\predict{') >= 1),
    (u'解题策略盒', lambda b, s: b.count(u'\\begin{strategy}') >= 1),
    (u'章末概念图', lambda b, s: len(re.findall(u'fig:cm-ch', b)) >= 1),
    (u'量纲自检', lambda b, s: b.count(u'\\item \\textbf{量纲自检}') >= 1),
    (u'要点回顾', lambda b, s: b.count(u'chapterreview') >= 1),
    (u'随堂自测+参考要点', lambda b, s: b.count(u'\\begin{miniquiz}') >= 1
     and b.count(u'\\quizans{') >= 1),
    (u'配图 ≥1', lambda b, s: b.count(u'\\includegraphics') >= 1
     or b.count(u'\\begin{tikzpicture}') >= 1),
    (u'习题精选有答案', lambda b, s: len(ANSWER.findall(s)) >= 3),
]

TOOL_CRITERIA = [
    (u'小节结构 ≥5', lambda c: len(ANYSEC.findall(c)) >= 5),
    (u'表格/清单 ≥1', lambda c: u'\\label{tab:' in c or u'\\begin{itemize}' in c
     or u'\\begin{longtable}' in c),
    (u'标签可检索', lambda c: len(re.findall(u'\\\\label\\{', c)) >= 1),
    (u'与正文互链', lambda c: u'\\ref{' in c or bool(re.search(u'第~?\\s*\\d+\\s*章', c))
     or bool(re.search(u'附录~?[A-Z]', c))),
]

TOPIC_CRITERIA = [
    (u'估算练习小节', lambda c: any(u'估算练习' in m for m in ANYSEC.findall(c))),
    (u'答案 ≥6', lambda c: len(ANSWER.findall(c)) >= 6),
    (u'出处声明', lambda c: u'出处：' in c),
    (u'适用边界（不复制）', lambda c: u'不复制' in c),
    (u'数据表 ≥1', lambda c: u'\\label{tab:' in c),
    (u'语义盒 ≥1', lambda c: c.count(u'\\begin{tip}') + c.count(u'\\begin{insight}')
     + c.count(u'\\begin{warning}') >= 1),
]


def read(path):
    return io.open(path, encoding='utf-8', newline='').read()


def chapter_spans(text):
    idx = [m.start() for m in re.finditer(u'\\\\chapter\\{', text)]
    idx.append(len(text))
    return [(text[idx[i]:idx[i + 1]], i + 1) for i in range(len(idx) - 1)]


def appendix_files(src):
    apps = sorted(f for f in os.listdir(src)
                  if re.match(r'appendix_([a-z]{1,2})_.*\.tex$', f))
    tools = [f for f in apps
             if re.match(r'appendix_([a-z]{1,2})_', f).group(1) in TOOL_APX]
    topics = [f for f in apps if f not in tools]
    return tools, topics


def audit(src, main_name='plasma_physics_textbook_v2.tex'):
    """返回 (章节结果, 工具型结果, 专题型结果, 缺口列表)。

    结果形如 [(对象名, [(判据名, 是否通过)…])…]；缺口形如 [(对象名, 判据名)…]。
    """
    main_tex = read(os.path.join(src, main_name))
    ch_rows, gaps = [], []
    for body, num in chapter_spans(main_tex):
        sp = os.path.join(src, 'ch%d_supplement.tex' % num)
        sup = read(sp) if os.path.exists(sp) else u''
        cells = []
        for name, fn in CH_CRITERIA:
            ok = bool(fn(body, sup))
            cells.append((name, ok))
            if not ok:
                gaps.append((u'第 %d 章' % num, name))
        ch_rows.append((u'第 %d 章' % num, cells))

    tools, topics = appendix_files(src)
    tool_rows, topic_rows = [], []
    for f in tools:
        c = read(os.path.join(src, f))
        cells = []
        for name, fn in TOOL_CRITERIA:
            ok = bool(fn(c))
            cells.append((name, ok))
            if not ok:
                gaps.append((f, name))
        tool_rows.append((f, cells))
    for f in topics:
        c = read(os.path.join(src, f))
        cells = []
        for name, fn in TOPIC_CRITERIA:
            ok = bool(fn(c))
            cells.append((name, ok))
            if not ok:
                gaps.append((f, name))
        topic_rows.append((f, cells))
    return ch_rows, tool_rows, topic_rows, gaps


def markdown(src, main_name='plasma_physics_textbook_v2.tex'):
    """生成审计报告正文（Markdown）。"""
    ch_rows, tool_rows, topic_rows, gaps = audit(src, main_name)
    L = [u'# 自学标准审计（自动生成）', u'',
         u'由 `tools/selfstudy_audit.py --write` 生成，判据在 `check/selfstudy_criteria.py`'
         u'（与 `check/full_check.py` 第 R 类守卫共用同一份标准）。', u'',
         u'三分对象各用各的标准：章节看「正文 + 本章习题精选」；工具型附录（A–L）看结构、'
         u'可检索与互链；专题附录（M–AZ）看练习、答案、出处、适用边界、数据表与语义盒。', u'',
         u'答案块的写法本书有三种——计算型用「解答」、案例型用「分析＋结论」——判据同时接受，'
         u'不强行统一。', u'']

    def table(rows, criteria):
        L.append(u'| 对象 | ' + u' | '.join(c[0] for c in criteria) + u' |')
        L.append(u'| --- |' + u' --- |' * len(criteria))
        for name, cells in rows:
            L.append(u'| %s | %s |' % (name, u' | '.join(u'✓' if ok else u'✗'
                                                         for _n, ok in cells)))
        L.append(u'')

    L.append(u'## 一、第 1–12 章')
    L.append(u'')
    table(ch_rows, CH_CRITERIA)
    L.append(u'## 二、工具型附录（A–L）')
    L.append(u'')
    table(tool_rows, TOOL_CRITERIA)
    L.append(u'## 三、专题附录（M–AZ）')
    L.append(u'')
    table(topic_rows, TOPIC_CRITERIA)
    L.append(u'## 四、结论')
    L.append(u'')
    if not gaps:
        L.append(u'**无缺口**：%d 章、%d 个工具型附录、%d 个专题附录全部满足各自标准。'
                 % (len(ch_rows), len(tool_rows), len(topic_rows)))
    else:
        for who, name in gaps:
            L.append(u'- %s 缺：%s' % (who, name))
    L.append(u'')
    return u'\n'.join(L)
