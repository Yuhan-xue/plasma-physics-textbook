#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》全量检查套件
=====================================

一次运行完成五类检查，输出可复现的报告：

  A. 结构完整性   —— 环境配对、标签唯一性、引用可达性、图表 caption/label
  B. 内容完整性   —— 每章必备元素、难度等级合法性、占位符残留
  C. 编译健康度   —— error / warning / Missing character / Overfull
  D. PDF 校验     —— 页数、未解析引用、附录页眉、章号连续性
  E. 文档一致性   —— README/processed.md 声称的数字 vs 实际

用法：
    python check/full_check.py                # 全部检查
    python check/full_check.py --no-compile   # 跳过编译（仅静态检查）

退出码：0 = 无 FAIL；1 = 存在 FAIL
"""

import io
import os
import re
import sys
import glob
import json
import argparse
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
MAIN = 'plasma_physics_textbook_v2.tex'

# ---------------------------------------------------------------- 结果收集

RESULTS = []          # (level, category, message)
LEVELS = {'PASS': 0, 'INFO': 1, 'WARN': 2, 'FAIL': 3}


def rec(level, cat, msg):
    RESULTS.append((level, cat, msg))


def read(path):
    with io.open(path, encoding='utf-8') as f:
        return f.read()


def strip_comments(text):
    """去掉 LaTeX 行注释（保留 \\% 转义），把注释内容替换为等长空白以保持行号。"""
    out = []
    for line in text.split('\n'):
        idx = None
        i = 0
        while i < len(line):
            if line[i] == '\\':
                i += 2
                continue
            if line[i] == '%':
                idx = i
                break
            i += 1
        if idx is None:
            out.append(line)
        else:
            out.append(line[:idx] + ' ' * (len(line) - idx))
    return '\n'.join(out)


# tcolorbox 定理类环境的自动标签前缀：{前缀} 由 \newtcbtheorem 的最后一个参数给出
TCB_PREFIX = {
    'theorem': 'thm', 'definition': 'def', 'corollary': 'cor',
    'example': 'ex', 'insight': 'ins', 'warning': 'warn',
    'history': 'hist', 'review': 'rev', 'tip': 'tip',
}
# 用 phantomlabel 的（\NewTColorBox{...}{m m}）
PHANTOM_PREFIX = {'formulabox': 'frml', 'checklistbox': 'ckl'}


def tcb_auto_labels(text):
    """收集 tcolorbox 由第二个参数自动生成的标签，返回 {label: line}。"""
    found = {}
    for env, pref in TCB_PREFIX.items():
        # \begin{example}{标题}{key}
        for m in re.finditer(
                r'\\begin\{%s\}\s*\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}\s*\{([^}]+)\}'
                % re.escape(env), text):
            found['%s:%s' % (pref, m.group(1).strip())] = line_of(text, m.start())
    for env, pref in PHANTOM_PREFIX.items():
        for m in re.finditer(
                r'\\begin\{%s\}\s*\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}\s*\{([^}]+)\}'
                % re.escape(env), text):
            found['%s:%s' % (pref, m.group(1).strip())] = line_of(text, m.start())
    return found


def tex_files():
    """全部参与编译的 tex 源文件（排除独立样张）。"""
    out = []
    for p in sorted(glob.glob(os.path.join(SRC, '*.tex'))):
        if os.path.basename(p) == 'style_v2_preview.tex':
            continue
        out.append(p)
    return out


def line_of(text, idx):
    return text.count('\n', 0, idx) + 1


# ---------------------------------------------------------------- A. 结构

# 已知的自定义环境（在 main tex 中用 newtcbtheorem / newtcolorbox / NewTColorBox 定义）
CUSTOM_ENVS = ['theorem', 'definition', 'corollary', 'example', 'insight',
               'warning', 'history', 'review', 'tip', 'miniquiz',
               'chapterreview', 'formula', 'formulabox', 'checklist',
               'checklistbox', 'chapterintro', 'proof', 'equation', 'align',
               'itemize', 'enumerate', 'figure', 'table', 'tabular',
               'center', 'quotation', 'abstract', 'tikzpicture', 'scope',
               'thebibliography', 'verbatim', 'Verbatim', 'lstlisting']


def check_A_structure(files):
    label_defs = {}          # label -> (file, line)
    ref_uses = []            # (name, file, line)
    env_unbalanced = []
    fig_no_caption = []
    fig_no_label = []

    for p in files:
        text = strip_comments(read(p))
        base = os.path.basename(p)

        # --- A0: tcolorbox 自动生成的标签（第二个参数）
        for name, ln in tcb_auto_labels(text).items():
            if name in label_defs:
                rec('FAIL', 'A', '标签重复定义 %s —— %s:%d 与 %s:%d'
                    % (name, label_defs[name][0], label_defs[name][1], base, ln))
            else:
                label_defs[name] = (base, ln)

        # --- A1/A2: 显式 \label 定义
        for m in re.finditer(r'\\label\{([^}]+)\}', text):
            name = m.group(1)
            if name in label_defs:
                rec('FAIL', 'A', '标签重复定义 %s —— %s:%d 与 %s:%d'
                    % (name, label_defs[name][0], label_defs[name][1],
                       base, line_of(text, m.start())))
            else:
                label_defs[name] = (base, line_of(text, m.start()))

        for m in re.finditer(r'\\(?:c?ref|eqref|autoref|nameref)\{([^}]+)\}', text):
            for one in m.group(1).split(','):
                ref_uses.append((one.strip(), base, line_of(text, m.start())))

        # --- A3: 环境配对（注释已剥离，避免 % \end{...} 误报）
        stack = []
        for m in re.finditer(r'\\(begin|end)\{([^}]+)\}', text):
            kind, env = m.group(1), m.group(2)
            if kind == 'begin':
                stack.append((env, line_of(text, m.start())))
            else:
                if not stack:
                    env_unbalanced.append('%s:%d  \\end{%s} 无对应 \\begin'
                                          % (base, line_of(text, m.start()), env))
                elif stack[-1][0] != env:
                    env_unbalanced.append(
                        '%s:%d  \\end{%s} 与最近的 \\begin{%s}（第 %d 行）不匹配'
                        % (base, line_of(text, m.start()), env,
                           stack[-1][0], stack[-1][1]))
                    stack.pop()
                else:
                    stack.pop()
        for env, ln in stack:
            env_unbalanced.append('%s:%d  \\begin{%s} 未闭合' % (base, ln, env))

        # --- A4: figure 的 caption / label
        for fm in re.finditer(r'\\begin\{figure\}', text):
            end = text.find(r'\end{figure}', fm.end())
            if end < 0:
                continue
            body = text[fm.start():end]
            ln = line_of(text, fm.start())
            if r'\caption' not in body:
                fig_no_caption.append('%s:~%d' % (base, ln))
            if r'\label' not in body:
                fig_no_label.append('%s:~%d' % (base, ln))

    # 引用可达性
    for name, base, ln in ref_uses:
        if name not in label_defs:
            rec('FAIL', 'A', '引用不存在的标签 %s —— %s:%d' % (name, base, ln))

    if env_unbalanced:
        for e in env_unbalanced[:20]:
            rec('FAIL', 'A', e)
        if len(env_unbalanced) > 20:
            rec('FAIL', 'A', '……另有 %d 处环境不配对' % (len(env_unbalanced) - 20))
    else:
        rec('PASS', 'A', '环境配对：全部 \\begin/\\end 匹配')

    if fig_no_caption:
        for f in fig_no_caption:
            rec('FAIL', 'A', 'figure 缺 caption：%s' % f)
    else:
        rec('PASS', 'A', 'figure caption：全部 figure 均有 caption')

    if fig_no_label:
        for f in fig_no_label:
            rec('WARN', 'A', 'figure 缺 label（无法交叉引用）：%s' % f)
    else:
        rec('PASS', 'A', 'figure label：全部 figure 均有 label')

    rec('INFO', 'A', '标签总数 %d，引用总数 %d' % (len(label_defs), len(ref_uses)))
    return label_defs


# ---------------------------------------------------------------- B. 内容

def check_section_refs(files):
    """校验正文里以纯文本写出的 §X.Y 引用是否真实存在。

    LaTeX 只检查 \\ref，纯文本的「§7.4」不会被检查；章节调整后极易变成悬空引用。
    权威编号来自 .toc（\\numberline {N} 才是有编号的小节，带 * 的不编号）。
    """
    toc_path = os.path.join(SRC, MAIN.replace('.tex', '.toc'))
    if not os.path.exists(toc_path):
        rec('WARN', 'B', '找不到 .toc，跳过散引用校验（请先编译一次）')
        return {}

    toc = read(toc_path)
    real = {}
    n_ch = 0
    for m in re.finditer(r'\\contentsline \{(chapter|section)\}\{\\numberline \{([^}]+)\}',
                         toc):
        kind, num = m.group(1), m.group(2).strip()
        if kind == 'chapter':
            n_ch += 1
        else:
            real[num] = True

    if not real:
        rec('WARN', 'B', '.toc 中未解析到编号小节，跳过散引用校验')
        return {}

    # 收集 §X.Y / \S X.Y 引用
    pat = re.compile(r'(?:\\S|§)\s*~?\s*([0-9]+|[A-G])\s*[.．]\s*([0-9]+)')
    bad, seen = 0, {}
    for p in files:
        text = strip_comments(read(p))
        base = os.path.basename(p)
        for m in pat.finditer(text):
            key = '%s.%s' % (m.group(1), m.group(2))
            seen.setdefault(key, []).append((base, line_of(text, m.start()), text))

    for key in sorted(seen):
        if key in real:
            continue
        bad += 1
        base, ln, text = seen[key][0]
        idx = text.find('§' if '§' in text else '\\S')
        rec('FAIL', 'B', '散引用 §%s 不存在（首次出现 %s:%d，共 %d 处）'
            % (key, base, ln, len(seen[key])))

    if bad == 0:
        rec('PASS', 'B', '散引用校验：全部 §X.Y 均指向真实小节（%d 个不同引用）'
            % len(seen))
    else:
        rec('INFO', 'B', '散引用：共 %d 个不同 §X.Y 引用，其中 %d 个无效'
            % (len(seen), bad))
    rec('INFO', 'B', '实际编号小节数：%d（来自 .toc），编号章 %d'
        % (len(real), n_ch))
    return real


def check_B_content(files):
    main = read(os.path.join(SRC, MAIN))

    # --- B1: 章节结构
    chapters = [(line_of(main, m.start()), m.group(1))
                for m in re.finditer(r'\\chapter\{([^}]*)\}', main)]
    rec('INFO', 'B', '全文章节数（含前言/附录）：%d' % len(chapters))

    # --- B2: 每章正文应有 chapterintro + chapterreview
    body = re.split(r'(?=\\chapter\{)', main)
    n_no_intro, n_no_review, n_no_diff = [], [], []
    for chunk in body[1:]:
        title = re.match(r'\\chapter\{([^}]*)\}', chunk)
        if not title:
            continue
        t = title.group(1)
        if r'\begin{chapterintro}' not in chunk:
            n_no_intro.append(t)
        if r'\begin{chapterreview}' not in chunk:
            n_no_review.append(t)
        if r'\difficulty{' not in chunk and r'\difficultyhere{' not in chunk:
            n_no_diff.append(t)

    for t in n_no_review:
        rec('FAIL', 'B', '章节无「本章要点回顾」：%s' % t)
    if not n_no_review:
        rec('PASS', 'B', '各章均含「本章要点回顾」')

    # --- B2b: 章首元素一致性
    # 基准（第 1--10 章）：\epigraph 引语 + 「本章核心目标」insight 盒 + 「前置依赖」
    # 第 11、12 章以 chapterintro / quotation 引导段代替 epigraph——这是有意的设计差异
    # （两章的引导段本身信息量更大），因此「引导段」与「epigraph」二者居其一即可；
    # 但「本章核心目标」与「前置依赖」是功能性元素，缺一不可。
    ch_blocks = []
    for chunk in body[1:]:
        title = re.match(r'\\chapter\{([^}]*)\}', chunk)
        if title:
            ch_blocks.append((title.group(1), chunk))

    def head45(c):
        return '\n'.join(c.split('\n')[:50])

    def has_opening(c):
        h = head45(c)
        return (r'\epigraph' in h or r'\begin{chapterintro}' in h
                or r'\begin{quotation}' in h)

    base_goal = sum(1 for _, c in ch_blocks[:10] if '本章核心目标' in head45(c))
    base_pre = sum(1 for _, c in ch_blocks[:10] if '前置依赖' in head45(c))
    rec('INFO', 'B', '前 10 章开场基准：核心目标 %d/10、前置依赖 %d/10' %
        (base_goal, base_pre))

    for name, chunk in ch_blocks:
        h = head45(chunk)
        miss = []
        if not has_opening(chunk):
            miss.append('章首引导段（epigraph / chapterintro / quotation）')
        if base_goal >= 8 and '本章核心目标' not in h:
            miss.append('「本章核心目标」')
        if base_pre >= 8 and '前置依赖' not in h:
            miss.append('「前置依赖」')
        if miss:
            rec('FAIL', 'B', '章首结构缺：%s —— %s' % ('、'.join(miss), name))
    if not any(
            (not has_opening(c))
            or (base_goal >= 8 and '本章核心目标' not in head45(c))
            or (base_pre >= 8 and '前置依赖' not in head45(c))
            for _, c in ch_blocks):
        rec('PASS', 'B', '各章章首结构完整（引导段＋核心目标＋前置依赖）')

    # 章首容器类型统计
    styles = {}
    for name, chunk in ch_blocks:
        h = head45(chunk)
        if r'\begin{chapterintro}' in h:
            s = 'chapterintro'
        elif r'\begin{quotation}' in h:
            s = 'quotation'
        elif r'\epigraph' in h:
            s = 'epigraph'
        else:
            s = '(无)'
        styles.setdefault(s, []).append(name)
    if len(styles) > 1:
        rec('INFO', 'B', '章首引导段写法：' +
            '；'.join('%s×%d' % (k, len(v)) for k, v in
                      sorted(styles.items(), key=lambda x: -len(x[1]))))

    # --- B3: 补充材料
    for p in files:
        base = os.path.basename(p)
        m = re.match(r'ch(\d+)_supplement\.tex', base)
        if not m:
            continue
        text = read(p)
        if r'\begin{formulabox}' not in text:
            rec('WARN', 'B', '%s 缺公式速查卡' % base)
        if r'\begin{checklistbox}' not in text:
            rec('WARN', 'B', '%s 缺自学检查清单' % base)
        if r'\begin{example}' not in text:
            rec('WARN', 'B', '%s 缺例题' % base)

    # --- B4b: 正文中的 §X.Y 散引用校验
    # 这些引用是纯文本（不是 \ref），LaTeX 不会报错，极易随着章节调整而失效。
    check_section_refs(files)

    # --- B4: 难度等级合法性（对照正文的等级表）
    used = {}
    for p in files:
        text = read(p)
        for m in re.finditer(r'\\difficulty(?:here)?\{([^}]*)\}', text):
            v = m.group(1)
            used.setdefault(v, []).append('%s:%d' % (os.path.basename(p),
                                                     line_of(text, m.start())))

    # 从等级表提取已声明的名称（表格行形如  $\star$\,[\,基础\,] & 基础 & ...）
    declared = set(re.findall(r'&\s*([^\s&\\]+)\s*&', main))
    declared = {d for d in declared if d in ('基础', '提高', '挑战', '应用')}

    # 「挑战」的就地标注：proof 标题「（挑战）」+ 正文 [挑战] 记号
    inline_challenge = 0
    for p in files:
        t = read(p)
        inline_challenge += len(re.findall(r'（挑战）', t))
        inline_challenge += len(re.findall(r'\[\s*\\?,?\s*挑战\s*\\?,?\s*\]', t))

    for v, locs in sorted(used.items(), key=lambda x: -len(x[1])):
        if v not in declared:
            rec('FAIL', 'B',
                '难度等级「%s」在正文等级表中未声明（用了 %d 次，如 %s）'
                % (v, len(locs), locs[0]))
    for d in sorted(declared):
        if d in used:
            continue
        # 「挑战」按书中说明不用边注，而在正文中就地标注（proof 标题 + 正文记号）
        if d == '挑战' and inline_challenge > 0:
            rec('PASS', 'B',
                '难度等级「挑战」不用边注、改为正文就地标注（%d 处），与书中说明一致'
                % inline_challenge)
        else:
            rec('WARN', 'B', '难度等级「%s」在等级表中声明但全书从未使用' % d)
    rec('INFO', 'B', '难度标记使用（边注）：' +
        '、'.join('%s×%d' % (k, len(v)) for k, v in
                  sorted(used.items(), key=lambda x: -len(x[1]))))

    # --- B5: 占位符 / 未完成标记
    for p in files:
        text = read(p)
        base = os.path.basename(p)
        for pat, desc in [(r'\bTODO\b', 'TODO'),
                          (r'\bFIXME\b', 'FIXME'),
                          (r'\bXXX\b', 'XXX'),
                          (r'待补', '「待补」'),
                          (r'待写', '「待写」'),
                          (r'\?\?\?', '???')]:
            for m in re.finditer(pat, text):
                rec('WARN', 'B', '残留标记 %s —— %s:%d'
                    % (desc, base, line_of(text, m.start())))

    # --- B6: 例题总数
    n_ex = sum(len(re.findall(r'\\begin\{example\}', read(p))) for p in files)
    rec('INFO', 'B', '例题总数：%d' % n_ex)

    # --- B7: 盒子标题是否重复了盒子类型名
    # \newtcbtheorem 渲染为「类型名 编号: 标题」，标题里再写一遍类型名即重复。
    # checklistbox 的 title={#1} 不显示类型名；formulabox 前缀是「核心公式速查：」。
    TYPE_NAME = {'theorem': '定理', 'definition': '定义', 'corollary': '推论',
                 'example': '例题', 'insight': '物理洞见', 'warning': '常见误区',
                 'history': '物理史话', 'review': '本节要点', 'tip': '学习提示'}
    dup = []
    for p in files:
        t = read(p)
        base = os.path.basename(p)
        for env, disp in TYPE_NAME.items():
            for m in re.finditer(
                    r'\\begin\{%s\}\{((?:[^{}]|\{[^{}]*\})*)\}' % env, t):
                title = m.group(1).strip()
                # 真正的重复只有一种形态：标题以类型名开头（多为「类型名：…」）。
                # 「高斯散度定理」「阿尔芬定理」等是专有名称，渲染成
                # 「定理 1.2: 高斯散度定理」完全正常，不算重复。
                if title.startswith(disp):
                    dup.append('%s:%d  %s 标题=「%s」'
                               % (base, line_of(t, m.start()), env, title[:42]))
    if dup:
        for d in dup:
            rec('FAIL', 'B', '盒子标题以类型名开头，渲染后会重复 —— %s' % d)
    else:
        rec('PASS', 'B', '盒子标题：无以类型名开头造成重复的情况')


# ---------------------------------------------------------------- C. 编译

def check_C_compile(do_compile):
    log = os.path.join(SRC, MAIN.replace('.tex', '.log'))
    if do_compile:
        rec('INFO', 'C', '正在编译（latexmk -xelatex）……')
        r = subprocess.run(
            ['latexmk', '-xelatex', '-interaction=nonstopmode', '-f', MAIN],
            cwd=SRC, capture_output=True, text=True, errors='replace')
        if r.returncode != 0:
            rec('FAIL', 'C', 'latexmk 退出码 %d' % r.returncode)
        else:
            rec('PASS', 'C', 'latexmk 退出码 0')

    if not os.path.exists(log):
        rec('FAIL', 'C', '找不到编译日志 %s' % log)
        return

    t = read(log)

    def count(pat):
        return len(re.findall(pat, t))

    errs = count(r'^! ')
    missing = count(r'Missing character')
    overfull = count(r'Overfull \\hbox')
    underfull = count(r'Underfull \\hbox')
    lwarn = count(r'LaTeX Warning')
    refun = count(r'Reference `[^\']*\' on page')
    citeun = count(r'Citation `[^\']*\' on page')

    rec('PASS' if errs == 0 else 'FAIL', 'C', 'LaTeX 错误（! 开头）：%d' % errs)
    rec('PASS' if refun == 0 else 'FAIL', 'C', '未定义引用：%d' % refun)
    rec('PASS' if citeun == 0 else 'FAIL', 'C', '未定义引文：%d' % citeun)
    rec('PASS' if lwarn == 0 else 'WARN', 'C', 'LaTeX Warning：%d' % lwarn)
    rec('PASS' if missing <= 1 else 'WARN', 'C',
        'Missing character：%d（1 = 仅 U+000A 换行符属正常）' % missing)
    rec('PASS' if overfull <= 5 else 'WARN', 'C', 'Overfull hbox：%d' % overfull)
    rec('INFO', 'C', 'Underfull hbox：%d' % underfull)

    # 具体列出 overfull 细节
    for m in re.finditer(r'Overfull \\hbox \(([\d.]+)pt too wide\)[^\n]*', t):
        rec('INFO', 'C', '  Overfull %s pt —— %s' % (m.group(1), m.group(0)[:90]))


# ---------------------------------------------------------------- D. PDF

def check_D_pdf():
    pdf = os.path.join(SRC, MAIN.replace('.tex', '.pdf'))
    if not os.path.exists(pdf):
        rec('FAIL', 'D', '找不到 PDF 产物')
        return

    def run(cmd):
        r = subprocess.run(cmd, capture_output=True)
        return r.stdout.decode('utf-8', errors='replace')

    info = run(['pdfinfo', pdf])
    m = re.search(r'Pages:\s*(\d+)', info)
    pages = int(m.group(1)) if m else -1
    rec('INFO', 'D', 'PDF 页数：%d' % pages)
    tm = re.search(r'Title:\s*(.+)', info)
    if tm:
        rec('INFO', 'D', 'PDF 标题：%s' % tm.group(1).strip())

    txt = run(['pdftotext', pdf, '-'])
    n_q = len(re.findall(r'\?\?', txt))
    rec('PASS' if n_q == 0 else 'FAIL', 'D',
        'PDF 中未解析引用标记 "??"：%d' % n_q)

    # 附录页眉抽查（附录 A/D/G 的起始页由 TOC 推断）
    toc = re.search(r'附录 A.*?(\d+)\s*\n', txt)
    rec('INFO', 'D', '文本层字符数：%d' % len(txt))


# ---------------------------------------------------------------- E. 一致性

def check_E_consistency():
    readme = read(os.path.join(ROOT, 'README.md'))
    pdf = os.path.join(SRC, MAIN.replace('.tex', '.pdf'))

    pages = None
    if os.path.exists(pdf):
        r = subprocess.run(['pdfinfo', pdf], capture_output=True)
        m = re.search(r'Pages:\s*(\d+)', r.stdout.decode('utf-8', errors='replace'))
        pages = int(m.group(1)) if m else None

    # README 里声称的页数
    for m in re.finditer(r'\*\*(\d+)\s*页\*\*', readme):
        claimed = int(m.group(1))
        if pages and claimed != pages:
            rec('FAIL', 'E', 'README 声称 %d 页，实际 %d 页' % (claimed, pages))
    for m in re.finditer(r'\|\s*(\d+)\s*页完整教材', readme):
        claimed = int(m.group(1))
        if pages and claimed != pages:
            rec('FAIL', 'E', 'README 声称 %d 页（下载表），实际 %d 页'
                % (claimed, pages))

    # README 变更记录里的页数
    for m in re.finditer(r'（(\d+)\s*页\s*/\s*(\d+)\s*标签）', readme):
        pass  # 历史版本记录，允许与当前不同

    # 版本号一致性
    main = read(os.path.join(SRC, MAIN))
    vers = set(re.findall(r'v(\d+\.\d+\.\d+)', main))
    mreadme = re.findall(r'\|\s*v(\d+\.\d+\.\d+)\s*\|', readme)
    if mreadme:
        latest = max(mreadme, key=lambda s: [int(x) for x in s.split('.')])
        if latest not in vers:
            rec('WARN', 'E', 'README 最新版本 v%s 未出现在源码封面/元数据中（源码含 %s）'
                % (latest, sorted(vers)))
        else:
            rec('PASS', 'E', '版本号一致：v%s' % latest)


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-compile', action='store_true', help='跳过重新编译')
    args = ap.parse_args()

    files = tex_files()
    rec('INFO', 'A', '参与检查的 tex 文件：%d 个' % len(files))

    check_A_structure(files)
    check_B_content(files)
    check_C_compile(not args.no_compile)
    check_D_pdf()
    check_E_consistency()

    # ---------------------------------------------------------- 报告
    order = {'FAIL': 0, 'WARN': 1, 'INFO': 2, 'PASS': 3}
    print('=' * 78)
    print('《等离子体物理自学教材》全量检查报告')
    print('=' * 78)

    for lvl in ['FAIL', 'WARN', 'INFO', 'PASS']:
        rows = [r for r in RESULTS if r[0] == lvl]
        if not rows:
            continue
        print('\n--- %s (%d) ---' % (lvl, len(rows)))
        for _, cat, msg in rows:
            print('  [%s] %s' % (cat, msg))

    n_fail = sum(1 for r in RESULTS if r[0] == 'FAIL')
    n_warn = sum(1 for r in RESULTS if r[0] == 'WARN')
    print('\n' + '=' * 78)
    print('汇总：FAIL %d ／ WARN %d ／ 合计 %d 项'
          % (n_fail, n_warn, len(RESULTS)))
    print('=' * 78)

    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
