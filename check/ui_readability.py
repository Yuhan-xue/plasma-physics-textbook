#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》可读性与版式（UI）审核
================================================

与 full_check.py 互补：full_check 查「对不对」，本模块查「好不好读、好不好看」。

  F. 可读性  —— 句长/段长分布、超长句段、散文中埋藏的枚举、标题层级、
                连续公式墙、术语一致性
  G. 版式 UI —— 颜色对比度（WCAG 计算）、盒子嵌套、浮动体与引用距离、
                边注溢出、页眉一致性、表格超宽

用法：
    python check/ui_readability.py
    python check/ui_readability.py --json      # 机器可读输出

退出码：0 = 无 FAIL；1 = 存在 FAIL
"""

import io
import os
import re
import sys
import glob
import json
import argparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
MAIN = 'plasma_physics_textbook_v2.tex'

RESULTS = []


def rec(level, cat, msg):
    RESULTS.append((level, cat, msg))


def read(p):
    return io.open(p, encoding='utf-8').read()


def strip_comments(t):
    out = []
    for line in t.split('\n'):
        i, idx = 0, None
        while i < len(line):
            if line[i] == '\\':
                i += 2
                continue
            if line[i] == '%':
                idx = i
                break
            i += 1
        out.append(line if idx is None else line[:idx])
    return '\n'.join(out)


def tex_files():
    return [p for p in sorted(glob.glob(os.path.join(SRC, '*.tex')))
            if os.path.basename(p) != 'style_v2_preview.tex']


def line_of(text, idx):
    return text.count('\n', 0, idx) + 1


# ---------------------------------------------------------------- 文本提取

DROP_ENVS = ['tabular', 'tabularx', 'tikzpicture', 'equation', 'equation*',
             'align', 'align*', 'gather', 'gather*', 'verbatim', 'Verbatim',
             'lstlisting', 'thebibliography', 'array', 'matrix', 'bmatrix',
             'pmatrix', 'vmatrix', 'cases', 'split', 'aligned']

# 「结构化小结」环境：其条目本就是「把整节压缩成一行」的密集摘要，
# 不是供逐句阅读的散文。全书 13 个 chapterreview 的 57 个条目**无一**使用
# 嵌套列表（中位 63 字，p90 173 字），说明密集单行是本书既有版式约定。
# 因此对它们单独统计，不与正文散文混用同一套句长阈值。
STRUCTURED_ENVS = ['chapterreview', 'checklistbox', 'checklist']


def _drop_env(t, env):
    out, i = [], 0
    b, e = r'\begin{%s}' % env, r'\end{%s}' % env
    while True:
        j = t.find(b, i)
        if j < 0:
            out.append(t[i:])
            break
        out.append(t[i:j])
        k = t.find(e, j)
        if k < 0:
            break
        i = k + len(e)
    return ''.join(out)


def strip_citations(t):
    """剔除文献条目。

    附录 E 的书单条目形如「\\item \\textbf{作者, \\textit{书名}, 出版社, 年份.}」，
    它们是格式化引用而非散文，计入句长统计会制造假的长句。
    """
    out = []
    for line in t.split('\n'):
        s = line.strip()
        if s.startswith(r'\item') and re.search(
                r'\\textit\{|\\emph\{', s) and re.search(
                r'(19|20)\d{2}', s) and re.search(
                r'Press|Wiley|Springer|Publishers?|University|Verlag|'
                r'Academic|Publishing|Inc\.|Ltd\.|DOI', s):
            out.append('')
            continue
        out.append(line)
    return '\n'.join(out)


def prose(t, drop_preamble=False, drop_structured=False):
    """把 LaTeX 源码还原成接近排版的散文，供可读性统计。

    关键：必须先剔除行间公式（\\[...\\]、$$...$$），否则公式会被当成超长句子。

    drop_structured=True 时再剔除 chapterreview/checklist 等「结构化小结」，
    用于单独统计正文散文——小结条目密集是版式约定，不该按散文句长判定。
    """
    if drop_preamble:
        j = t.find(r'\begin{document}')
        if j >= 0:
            t = t[j:]
    if drop_structured:
        for env in STRUCTURED_ENVS:
            for _ in range(3):
                if r'\begin{%s}' % env in t:
                    t = _drop_env(t, env)
    t = re.sub(r'\\\[.*?\\\]', '\n', t, flags=re.S)
    t = re.sub(r'\$\$.*?\$\$', '\n', t, flags=re.S)
    t = re.sub(r'\\\(.*?\\\)', ' ', t, flags=re.S)
    t = re.sub(r'\$[^$]*\$', ' ', t)
    for env in DROP_ENVS:
        for _ in range(3):
            if r'\begin{%s}' % env in t:
                t = _drop_env(t, env)
    t = strip_citations(t)
    t = t.replace('\\\\', '\n')
    t = re.sub(r'\\begin\{[^}]*\}(\[[^\]]*\])?(\{[^{}]*\})*', '\n', t)
    t = re.sub(r'\\end\{[^}]*\}', '\n', t)
    t = re.sub(r'\\(label|ref|cref|Cref|eqref|cite|index|addcontentsline|'
               r'includegraphics|difficulty|difficultyhere|epigraph|'
               r'tcblower|centering|toprule|midrule|bottomrule|hline)'
               r'(\[[^\]]*\])?(\{[^{}]*\})*', ' ', t)
    t = re.sub(r'\\(section|subsection|subsubsection|chapter|item|caption)\*?'
               r'(\[[^\]]*\])?', '\n', t)
    t = re.sub(r'\\[a-zA-Z@]+\*?(\[[^\]]*\])?', ' ', t)
    # 组内设置（如 {\emergencystretch=3em ...}）：命令名已被上一步删掉，
    # 残留的 "=3em" 是解析产物而非正文，必须一并清掉。
    t = re.sub(r'=\s*[\d.]+\s*(?:em|ex|pt|pc|cm|mm|in|sp|bp|dd|cc)\b', ' ', t)
    t = re.sub(r'[{}$&\\]', ' ', t)
    return t


def paragraphs(body):
    out, cur = [], []
    for line in body.split('\n'):
        s = re.sub(r'\s+', ' ', line).strip()
        if s:
            cur.append(s)
        else:
            if cur:
                out.append(''.join(cur))
                cur = []
    if cur:
        out.append(''.join(cur))
    return [p for p in out if len(p) >= 24]


def sentences(par):
    return [s.strip() for s in re.split(r'[。！？]', par) if len(s.strip()) >= 10]


# ================================================================ F 可读性

def check_F_readability(files):
    all_s, all_p = [], []
    summ_s = []
    for p in files:
        raw = strip_comments(read(p))
        body = prose(raw, drop_preamble=True, drop_structured=True)
        ps = paragraphs(body)
        ss = []
        for x in ps:
            ss += sentences(x)
        all_p += ps
        all_s += ss
        # 结构化小结（本章要点回顾/自查清单）单独取，单独定标准
        for env in STRUCTURED_ENVS:
            for m in re.finditer(r'(?s)\\begin\{%s\}(.*?)\\end\{%s\}'
                                 % (env, env), raw):
                for x in paragraphs(prose(m.group(1))):
                    summ_s += sentences(x)

    def stats(v):
        v = sorted(len(x) for x in v)
        n = len(v)
        if not n:
            return {}
        return {'n': n,
                'p50': v[n // 2], 'p75': v[int(n * .75)], 'p90': v[int(n * .9)],
                'p95': v[int(n * .95)], 'p99': v[min(int(n * .99), n - 1)],
                'max': v[-1]}

    ss, ps = stats(all_s), stats(all_p)
    rec('INFO', 'F', '正文句长(字)：中位 %d，p75 %d，p90 %d，p95 %d，p99 %d，最长 %d'
        % (ss['p50'], ss['p75'], ss['p90'], ss['p95'], ss['p99'], ss['max']))
    rec('INFO', 'F', '正文段长(字)：中位 %d，p75 %d，p90 %d，p95 %d，p99 %d，最长 %d'
        % (ps['p50'], ps['p75'], ps['p90'], ps['p95'], ps['p99'], ps['max']))

    if summ_s:
        ms = stats(summ_s)
        rec('INFO', 'F', '要点小结条目(字)：中位 %d，p90 %d，最长 %d（%d 条，'
            '密集单行是本书版式约定）'
            % (ms['p50'], ms['p90'], ms['max'], ms['n']))

    for th, lvl in ((100, 'WARN'), (60, 'INFO')):
        c = sum(1 for x in all_s if len(x) > th)
        rec(lvl, 'F', '正文句长 >%d 字：%d 句（%.1f%%）'
            % (th, c, 100.0 * c / len(all_s)))

    # 可读性总评（中文技术教材的经验区间：p90 ≤ 60 属良好）
    if ss['p90'] <= 60:
        rec('PASS', 'F', '正文句长良好（p90=%d ≤ 60 字），适合自学阅读' % ss['p90'])
    elif ss['p90'] <= 80:
        rec('WARN', 'F', '正文句长偏长（p90=%d），部分句子建议拆分' % ss['p90'])
    else:
        rec('FAIL', 'F', '正文句长过长（p90=%d），显著影响自学可读性' % ss['p90'])

    if ps['p90'] <= 200:
        rec('PASS', 'F', '段长总体良好（p90=%d 字）' % ps['p90'])
    else:
        rec('WARN', 'F', '段长偏长（p90=%d 字），建议拆段' % ps['p90'])

    # ---- 超长句 / 超长段 清单
    long_s = [(len(s), i) for i, s in enumerate(all_s) if len(s) > 100]
    for ln, i in sorted(long_s, reverse=True)[:12]:
        rec('WARN', 'F', '超长句 %d 字：%s…' % (ln, all_s[i][:70]))
    long_p = [(len(x), x) for x in all_p if len(x) > 300]
    for ln, x in sorted(long_p, reverse=True)[:8]:
        rec('WARN', 'F', '超长段 %d 字：%s…' % (ln, x[:70]))

    # ---- 散文中埋藏的枚举（应改为列表）
    buried = 0
    for p in files:
        body = prose(strip_comments(read(p)), drop_preamble=True)
        for x in paragraphs(body):
            for s in sentences(x):
                marks = len(re.findall(r'\(i+\)|\([1-9]\)|（[1-9]）|'
                                       r'[①②③④⑤]', s))
                if marks >= 3 and len(s) > 80:
                    buried += 1
                    rec('WARN', 'F', '单句内埋了 %d 项枚举，建议改为列表：%s…'
                        % (marks, s[:70]))
    if buried == 0:
        rec('PASS', 'F', '未发现「单句内埋多项枚举」的情况')

    # ---- 标题层级跳跃
    main = strip_comments(read(os.path.join(SRC, MAIN)))
    levels, skips = [], []
    for m in re.finditer(r'\\(chapter|section|subsection|subsubsection)\*?\{',
                         main):
        lv = {'chapter': 0, 'section': 1, 'subsection': 2,
              'subsubsection': 3}[m.group(1)]
        if levels and lv > levels[-1] + 1:
            skips.append((line_of(main, m.start()), levels[-1], lv))
        levels.append(lv)
    if skips:
        for ln, a, b in skips:
            rec('FAIL', 'F', '标题层级跳跃（%s → %s）—— 主文件:%d'
                % (a, b, ln))
    else:
        rec('PASS', 'F', '标题层级：无跳级')

    # ---- 连续公式墙（≥3 个行间公式之间无解释文字）
    walls = 0
    for p in files:
        txt = strip_comments(read(p))
        lines = txt.split('\n')
        run, start = 0, 0
        for i, ln in enumerate(lines):
            s = ln.strip()
            is_eq = s.startswith('\\[') or s.startswith('\\begin{equation') \
                or s.startswith('\\begin{align') or s.startswith('$$')
            istxt = bool(s) and not is_eq and not s.startswith('\\end{') \
                and len(re.sub(r'\\[a-zA-Z]+(\[[^\]]*\])?', '', s).strip()) > 12
            if is_eq:
                if run == 0:
                    start = i + 1
                run += 1
            elif istxt:
                if run >= 3:
                    walls += 1
                    rec('WARN', 'F', '连续 %d 个行间公式无解释文字 —— %s:%d'
                        % (run, os.path.basename(p), start))
                run = 0
    if walls == 0:
        rec('PASS', 'F', '未发现连续公式墙（公式间均有解释文字）')

    return None


# ================================================================ G 版式 UI

def _srgb(c):
    c = c / 255.0
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def luminance(rgb):
    r, g, b = (_srgb(x) for x in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(rgb1, rgb2):
    l1, l2 = luminance(rgb1), luminance(rgb2)
    hi, lo = max(l1, l2), min(l1, l2)
    return (hi + 0.05) / (lo + 0.05)


def blend(rgb, pct):
    """tcolorbox 的 color!pct!black：按 pct% 与原色混合黑。"""
    f = pct / 100.0
    return tuple(int(round(c * f)) for c in rgb)


def check_G_ui(files):
    main = strip_comments(read(os.path.join(SRC, MAIN)))

    # ---- 颜色定义
    colors = {}
    for m in re.finditer(r'\\definecolor\{(\w+)\}\{RGB\}\{(\d+)\s*,\s*(\d+)\s*,'
                         r'\s*(\d+)\}', main):
        colors[m.group(1)] = tuple(int(m.group(i)) for i in (2, 3, 4))

    # ---- 各语义盒的边框色（= 标题栏底色）与标题文字色
    #      \newtcbtheorem{env}{名称}{colback=..,colframe=X,..}{prefix}
    boxes = []
    for m in re.finditer(r'\\newtcbtheorem(?:\[[^\]]*\])?\{(\w+)\}\{([^}]*)\}'
                         r'\{(.*?)\}\{(\w+)\}', main, flags=re.S):
        env, disp, opts = m.group(1), m.group(2), m.group(3)
        cf = re.search(r'colframe\s*=\s*([\w!]+)', opts)
        cbt = re.search(r'colbacktitle\s*=\s*([\w!]+)', opts)
        ct = re.search(r'coltitle\s*=\s*([\w!]+)', opts)
        boxes.append((env, disp, cf.group(1) if cf else None,
                      cbt.group(1) if cbt else None,
                      ct.group(1) if ct else None))

    def resolve(spec):
        """把 'warnyellow' / 'insightred!80' / 'warnyellow!60!black' 解析成 RGB。"""
        if not spec:
            return None
        parts = spec.split('!')
        base = colors.get(parts[0])
        if base is None:
            return None
        if len(parts) == 1:
            return base
        if len(parts) == 2:
            return blend(base, float(parts[1]))
        if len(parts) == 3 and parts[2] == 'black':
            return blend(base, float(parts[1]))
        if len(parts) == 3 and parts[2] == 'white':
            f = float(parts[1]) / 100.0
            return tuple(int(round(c * f + 255 * (1 - f))) for c in base)
        return base

    WHITE = (255, 255, 255)
    fails, warns, passes = [], [], []
    for env, disp, cf, cbt, ct in boxes:
        bg = resolve(cbt) or resolve(cf)
        if bg is None:
            continue
        fg = resolve(ct) or WHITE      # tcolorbox 默认 coltitle=white
        cr = contrast(fg, bg)
        row = '%-12s %-8s 标题栏对比度 %.2f:1' % (env, disp, cr)
        if cr < 3.0:
            fails.append(row + '  ← 低于 WCAG AA 大字号下限 3:1')
        elif cr < 4.5:
            warns.append(row + '  ← 低于 WCAG AA 正文下限 4.5:1（大字达标）')
        else:
            passes.append(row)

    for r in fails:
        rec('FAIL', 'G', r)
    for r in warns:
        rec('WARN', 'G', r)
    rec('PASS' if not fails else 'INFO', 'G',
        '颜色对比度：%d 个盒达标（≥4.5:1）、%d 个仅达大字标准、%d 个不达标'
        % (len(passes), len(warns), len(fails)))

    # ---- 盒子嵌套
    nest = 0
    BOX = list({b[0] for b in boxes} |
               {'miniquiz', 'chapterreview', 'formula', 'formulabox',
                'checklist', 'checklistbox', 'chapterintro'})
    for p in files:
        txt = strip_comments(read(p))
        stack = []
        for m in re.finditer(r'\\(begin|end)\{([^}]+)\}', txt):
            k, env = m.group(1), m.group(2)
            if env not in BOX:
                continue
            if k == 'begin':
                if stack:
                    nest += 1
                    rec('FAIL', 'G', '盒子嵌套：%s 内又开 %s —— %s:%d'
                        % (stack[-1], env, os.path.basename(p),
                           line_of(txt, m.start())))
                stack.append(env)
            elif stack and stack[-1] == env:
                stack.pop()
    if nest == 0:
        rec('PASS', 'G', '盒子嵌套：无语义盒相互嵌套')

    # ---- 语义盒标题缺失
    # tcolorbox 在标题为空时**不渲染冒号**，PDF 里会直接显示成「物理洞见 11.5」
    # 再接正文，读者看不到任何标题。空 key 无害（只是不生成 \label），只查标题。
    SEM = {'insight': '物理洞见', 'warning': '常见误区', 'history': '物理史话',
           'definition': '定义', 'theorem': '定理', 'corollary': '推论',
           'example': '例题', 'tip': '学习提示', 'review': '本节要点'}
    noTitle = []
    for p in files:
        txt = strip_comments(read(p))
        for env, cn in SEM.items():
            for m in re.finditer(r'\\begin\{%s\}\{([^}]*)\}' % env, txt):
                if m.group(1).strip() == '':
                    noTitle.append((os.path.basename(p),
                                    line_of(txt, m.start()), cn))
    if noTitle:
        for fn, ln, cn in noTitle[:15]:
            rec('FAIL', 'G', '语义盒标题为空（PDF 中只见「%s 编号」）—— %s:%d'
                % (cn, fn, ln))
        if len(noTitle) > 15:
            rec('FAIL', 'G', '……另有 %d 处标题为空' % (len(noTitle) - 15))
    else:
        rec('PASS', 'G', '语义盒标题：全部非空（%d 类盒均已检查）' % len(SEM))

    # ---- 浮动体与首次引用的距离（图）
    fig_pos = {}
    for m in re.finditer(r'\\begin\{figure\}', main):
        end = main.find(r'\end{figure}', m.end())
        body = main[m.start():end]
        lm = re.search(r'\\label\{([^}]+)\}', body)
        if lm:
            fig_pos[lm.group(1)] = line_of(main, m.start())
    far = 0
    for lab, fline in fig_pos.items():
        for m in re.finditer(r'\\(?:c?ref|autoref)\{%s\}' % re.escape(lab), main):
            d = abs(line_of(main, m.start()) - fline)
            if d > 250:
                far += 1
                rec('WARN', 'G', '图 %s 与其引用相距 %d 行（引用在 %d，图在 %d）'
                    % (lab, d, line_of(main, m.start()), fline))
            break
    if far == 0 and fig_pos:
        rec('PASS', 'G', '浮动体距离：所有图的首次引用均在 250 行内（同章可见）')

    # ---- 边注溢出（difficulty marginnote）
    log = os.path.join(SRC, MAIN.replace('.tex', '.log'))
    if os.path.exists(log):
        lt = read(log)
        mo = len(re.findall(r'Marginpar.*too (?:wide|high)', lt))
        ov = len(re.findall(r'Overfull \\hbox', lt))
        un = len(re.findall(r'Underfull \\hbox', lt))
        rec('PASS' if mo == 0 else 'FAIL', 'G',
            '边注（难度标记）溢出：%d' % mo)
        rec('PASS' if ov <= 5 else 'WARN', 'G', 'Overfull hbox：%d' % ov)
        rec('INFO', 'G', 'Underfull hbox：%d' % un)
        # 宽表格导致的 overfull
        for m in re.finditer(r'Overfull \\hbox \(([\d.]+)pt too wide\)'
                             r'[^\n]*', lt):
            rec('INFO', 'G', '  溢出 %s pt —— %s' % (m.group(1), m.group(0)[:80]))

    # ---- 页眉一致性（从 PDF 检查附录页眉）
    pdf = os.path.join(SRC, MAIN.replace('.tex', '.pdf'))
    if os.path.exists(pdf):
        import subprocess
        r = subprocess.run(['pdftotext', '-layout', pdf, '-'],
                           capture_output=True)
        txt = r.stdout.decode('utf-8', errors='replace')
        ap = len(re.findall(r'附录 [A-G]', txt))
        rec('PASS' if ap > 0 else 'WARN', 'G',
            '附录页眉：检测到 %d 处「附录 X」标记' % ap)


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    files = tex_files()
    check_F_readability(files)
    check_G_ui(files)

    if args.json:
        print(json.dumps([{'level': a, 'cat': b, 'msg': c}
                          for a, b, c in RESULTS],
                         ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》可读性与版式（UI）审核')
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
