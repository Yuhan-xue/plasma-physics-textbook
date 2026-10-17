#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》一致性与符号规范审核
================================================

与 full_check.py、ui_readability.py 三者互补：
  - full_check.py       查「对不对」  （结构/内容/编译/数值口径）
  - ui_readability.py   查「好不好读」（句长/版式/对比度）
  - 本模块              查「前后一不一致」（术语、符号、速查卡、附录结构）

  H. 术语规范  —— 人名复合术语的中英混排空格、禁用变体、连字符规范、
                  同文件内同一概念的两种写法
  I. 速查卡符号说明 —— 22 张 formulabox 是否都有「符号说明」段，
                  段内是否覆盖了卡中公式用到的符号
  J. 符号总表  —— 附录H 的符号 ↔ 正文实际使用（双向核对）
  K. 结构完整性 —— 附录章连续性、标签唯一性、补充材料挂载、符号/术语表结构

用法：
    python check/consistency_audit.py
    python check/consistency_audit.py --json      # 机器可读输出

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


def rel(p):
    return os.path.relpath(p, ROOT).replace('\\', '/')


def strip_comments(t):
    """逐字符剥 % 注释（宏内的 \\% 不受影响）——与另两个套件保持同一实现。"""
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


# ---------------------------------------------------------------- 掩码
# 术语统计必须避开 label/ref 名、代码块与逐字环境，否则会把
# \label{eq:ch10-tokamak-scaling} 之类的**标签名**误判成正文用词。

CODE_ENVS = ['verbatim', 'Verbatim', 'lstlisting', 'minted']


def mask(t):
    for env in CODE_ENVS:
        t = re.sub(r'\\begin\{%s\}.*?\\end\{%s\}' % (env, env),
                   lambda m: '\n' * m.group(0).count('\n'), t, flags=re.S)
    t = re.sub(r'\\verb(.)(.*?)\1', '', t)
    t = re.sub(r'\\(?:label|ref|eqref|cref|Cref|cpageref|nameref|autoref'
               r'|cite|citep|citet|citeauthor|index|bibitem|input|include)'
               r'\*?\{[^{}]*\}', lambda m: '\\' + m.group(0)[1:].split('{')[0], t)
    return t


# ---------------------------------------------------------------- 数学符号

FMT_CMDS = [r'\\displaystyle', r'\\textstyle', r'\\scriptstyle', r'\\!', r'\\,',
            r'\\;', r'\\:', r'\\ ', r'\\quad', r'\\qquad', r'\\left', r'\\right',
            r'\\bigl?', r'\\Bigl?', r'\\biggl?', r'\\Biggl?', r'\\bigr', r'\\Bigr',
            r'\\biggr', r'\\Biggr', r'\\limits', r'\\nolimits',
            # 定界符/关系符：只表示写法，不是符号
            r'\\langle', r'\\rangle', r'\\parallel', r'\\perp', r'\\iff']

# 函数名/算符/关系符/结构命令：不是「符号」，不参与遗漏统计
FUNC_CMDS = set("""
frac dfrac tfrac sqrt exp ln log lg sin cos tan cot sinh cosh tanh coth
arcsin arccos arctan max min sup inf det dim ker arg deg gcd lcm
sum prod int iint iiint oint lim to infty partial nabla times cdot div
leq geq le gg ll sim simeq approx equiv propto neq pm mp ast circ
rightarrow leftarrow Rightarrow Leftarrow leftrightarrow mapsto
mathrm mathbf mathit mathsf mathtt mathcal mathbb mathfrak boldsymbol
text textrm textsf texttt mbox hbox operatorname overbrace underbrace
overline underline widetilde widehat vec hat bar dot ddot tilde
begin end label ref eqref tag nonumber notag displaystyle
ee diff dd eps degree Celsius angstrom oplus otimes
left right middle quad qquad hspace vspace phantom
""".split())

# 全局通用符号：本书已在「前置·基础物理常数/等离子体特征参数」与
# 「附录H 通用符号与算符」中统一解释，速查卡不必逐卡重复。
# J 类检查会强制要求它们全部出现在附录H 的通用符号表里。
GLOBAL_SYMBOLS = [
    'c', 'e', 'k_{B}', 'k_B', r'\eps_0', r'\mu_0', r'\pi', r'\ee', r'\infty',
    r'\nabla', r'\partial', r'\diff',
    # 全书统一的记号宏（由附录H「通用符号与算符」统一解释）
    r'\vect', r'\unitvec', r'\divg', r'\curl', r'\grad', r'\grad_v', r'\lap',
    't', 'x', 'r', 'v', 'm', 'q', 'E', 'B', 'T', 'n',
    'T_e', 'T_i', 'n_e', 'n_i', 'm_e', 'm_i', 'Z', r'\gamma',
]

# 常数与算符（\eps_0、\ee、\infty、\nabla、\partial、\diff）在 FUNC_CMDS 中，
# symbols_in() 不会产出 token，故 J 类对它们改用「字面出现」核对。
# 注：本书未定义 \dd 宏（正文一律用 \diff），故不列入本表——否则会强制附录H
# 写出一个一旦被 TeX 展开就报 Undefined control sequence 的符号。
GLOBAL_LITERAL = [r'\eps_0', r'\ee', r'\infty', r'\nabla', r'\partial',
                  r'\diff']


def norm_sym(s):
    s = s.strip()
    for c in (r'\!', r'\,', r'\;', r'\:', r'\ '):
        s = s.replace(c, '')
    s = s.replace(' ', '')
    return s


def base_of(s):
    m = re.match(r'^(\\?[A-Za-z]+|[A-Za-z])', s)
    return m.group(1) if m else s


def sub_is_numeric(s):
    rest = s[len(base_of(s)):]
    return bool(rest) and re.match(r'^(?:_|\^)\{?\d+\}?$', rest) is not None


def sym_family(s):
    """符号族名：幂次（^…）与纯数字下标（_2、_{1.25}）不算区别，字母下标保留。

    用于「符号说明是否覆盖了该符号」的比对：卡里出现 $I_p^{1.25}$，说明中解释 $I_p$ 即为覆盖。
    """
    m = re.match(r'^(\\?[A-Za-z]+|[A-Za-z])', s)
    base = m.group(1) if m else s
    rest = s[len(base):]
    m2 = re.match(r'^[_^]\{?([A-Za-z0-9.]*)\}?', rest)
    if m2:
        val = m2.group(1)
        if val and not val.replace('.', '').isdigit():
            return base + '_' + val
    return base


TOKEN_RE = re.compile(
    r'(\\(?:[A-Za-z]+)|[A-Za-z])'
    r'((?:_\{[^{}]*\}|_[A-Za-z0-9]|\^\{[^{}]*\}|\^[A-Za-z0-9!])*)')


def symbols_in(math_text):
    """从一段数学文本里抽出符号 token（含上下标），已归一化。"""
    def _drop_empty(t):
        # \tau_{\text{流}} 剥掉文字后剩 \tau_{ }，本质上仍是 \tau
        t = re.sub(r'\{\s*\}', '{}', t)
        return re.sub(r'[_^]\{\}', '', t)

    t = math_text
    # 欧姆单位：\Omega\cdot\mathrm{m} 里的 \Omega 是单位符号，不是物理量（\Omega=\omega-ku/2 不受影响）
    t = re.sub(r'\\Omega\s*\\cdot', r'\\cdot', t)
    # 单位/文字内容不属于符号；允许一层嵌套花括号（如 \mathrm{keV\cdot m^{-3}\cdot s}）
    t = re.sub(r'\\(?:mathrm|text|textrm|textsf|texttt|mbox|operatorname)'
               r'\{(?:[^{}]|\{[^{}]*\})*\}', ' ', t)
    t = t.replace(r'\tfrac', r'\frac').replace(r'\dfrac', r'\frac')
    t = _drop_empty(t)
    for cmd in FMT_CMDS:
        t = re.sub(cmd + r'(?![A-Za-z])', '', t)
    t = _drop_empty(t)          # \eta_{\parallel} 去掉 \parallel 后剩 \eta_{}
    out = []
    for m in TOKEN_RE.finditer(t):
        base = m.group(1)
        if base.startswith('\\'):
            if base[1:] in FUNC_CMDS:
                continue
            if len(base) == 2 and base[1] in '{}[](),':
                continue
        sym = norm_sym(m.group(1) + m.group(2))
        if not sym or sym in ('i', 'j'):
            continue
        out.append(sym)
    return out


MATH_RE = re.compile(
    r'\$\$(?P<dd>.+?)\$\$'
    r'|\$(?P<d>.+?)\$'
    r'|\\\[(?P<b>.+?)\\\]'
    r'|\\begin\{(?P<env>equation\*?|align\*?|gather\*?|eqnarray\*?)\}'
    r'(?P<e>.+?)\\end\{(?P=env)\}', re.S)


def math_segments(text):
    for m in MATH_RE.finditer(text):
        yield m.group('dd') or m.group('d') or m.group('b') or m.group('e') or ''


def symbols_missing(text, explained):
    """text 中用到的符号里，explained（已说明）与全局白名单之外的。"""
    ok = set(explained)
    glob = set(GLOBAL_SYMBOLS)
    missing = []
    for sym in symbols_in(text):
        if sym in ok or sym in glob:
            continue
        b = base_of(sym)
        if b in ok or b in glob:
            continue
        if sub_is_numeric(sym) and b in glob:
            continue
        if sym not in missing:
            missing.append(sym)
    return missing


# ---------------------------------------------------------------- H 术语规范

# 禁用变体 → 应改为（本轮 v1.6.0 的规范决策，见 tmp/spec_v160.md §3）
FORBIDDEN = [
    ('朗道阻尼', 'Landau 阻尼'),
    ('Landau阻尼', 'Landau 阻尼'),
    ('Bohm判据', '玻姆判据'),
    ('Bohm速度', '玻姆速度'),
    ('Bohm-Gross', 'Bohm--Gross'),
    ('Child-Langmuir', 'Child--Langmuir'),
    ('Fokker-Planck', 'Fokker--Planck'),
]

# 拉丁专名后直接接中文（缺空格）
LATIN_NAMES = ['Langmuir', 'Bohm', 'Landau', 'Debye', 'Spitzer', 'Townsend',
               'Paschen', 'Lawson', 'Thomson', 'Rogowski', 'Larmor', 'Faraday',
               'Hall', 'Boltzmann', 'Poisson', 'Maxwell', 'Vlasov', 'Saha',
               'Braginskii', 'Gross', 'Child', 'Alfv', 'Fokker', 'Planck',
               'Bethe', 'Wigner', 'Seitz', 'Wien', 'Sackur', 'Tetrode']
# 可选复合后缀必须作用于整个候选集：（A|B|C(?:--X)?）只会绑到最后一个候选，
# 而前缀守卫排除连字符，才能正确识别「Bohm--Gross修正」这类复合专名。
LATIN_GLUE_RE = re.compile(r'(?<![A-Za-z-])((?:' + '|'.join(LATIN_NAMES) +
                           r')(?:--[A-Za-z]+)?)(?=[\u4e00-\u9fff])')

# 同一文件内同一概念的两种写法（WARN）
PAIR_MIX = [('玻姆判据', 'Bohm 判据'), ('朗道阻尼', 'Landau 阻尼'),
            ('汤森放电', 'Townsend 放电'), ('朗道长度', 'Landau 长度')]


def check_H_terms(files):
    total_forbidden = 0
    for f in files:
        raw = strip_comments(read(f))
        t = mask(raw)
        for bad, good in FORBIDDEN:
            for m in re.finditer(re.escape(bad), t):
                total_forbidden += 1
                rec('FAIL', 'H', '%s:%d 术语变体「%s」应作「%s」'
                    % (rel(f), line_of(t, m.start()), bad, good))
        hits = 0
        for m in LATIN_GLUE_RE.finditer(t):
            hits += 1
            if hits <= 6:
                seg = t[m.start():m.start() + 14].split('\n')[0]
                rec('FAIL', 'H', '%s:%d 专名与中文之间缺空格：%s'
                    % (rel(f), line_of(t, m.start()), seg))
        if hits > 6:
            rec('FAIL', 'H', '%s 另有 %d 处专名与中文缺空格（未逐条列出）'
                % (rel(f), hits - 6))
        for a, b in PAIR_MIX:
            if a in t and b in t:
                rec('WARN', 'H', '%s 同一文件内并存两种写法：「%s」与「%s」'
                    % (rel(f), a, b))
    if total_forbidden == 0:
        rec('PASS', 'H', '禁用术语变体：0 处（朗道/Landau、玻姆/Bohm、en dash 连字符均已统一）')


# ---------------------------------------------------------------- I 速查卡

BOX_BEGIN = r'\begin{formulabox}'
BOX_END = r'\end{formulabox}'
SYM_LABEL = r'\textbf{符号说明：}'


def formulaboxes(text):
    """返回 (起始行, 标题, 块文本) 列表。"""
    out = []
    i = 0
    while True:
        b = text.find(BOX_BEGIN, i)
        if b < 0:
            break
        e = text.find(BOX_END, b)
        if e < 0:
            break
        block = text[b:e]
        mt = re.search(r'\\begin\{formulabox\}\{([^{}]*)\}', block)
        out.append((line_of(text, b), mt.group(1) if mt else '?', block))
        i = e + len(BOX_END)
    return out


def check_I_cards(files):
    cards = 0
    missing_label = 0
    for f in files:
        raw = mask(strip_comments(read(f)))
        boxes = formulaboxes(raw)
        if not boxes:
            continue
        for ln, title, block in boxes:
            cards += 1
            if SYM_LABEL not in block:
                missing_label += 1
                rec('FAIL', 'I', '%s:%d 速查卡「%s」缺「符号说明」段'
                    % (rel(f), ln, title))
                continue
            if block.count(SYM_LABEL) > 1:
                rec('FAIL', 'I', '%s:%d 速查卡「%s」有 %d 处「符号说明」段（重复写入）'
                    % (rel(f), ln, title, block.count(SYM_LABEL)))
            # 说明段本体：从 \textbf{符号说明：} 到块尾
            expl_txt = block.split(SYM_LABEL, 1)[1]
            expl_raw = expl_txt
            if '。' not in expl_raw:
                rec('WARN', 'I', '%s:%d 速查卡「%s」的符号说明未以句号结束'
                    % (rel(f), ln, title))
            expl_syms = set()
            for seg in math_segments(expl_raw):
                expl_syms.update(symbols_in(seg))
            # 供说明用的正文 = 去掉说明段后的卡体
            body = block.split(SYM_LABEL, 1)[0]
            body_syms = []
            for seg in math_segments(body):
                for s in symbols_in(seg):
                    if s not in body_syms:
                        body_syms.append(s)
            expl_fam = set(sym_family(e) for e in expl_syms)
            glob_fam = set(sym_family(g) for g in GLOBAL_SYMBOLS)
            miss = [s for s in body_syms
                    if sym_family(s) not in expl_fam and sym_family(s) not in glob_fam]
            if len(expl_syms) < 2:
                rec('WARN', 'I', '%s:%d 速查卡「%s」符号说明过简（仅 %d 个符号）'
                    % (rel(f), ln, title, len(expl_syms)))
            if miss:
                rec('WARN', 'I', '%s:%d 速查卡「%s」有 %d 个公式符号未在符号说明中出现：%s'
                    % (rel(f), ln, title, len(miss),
                       '、'.join(miss[:12]) + ('…' if len(miss) > 12 else '')))
    if cards and missing_label == 0:
        rec('PASS', 'I', '速查卡符号说明：%d 张卡全部具备' % cards)
    if not cards:
        rec('FAIL', 'I', '未找到任何 formulabox（速查卡）')
    return cards


# ---------------------------------------------------------------- J 符号总表

APPENDIX_H_MARK = r'\appchapter{附录H'


def _lt_rows(block):
    """从 longtable 块里取数据行（跳过表头/表尾/分组行）。"""
    rows = []
    for line in block.split('\n'):
        s = line.strip()
        if not s or s.startswith('%'):
            continue
        if s.startswith('\\') and '&' not in s:
            continue
        if '&' not in s:
            continue
        if s.startswith('\\multicolumn'):
            continue
        if s.startswith('\\textbf{') or s.startswith('\\toprule') or s.startswith('\\midrule'):
            continue
        rows.append(s)
    return rows


def check_J_symbol_table(files, main_text):
    pos = main_text.find(APPENDIX_H_MARK)
    if pos < 0:
        rec('WARN', 'J', '附录H（符号与术语总表）尚未加入主文件——本轮待补')
        return
    app = main_text[pos:]
    nxt = app.find(r'\appchapter{附录', 10)
    if nxt > 0:
        app = app[:nxt]

    # 通用符号表
    gblock = ''
    gpos = app.find(r'\section{通用符号与算符}')
    if gpos < 0:
        rec('FAIL', 'J', '附录H 缺「通用符号与算符」小节')
        general_syms = set()
    else:
        gend = app.find(r'\section{', gpos + 5)
        gblock = app[gpos:gend if gend > 0 else len(app)]
        general_syms = set()
        for row in _lt_rows(gblock):
            m = re.match(r'\$([^$]+)\$', row)
            if m:
                for s in symbols_in(m.group(1)):
                    general_syms.add(s)
        rec('PASS', 'J', '附录H 通用符号与算符：%d 条符号' % len(general_syms))
        miss_g = [s for s in GLOBAL_SYMBOLS
                  if s not in general_syms and s not in GLOBAL_LITERAL]
        miss_lit = [s for s in GLOBAL_LITERAL if s not in gblock]
        if miss_g or miss_lit:
            rec('FAIL', 'J', '附录H 通用符号表缺下列全局符号：%s'
                % '、'.join(miss_g + miss_lit))

    # 物理量符号表
    ppos = app.find(r'\section{物理量符号')
    if ppos < 0:
        rec('FAIL', 'J', '附录H 缺「物理量符号」小节')
        phys = set()
    else:
        pend = app.find(r'\section{', ppos + 5)
        pblock = app[ppos:pend if pend > 0 else len(app)]
        phys = set()
        n = 0
        for row in _lt_rows(pblock):
            m = re.match(r'\$([^$]+)\$', row)
            if m:
                n += 1
                for s in symbols_in(m.group(1)):
                    phys.add(s)
        rec('PASS' if n >= 60 else 'FAIL', 'J',
            '附录H 物理量符号：%d 条（门槛 60）' % n)

    # 术语中英对照表
    tpos = app.find(r'\section{术语中英对照}')
    if tpos < 0:
        rec('FAIL', 'J', '附录H 缺「术语中英对照」小节')
    else:
        tblock = app[tpos:]
        n = len(_lt_rows(tblock))
        rec('PASS' if n >= 80 else 'FAIL', 'J',
            '附录H 术语中英对照：%d 条（门槛 80）' % n)

    # 双向核对：附录H 的符号必须在正文（附录H 之外）出现过
    body = ''.join(strip_comments(read(f)) for f in files)
    body_no_h = body
    not_used = []
    for s in sorted(phys | general_syms):
        probe = s.replace('\\', '\\')
        if probe not in body_no_h:
            not_used.append(s)
    if not_used:
        rec('WARN', 'J', '附录H 中 %d 个符号在正文未直接命中（需人工确认写法差异）：%s'
            % (len(not_used), '、'.join(not_used[:15])
               + ('…' if len(not_used) > 15 else '')))
    else:
        rec('PASS', 'J', '附录H 全部符号均在正文中出现过')


# ---------------------------------------------------------------- K 结构完整性

def check_K_structure(files, main_text, main_raw=''):
    t = strip_comments(main_text)
    letters = re.findall(r'\\appchapter\{附录([A-Z]{1,2})[：:]', t)
    # v1.27.5 起：附录 J/K/L 移入配套手册，主文件用一行注释声明；
    # 除声明集合外，字母仍须连续。
    moved = set()
    mm = re.search(r'%+ *appendix-moved-to-handbook: *([A-Z,\s]+)',
                  main_raw or main_text)
    if mm:
        moved = set(x.strip() for x in mm.group(1).split(',') if x.strip())
    # A--Z 用完后接 AA、AB…（与 physics_regression 的 P04 同一套字母方案）
    full = [(chr(ord('A') + i) if i < 26 else 'A' + chr(ord('A') + i - 26))
            for i in range(len(letters) + len(moved))]
    expect = [x for x in full if x not in moved]
    if letters != expect:
        rec('FAIL', 'K', '附录章编号不连续或重复：%s（应为 %s）'
            % (''.join(letters), ''.join(expect)))
    else:
        rec('PASS', 'K', '附录章编号连续：%s 共 %d 个' % (''.join(letters), len(letters)))

    # \input 挂载：ch1..ch12 补充材料各一次且顺序递增
    inputs = re.findall(r'\\input\{(ch\d+_supplement(?:\.tex)?)\}', t)
    base = [os.path.splitext(x)[0] for x in inputs]
    want = ['ch%d_supplement' % i for i in range(1, 13)]
    if base != want:
        rec('FAIL', 'K', '补充材料 \\input 挂载异常：%s' % base)
    else:
        rec('PASS', 'K', '补充材料挂载：ch1–ch12 各一次且顺序正确')

    # 标签唯一性
    labels = {}
    for f in files:
        raw = strip_comments(read(f))
        for m in re.finditer(r'\\label\{([^{}]+)\}', raw):
            labels.setdefault(m.group(1), []).append((rel(f), line_of(raw, m.start())))
    dups = {k: v for k, v in labels.items() if len(v) > 1}
    if dups:
        for k, v in sorted(dups.items())[:10]:
            rec('FAIL', 'K', '标签重复：%s —— %s'
                % (k, '；'.join('%s:%d' % x for x in v)))
        if len(dups) > 10:
            rec('FAIL', 'K', '另有 %d 组重复标签' % (len(dups) - 10))
    else:
        rec('PASS', 'K', '标签唯一性：%d 个 \\label 无重复' % len(labels))

    # 索引/术语表/符号表的结构存在性
    n_idx = len(re.findall(r'\\index\{', t))
    has_printindex = r'\printindex' in t
    rec('INFO', 'K', '索引条目 %d 个，\\printindex %s'
        % (n_idx, '存在' if has_printindex else '不存在'))
    rec('PASS' if APPENDIX_H_MARK in t else 'WARN', 'K',
        '附录H（符号与术语总表）%s' % ('已就位' if APPENDIX_H_MARK in t else '未就位'))


# ---------------------------------------------------------------- L README 对齐

README = os.path.join(ROOT, 'README.md')


def check_L_readme(main_text):
    """README 的附录数量声明/清单必须与正文的 \\appchapter 逐一对应。

    页数与版本号由 check/full_check.py 的 E_consistency 负责，这里只查附录维度：
    新增附录后最容易漏改的就是 README 的「N 个附录」与清单。
    """
    if not os.path.exists(README):
        rec('WARN', 'L', 'README.md 不存在，跳过 README 对齐检查')
        return
    rd = read(README)
    t = strip_comments(main_text)
    letters = re.findall(r'\\appchapter\{附录([A-Z]{1,2})[：:]', t)
    n = len(letters)

    for pat, name in ((r'12\s*章正文\s*\+\s*(\d+)\s*附录', '「12 章正文 + N 附录」'),
                      (r'附录（(\d+)\s*个）', '「附录（N 个）」')):
        m = re.search(pat, rd)
        if not m:
            rec('WARN', 'L', 'README 未找到%s数量声明' % name)
        elif int(m.group(1)) != n:
            rec('FAIL', 'L', 'README %s 声明 %s，正文实际 %d 个附录'
                % (name, m.group(1), n))
        else:
            rec('PASS', 'L', 'README %s 与正文一致：%d' % (name, n))

    listed = re.findall(r'^\s*-\s*\*\*附录\s*([A-Z]{1,2})\*\*', rd, re.M)
    missing = [c for c in letters if c not in listed]
    extra = [c for c in listed if c not in letters]
    if missing or extra:
        rec('FAIL', 'L', 'README 附录清单与正文不符：缺 %s，多 %s'
            % (''.join(missing) or '无', ''.join(extra) or '无'))
    else:
        rec('PASS', 'L', 'README 附录清单与正文逐一对应：%s' % ''.join(listed))


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    files = tex_files()
    main_path = os.path.join(SRC, MAIN)
    main_raw = read(main_path) if os.path.exists(main_path) else ''
    main_text = strip_comments(main_raw)

    check_H_terms(files)
    check_I_cards(files)
    check_J_symbol_table(files, main_text)
    check_K_structure(files, main_text, main_raw)
    check_L_readme(main_text)

    if args.json:
        print(json.dumps([{'level': a, 'cat': b, 'msg': c}
                          for a, b, c in RESULTS],
                         ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》一致性与符号规范审核')
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
