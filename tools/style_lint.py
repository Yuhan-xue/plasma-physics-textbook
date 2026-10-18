# -*- coding: utf-8 -*-
"""文风体检：量「套话密度」与「句子长度」，给出改前改后的可比数字。

为什么要有这个脚本：去 AI 味是主观判断，但"套话几处、破折号多密、超长句几句"
是客观的。写作前跑一次 `--rev HEAD` 存底，改完再跑一次，就知道这轮到底改掉了
多少，而不是"感觉干净了"。

口径（与 docs/style.md 一致）：
  - 套话：docs/style.md §2 禁用词表里的空转连接词与自我表演句；
  - 破折号密度：每千个正文字符里出现几个「——」（中文破折号滥用是典型的生成痕迹）；
  - 超长句：去掉数学与命令后，单句 > 100 字；
  - 句长 p90：中文技术教材的经验良好线是 ≤ 60 字。

用法：
    python -X utf8 tools/style_lint.py              # 扫当前工作区 src/*.tex
    python -X utf8 tools/style_lint.py --rev HEAD   # 扫 git 里的某一版（改前存底）
    python -X utf8 tools/style_lint.py --file src/ch3_supplement.tex --list
    python -X utf8 tools/style_lint.py --json

退出码：0 = 只报告；1 = 句长 p90 超 60 或超长句占比 > 3%（可当门禁用）。
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
MAIN = 'plasma_physics_textbook_v2.tex'

# docs/style.md §2 的禁用词表（只收「删掉不影响意思」的那批）
CLICHES = [
    '值得注意的是', '需要指出的是', '不难发现', '不难看出', '众所周知',
    '毫无疑问', '显而易见', '综上所述', '总而言之', '总的来说', '总之',
    '换言之', '简单来说', '从某种意义上说', '在某种意义上', '事实上',
    '实际上', '与此同时', '进一步地', '更重要的是', '由此可见',
    '我们可以看到', '让我们', '你会发现', '本节将', '本书将',
    '起到重要作用', '具有重要意义', '奠定了基础', '开启新篇章',
    '赋能', '助力', '抓手', '闭环', '范式',
]

# 三句同构的排比痕迹（先...其次...最后 的机械枚举）
PATTERNS = [
    ('三段枚举', re.compile(r'首先[^。]{0,60}其次[^。]{0,60}(最后|再次)')),
    ('不是X而是Y', re.compile(r'不是[^。；，]{2,20}而是')),
    ('并非X而是Y', re.compile(r'并非[^。；，]{2,20}而是')),
    ('既X又Y', re.compile(r'既[^。；]{1,18}[也又][^。；]{1,18}')),
]

EMDASH = '——'


def tex_files():
    names = sorted(n for n in os.listdir(SRC) if n.endswith('.tex'))
    return [os.path.join(SRC, n) for n in names]


def read_text(path, rev=None):
    if rev:
        rel = os.path.relpath(path, ROOT).replace('\\', '/')
        try:
            out = subprocess.run(['git', 'show', '%s:%s' % (rev, rel)],
                                 cwd=ROOT, capture_output=True)
        except OSError:
            return None
        if out.returncode != 0:
            return None
        return out.stdout.decode('utf-8', 'replace')
    with io.open(path, encoding='utf-8', newline='') as fh:
        return fh.read()


def clean(line):
    """去掉注释、行内数学与 LaTeX 命令，只留可读中文/英文。"""
    s = re.sub(r'(?<!\\)%.*$', '', line)
    s = re.sub(r'\$[^$]*\$', ' ', s)
    s = re.sub(r'\\[a-zA-Z@]+\*?(\[[^\]]*\])?(\{[^{}]*\})?', ' ', s)
    return s.strip()


def prose_of(raw):
    """→ (正文字符数, 句子列表)。跳过导言区、绘图区、表格与纯命令行。

    表格与 `\\input` 清单不是散文：不剔掉的话，一张二十行的表会被当成
    一个几百字的"句子"，p90 直接失真。
    """
    if not raw:
        return 0, []
    body = raw.split('\\begin{document}', 1)[-1]
    body = re.sub(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}', ' ',
                  body, flags=re.S)
    body = re.sub(r'\\begin\{(tabular|longtable|tabularx)\*?\}.*?'
                  r'\\end\{\1\*?\}', ' ', body, flags=re.S)
    body = re.sub(r'\\begin\{(equation|align|gather|eqnarray)\*?\}.*?'
                  r'\\end\{\1\*?\}', ' ', body, flags=re.S)
    body = re.sub(r'\\begin\{(verbatim|lstlisting|Verbatim)\}.*?'
                  r'\\end\{\1\}', ' ', body, flags=re.S)
    keep = []
    for line in body.split('\n'):
        s = line.strip()
        if re.match(r'^\\(input|include|chapter|section|subsection'
                    r'|subsubsection|label|index|appchapter|addcontentsline'
                    r'|phantomsection|newpage|clearpage|tableofcontents)\b', s):
            continue
        if '&' in s or s.endswith('\\\\'):
            continue
        keep.append(line)
    text = '\n'.join(clean(l) for l in keep)
    chars = len(re.sub(r'\s', '', text))
    sentences = [s.strip() for s in re.split(r'[。！？]', text) if s.strip()]
    return chars, sentences


def pct(vals, q):
    if not vals:
        return 0
    xs = sorted(vals)
    i = min(len(xs) - 1, int(round(q * (len(xs) - 1))))
    return xs[i]


def scan(paths, rev=None, list_hits=False):
    total_chars = 0
    all_lens = []
    dash = 0
    long_sents = 0
    per_file = []
    hits = []
    cliche_counts = {k: 0 for k in CLICHES}
    pattern_counts = {k: 0 for k, _ in PATTERNS}

    for path in paths:
        raw = read_text(path, rev)
        if raw is None:
            continue
        chars, sents = prose_of(raw)
        total_chars += chars
        all_lens.extend(len(s) for s in sents)
        text = '\n'.join(sents)
        dash += text.count(EMDASH)
        n_long = sum(1 for s in sents if len(s) > 100)
        long_sents += n_long
        local = 0
        for k in CLICHES:
            c = text.count(k)
            if c:
                cliche_counts[k] += c
                local += c
                if list_hits:
                    for i, l in enumerate(raw.split('\n'), 1):
                        if k in clean(l):
                            hits.append((os.path.basename(path), i, k,
                                         clean(l)[:70]))
        for name, rx in PATTERNS:
            ms = rx.findall(text)
            if ms:
                pattern_counts[name] += len(ms)
                local += len(ms)
        per_file.append({
            'file': os.path.basename(path),
            'chars': chars,
            'p90': pct([len(s) for s in sents], 0.90),
            'long': n_long,
            'cliche': local,
        })

    p90 = pct(all_lens, 0.90)
    return {
        'chars': total_chars,
        'sentences': len(all_lens),
        'p90': p90,
        'p50': pct(all_lens, 0.50),
        'max': max(all_lens) if all_lens else 0,
        'long': long_sents,
        'long_ratio': (100.0 * long_sents / len(all_lens)) if all_lens else 0,
        'dash': dash,
        'dash_per_k': (1000.0 * dash / total_chars) if total_chars else 0,
        'cliche': sum(cliche_counts.values()),
        'cliche_per_k': (1000.0 * sum(cliche_counts.values()) / total_chars)
                        if total_chars else 0,
        'cliche_detail': {k: v for k, v in cliche_counts.items() if v},
        'pattern_detail': {k: v for k, v in pattern_counts.items() if v},
        'per_file': sorted(per_file, key=lambda d: -d['cliche']),
        'hits': hits,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rev', help='从 git 读某一版，例如 HEAD（改前存底）')
    ap.add_argument('--file', help='只扫一个文件（相对路径）')
    ap.add_argument('--list', action='store_true', help='逐条列出套话')
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    paths = [os.path.join(ROOT, args.file)] if args.file else tex_files()
    r = scan(paths, args.rev, args.list)

    if args.json:
        print(json.dumps(r, ensure_ascii=False, indent=2))
        return 0

    tag = ('git %s' % args.rev) if args.rev else '工作区'
    print('== 文风体检（%s）==' % tag)
    print('正文可读字符 %d，句 %d' % (r['chars'], r['sentences']))
    print('句长：中位 %d，p90 %d，最长 %d' % (r['p50'], r['p90'], r['max']))
    print('超长句(>100 字) %d 句（%.1f%%）' % (r['long'], r['long_ratio']))
    print('中文破折号 —— %d 处（每千字 %.2f）' % (r['dash'], r['dash_per_k']))
    print('套话 %d 处（每千字 %.2f）' % (r['cliche'], r['cliche_per_k']))
    if r['cliche_detail']:
        print('  ' + '，'.join('%s×%d' % kv for kv in
                              sorted(r['cliche_detail'].items(),
                                     key=lambda kv: -kv[1])))
    if r['pattern_detail']:
        print('  ' + '，'.join('%s×%d' % kv for kv in
                              sorted(r['pattern_detail'].items())))
    print('\n-- 套话最密的 12 个文件 --')
    for d in r['per_file'][:12]:
        if d['cliche']:
            print('  %-34s 套话 %2d  p90 %3d  超长 %2d  %6d 字'
                  % (d['file'], d['cliche'], d['p90'], d['long'], d['chars']))
    if args.list and r['hits']:
        print('\n-- 逐条 --')
        for f, ln, k, t in r['hits']:
            print('  %-30s L%-5d %-8s %s' % (f, ln, k, t))

    bad = r['p90'] > 60 or r['long_ratio'] > 3.0
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
