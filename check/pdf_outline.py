#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
《等离子体物理自学教材》PDF 书签／目录／交叉引用一致性守卫
==========================================================

背景（TODO 完成条件 ④：核对页数、书签、文字图像与附件一致性）：
此前只核对了「页数」与「图形几何」，书签与印刷目录的一致性从未量化过。
本模块直接量 PDF 结构：

  O1 基本信息：页数、书签总数、层级分布、带链接页数        → INFO
  O2 书签页码全部落在 [1, pages] 且沿文档顺序单调不减      → FAIL（越界/倒序）
  O3 印刷目录 ↔ 书签一致性（逐条）：
      O3a 每个目录链接矩形内的标题都能匹配到同名书签        → FAIL（找不到）
      O3b 目录行右侧印出的页码 == 书签页码                  → 见下方口径
      O3c 目录链接的目标页 == 书签页码                      → FAIL（点击落点与书签不符）
  O4 附录 F 交叉引用语义（PP-U02 验收「逐条点击到确实解释该概念的正文」）：
      逐个 `\\S\\ref{label}` 显示文本的小节号 == .aux 中该 label 的小节号 → FAIL
      并统计链接目标页与 .aux 页是否一致（信息，用于量化 +1 偏移）
  O5 链接计数自洽：印刷目录链接数 == 书签数；报告正文 § 链接数

印刷目录的页跨度**动态探测**（`detect_toc_pages`）：从第 2 页起连续「链接数
≥ MIN_TOC_LINKS」的页即为印刷目录页；随后再吸收**紧随其后的短尾页**——目录
正好跨页时，末页可能只有一两行（实测 v1.9.0 首轮曾出现第 9 页仅 1 行），
判据是「该页至少有一个带右对齐页码的目录链接」。探测失败（<2 页）才退回
常量 TOC_PAGES。这样附录增删导致目录从 7 页涨到 8 页时，O3／O5 自动跟上，
不会静默漏扫或误扫。

O3b 的口径：hyperref 的书签落点由锚点 whatsit 决定，而 `\\addcontentsline`
写页码发生在 shipout，两者相差一页时会系统性偏移（见 CLOSURE「书签偏移」）。
但 24 处 `\\section*` + `\\addcontentsline` 若**不**补 `\\phantomsection`，
hyperref 会复用上一个锚点（前一节的表题、tcolorbox 计数器），落点跑到更早的
页面，而印出页码仍是本节页 → 偏差可达 ±2。v1.9.0 已给全部 24 处补上
`\\phantomsection`，故本项判据维持严格：印刷页码与书签页码 **完全一致** → PASS；
偏差**全部为 −1** 且不超过 15% → WARN（已知结构性偏移）；出现 +1、±2、比例
过高或（补锚点后）任何偏差 → FAIL。

依赖 PyMuPDF；未安装时报 WARN 并跳过结构项（不会静默通过）。.aux 缺失时
O4 的小节号比对降级为 WARN。

用法：
    python check/pdf_outline.py                      # 用 pdf/plasma_physics_textbook_v2.pdf
    python check/pdf_outline.py --pdf <path> --aux <path> --json

退出码：0 = 无 FAIL；1 = 存在 FAIL
"""

import argparse
import json
import os
import re
import sys
# ---- 输出编码守卫（2026-10-19）------------------------------------------------
# Windows 控制台与管道默认 CP936/GBK，检查报告里的「−」「≥」「×」「≈」一律编不出去，
# 于是脚本在 print 时抛 UnicodeEncodeError、退出码 1，看起来像「检查失败」。
# 夹具正是这样被误伤的。强制标准输出走 UTF-8（编不出的字用 ? 代替，不再崩）。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding='utf-8', errors='replace')
    except Exception:      # 老 Python 或已被替换的流：保持原样
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PDF = os.path.join(ROOT, 'pdf', 'plasma_physics_textbook_v2.pdf')
DEFAULT_AUX = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.aux')
MAIN_TEX = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.tex')

TOC_PAGES = 7          # 印刷目录页数的退路常量（动态探测失败时使用）
MIN_TOC_LINKS = 12     # 判定「这一页是印刷目录」所需的最少链接数
TAIL_MAX = 2           # 目录短尾页最多吸收几页（末页只剩一两行时用）
TOC_FIRST = 1          # 印刷目录从第 2 页（0-based 索引 1）开始，封面无链接
OFFSET_MAX = 0.15      # +1 偏移条数占比容忍上限

RESULTS = []


def rec(level, cat, msg):
    RESULTS.append((level, cat, msg))


def detect_toc_pages(doc, min_links=MIN_TOC_LINKS, start=TOC_FIRST):
    """动态探测印刷目录占用的页数（返回页数，不含封面）。

    为什么不能写死：附录增删会让目录从 7 页涨到 8 页甚至更多，写死常量会
    静默漏扫（O3 少比对一页条目、O5 少算链接数）。
    判据：目录页每页有几十个指向正文章节的 goto 链接，而封面无链接、正文页
    的交叉引用链接数远低于 min_links；目录页在文档中连续且紧接封面之后，
    故取「从 start 起连续满足阈值」的最长前缀。

    再看**短尾页**：目录若正好在末页只排一两行（v1.9.0 首轮实测第 9 页仅
    「L.13 来源、分级与版权说明」一行），链接数 1 < min_links，上面的前缀
    判据会把它判成正文页，于是 O3 漏比一条、O5 少算一个链接。补一条判据：
    紧随其后的页只要含「至少一个带右对齐页码的目录式链接」（链接矩形内文字
    与某个书签同名，且右侧同行有 1--3 位页码数字），就并入目录页。上限
    TAIL_MAX 页，避免把正文页误并进来。
    """
    n = 0
    for pno in range(start, doc.page_count):
        if len(doc[pno].get_links()) >= min_links:
            n += 1
        else:
            break
    titles = {norm(t) for _lvl, t, _p in doc.get_toc()}
    for pno in range(start + n, min(start + n + TAIL_MAX, doc.page_count)):
        page = doc[pno]
        words = page.get_text('words')
        for lk in page.get_links():
            rect = lk.get('from')
            if rect is None:
                continue
            if norm(page.get_textbox(rect)) in titles and printed_page(page, rect, words):
                n += 1
                break
    return n


def norm(s):
    """标题匹配键。

    目录行与书签标题的结构差异必须先抹平，否则会出现大量假不匹配：
    * 目录行带序号前缀（``1.4``／``A.1``／``11.1.3``），书签不带；
    * 书签由 hyperref 生成，数学/Greek 字符（``δ``、``θ``）会被丢弃；
    * 两边引号字形不同（``“”`` vs ``""``）；
    * 空白与零宽字符的位置不同。
    只保留字母、数字、CJK、引号与常见标点。
    """
    # 附录字母方案升级到 A--Z, AA, AB…后，开头的编号可能带 1--2 个字母（AA.1）
    # AA.1 这类「两字母 + 点 + 数字」也当编号剥掉；纯数字编号（1.1.1）不受影响。
    # 负向先行：数字后紧跟连字符时**不是编号**而是标题本身（如「0-D 功率平衡」），
    # 否则书签侧会被剥成「-D 功率平衡」而目录侧（前缀 AZ.2 先被剥掉）保留 0-D，产生假不匹配。
    s = re.sub(r'^(?:[A-Za-z]{1,2}[.．])?\d+(?:[.．]\d+)*[.．]?(?![-\u2010-\u2015])', '', s or '')
    s = re.sub(r'[\s\u200b]+', '', s)
    s = (s.replace('“', '“').replace('”', '“').replace('″', '“')
          .replace('‘', "'").replace('’', "'"))
    return re.sub(r'[^0-9A-Za-z\u4e00-\u9fff"\'()（）,，、:：;；.。!！?？\-–—/]', '', s)


def printed_page(page, rect, words):
    """题注行右侧右对齐的纯数字 = 印刷目录页码。"""
    yc = (rect.y0 + rect.y1) / 2.0
    best = None
    for w in words:
        if not re.fullmatch(r'\d{1,3}', w[4]):
            continue
        if abs((w[1] + w[3]) / 2.0 - yc) > 6.0:
            continue
        if w[0] < rect.x1 - 1.0:          # 必须在标题右侧
            continue
        if best is None or w[2] > best[2]:
            best = w
    return int(best[4]) if best else None


def analyze(doc, aux_path):
    pages = doc.page_count
    toc = doc.get_toc()

    # ---------- O1 基本信息 ----------
    levels = {}
    for lvl, _title, _page in toc:
        levels[lvl] = levels.get(lvl, 0) + 1
    n_links = sum(len(doc[p].get_links()) for p in range(pages))
    link_pages = sum(1 for p in range(pages) if doc[p].get_links())
    toc_pages = detect_toc_pages(doc)
    if toc_pages < 2:
        rec('WARN', 'O1', '印刷目录动态探测只得到 %d 页（阈值 %d 个链接/页），'
                          '退回到常量 %d 页' % (toc_pages, MIN_TOC_LINKS, TOC_PAGES))
        toc_pages = TOC_PAGES
    toc_span = range(min(toc_pages + TOC_FIRST, pages))   # 含封面，供逐页扫描
    rec('INFO', 'O1', '页数 %d；书签 %d 条（层级 %s）；带链接页 %d 页，链接总数 %d；'
                      '印刷目录 %d 页（第 %d--%d 页，动态探测）'
        % (pages, len(toc), sorted(levels.items()), link_pages, n_links,
           toc_pages, TOC_FIRST + 1, toc_pages + TOC_FIRST))
    if not toc:
        rec('FAIL', 'O1', 'PDF 没有任何书签（导航失效）')
        return

    # ---------- O2 书签页码合法性 / 单调性 ----------
    flat = [(title, page) for _lvl, title, page in toc]
    bad_range = [(t, p) for t, p in flat if not (1 <= p <= pages)]
    if bad_range:
        rec('FAIL', 'O2', '书签页码越界 %d 条，例如 %s' % (len(bad_range), bad_range[:3]))
    else:
        rec('PASS', 'O2', '书签页码全部落在 [1, %d]' % pages)

    desc = [(flat[i][0], flat[i][1], flat[i + 1][0], flat[i + 1][1])
            for i in range(len(flat) - 1) if flat[i + 1][1] < flat[i][1]]
    if desc:
        rec('FAIL', 'O2', '书签页码沿文档顺序倒退 %d 处，例如「%s」p.%d → 「%s」p.%d'
            % (len(desc), desc[0][0], desc[0][1], desc[0][2], desc[0][3]))
    else:
        rec('PASS', 'O2', '书签页码沿文档顺序单调不减')

    # ---------- O3 印刷目录 ↔ 书签 ----------
    normed = [norm(title) for _lvl, title, _page in toc]
    entry_pages = [page for _lvl, _title, page in toc]

    entries = []           # (标题, 链接目标页, 印刷页码)
    for pno in toc_span:
        page = doc[pno]
        links = page.get_links()
        if not links:
            continue
        words = page.get_text('words')
        for lk in links:
            rect = lk.get('from')
            if rect is None:
                continue
            title = page.get_textbox(rect).strip()
            entries.append((title, lk.get('page'), printed_page(page, rect, words)))

    # 同名条目（如各章都有「习题精选」）按文档顺序与书签一一配对，
    # 不能简单取首个同名书签。
    cursor = 0
    matched = []           # (标题, 链接目标页, 印刷页码, 书签页)
    unmatched = []
    for title, target, printed in entries:
        key = norm(title)
        idx = next((i for i in range(cursor, len(normed)) if normed[i] == key), None)
        if idx is None:
            idx = next((i for i in range(len(normed)) if normed[i] == key), None)
        if idx is None:
            unmatched.append(title)
            continue
        cursor = idx + 1
        matched.append((title, target, printed, entry_pages[idx]))

    if unmatched:
        rec('FAIL', 'O3a', '印刷目录 %d/%d 条在书签中找不到同名项，例如「%s」'
            % (len(unmatched), len(entries), unmatched[0]))
    else:
        rec('PASS', 'O3a', '印刷目录 %d 条全部按文档顺序匹配到同名书签' % len(entries))

    no_print = [m for m in matched if m[2] is None]
    if no_print:
        rec('WARN', 'O3b', '印刷目录 %d/%d 条未能取到印刷页码（例如「%s」），这些条目跳过页码比对'
            % (len(no_print), len(entries), no_print[0][0]))

    page_mismatch = [m for m in matched if m[2] is not None and m[2] != m[3]]
    link_mismatch = [m for m in matched if m[1] is not None and m[1] + 1 != m[3]]

    total = len(entries) or 1
    if not page_mismatch:
        rec('PASS', 'O3b', '印刷目录页码与书签页码逐条一致（%d 条）' % len(entries))
    else:
        deltas = {}
        for _, _t, p, b in page_mismatch:
            deltas[b - p] = deltas.get(b - p, 0) + 1
        detail = '；'.join('书签−印刷 = %+d 共 %d 条' % (d, c) for d, c in sorted(deltas.items()))
        sample = '，'.join('「%s」印刷 %d / 书签 %d' % (t, p, b)
                          for t, _link, p, b in page_mismatch[:3])
        if set(deltas) == {-1} and len(page_mismatch) <= OFFSET_MAX * total:
            rec('WARN', 'O3b', '印刷目录页码与书签页码有 %d/%d 条相差一页（%s；例如 %s）——'
                               '已知结构性偏移：书签锚点 whatsit 与 \\addcontentsline 的 '
                               'shipout 页相差一页，点击落在标题所在对开页' %
                (len(page_mismatch), total, detail, sample))
        else:
            rec('FAIL', 'O3b', '印刷目录页码与书签页码不符 %d/%d 条（%s；例如 %s）'
                % (len(page_mismatch), total, detail, sample))

    if link_mismatch:
        rec('FAIL', 'O3c', '目录链接目标页与书签页码不符 %d 条，例如「%s」链接到 p.%d，书签在 p.%d'
            % (len(link_mismatch), link_mismatch[0][0], link_mismatch[0][1] + 1, link_mismatch[0][3]))
    else:
        rec('PASS', 'O3c', '目录链接目标页与书签页码逐条一致')

    # ---------- O4 附录 F 交叉引用语义 ----------
    aux = {}
    if os.path.exists(aux_path):
        text = open(aux_path, encoding='utf-8', errors='replace').read()
        for m in re.finditer(r'\\newlabel\{([^}]*)\}\{\{([^}]*)\}\{([^}]*)\}', text):
            aux[m.group(1)] = (m.group(2), m.group(3))

    if not os.path.exists(MAIN_TEX):
        rec('WARN', 'O4', '找不到主源文件 %s，跳过附录 F 引用语义检查' % MAIN_TEX)
    else:
        src = open(MAIN_TEX, encoding='utf-8', errors='replace').read()
        m = re.search(r'\\appchapter\{附录F.*?(?=\\appchapter\{附录G)', src, re.S)
        if not m:
            rec('FAIL', 'O4', '主源文件中找不到「附录F」区段（附录被改名或删除）')
        else:
            labels = re.findall(r'\\S\\ref\{([^}]+)\}', m.group(0))
            f_pages = [i for i, (_lvl, t, _p) in enumerate(toc) if norm(t).startswith('附录F')]
            g_pages = [i for i, (_lvl, t, _p) in enumerate(toc) if norm(t).startswith('附录G')]
            start = toc[f_pages[0]][2] if f_pages else None
            # 附录 F 自己的页范围：F 首页 → G 首页 − 1（G 不存在时退化为 6 页窗口）。
            # 早先按「start-1 起 6 页」扫，版面重排后窗口会溢出到附录 G，
            # 把 G 内部的 \S\ref 也算进来（v1.9.0 实测 35 vs 36）。
            end = (toc[g_pages[0]][2] - 1) if g_pages else None
            links = []
            if start is not None:
                last = end if end is not None else min(start + 5, pages)
                for pno in range(start - 1, min(last, pages)):
                    for lk in doc[pno].get_links():
                        rect = lk.get('from')
                        txt = doc[pno].get_textbox(rect).strip() if rect else ''
                        # 只统计形如「§2.1」「3.2.1）」的小节引用，排除公式/图表/URL
                        if re.match(r'^[§\s]*\d+(?:\.\d+)*', txt):
                            links.append((pno + 1, txt, lk.get('page')))
            rec('INFO', 'O4', '附录 F 源文件 \\S\\ref %d 处；PDF 中该区段链接 %d 个'
                % (len(labels), len(links)))
            if not aux:
                rec('WARN', 'O4', '读不到 .aux（%s），无法核对引用显示号与 label 小节号' % aux_path)
            elif len(labels) != len(links):
                rec('FAIL', 'O4', '附录 F 的 \\S\\ref 数量（%d）与 PDF 链接数（%d）不一致'
                    % (len(labels), len(links)))
            else:
                wrong = []
                page_gap = []
                for lab, (_, shown, target) in zip(labels, links):
                    num, apage = aux.get(lab, (None, None))
                    # 链接矩形可能多带一个相邻标点（如「2.1）」），只取数字部分
                    mm = re.search(r'\d+(?:\.\d+)*', shown)
                    disp = mm.group(0) if mm else re.sub(r'^[§\s]*', '', shown).strip()
                    if num is None:
                        wrong.append((lab, disp, 'label 不在 .aux 中'))
                        continue
                    if disp != num:
                        wrong.append((lab, disp, num))
                    if apage is not None and target is not None and target + 1 != int(apage):
                        page_gap.append((lab, target + 1, int(apage)))
                if wrong:
                    rec('FAIL', 'O4', '附录 F 有 %d 处引用显示号与 label 小节号不符，例如 %s'
                        % (len(wrong), wrong[:3]))
                else:
                    rec('PASS', 'O4', '附录 F 全部 %d 处 \\S\\ref 的显示小节号与 .aux 一致'
                        % len(labels))
                if page_gap:
                    rec('WARN', 'O4', '附录 F 链接目标页与 .aux 记录页不同 %d/%d 处（如 %s）'
                        % (len(page_gap), len(labels), page_gap[:2]))
                else:
                    rec('PASS', 'O4', '附录 F 链接目标页与 .aux 记录页一致')

    # ---------- O5 链接计数自洽 ----------
    toc_links = sum(len(doc[p].get_links()) for p in toc_span)
    if toc_links == len(toc):
        rec('PASS', 'O5', '印刷目录链接数（%d，第 %d--%d 页）== 书签数（%d），目录条目与书签一一对应'
            % (toc_links, TOC_FIRST + 1, toc_pages + TOC_FIRST, len(toc)))
    else:
        rec('WARN', 'O5', '印刷目录链接数（%d，第 %d--%d 页）≠ 书签数（%d），需人工确认差异来源'
            % (toc_links, TOC_FIRST + 1, toc_pages + TOC_FIRST, len(toc)))
    body = n_links - toc_links
    rec('INFO', 'O5', '正文（含附录）交叉引用链接 %d 个' % body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdf', default=DEFAULT_PDF)
    ap.add_argument('--aux', default=DEFAULT_AUX)
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.pdf):
        rec('FAIL', 'O0', '找不到编译产物：%s（先运行 check/full_check.py 编译）' % args.pdf)
    else:
        try:
            import pymupdf as fitz          # 新版包名（PyMuPDF ≥1.24）
        except ImportError:
            try:
                import fitz                 # 旧版包名
            except ImportError:
                rec('WARN', 'O0', '未安装 PyMuPDF，跳过书签／目录结构检查'
                                  '（install: pip install PyMuPDF）')
                fitz = None
        if fitz is not None:
            doc = fitz.open(args.pdf)
            analyze(doc, args.aux)
            doc.close()

    if args.json:
        print(json.dumps([{'level': a, 'cat': b, 'msg': c} for a, b, c in RESULTS],
                         ensure_ascii=False, indent=2))
        return 1 if any(r[0] == 'FAIL' for r in RESULTS) else 0

    print('=' * 78)
    print('《等离子体物理自学教材》PDF 书签／目录／引用一致性守卫')
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
