# -*- coding: utf-8 -*-
"""可证伪性夹具：用合成 PDF 检验 check/figure_layout.py 真的能抓到图形缺陷。

五张合成页（PyMuPDF 现场生成，不依赖仓库里的历史 PDF）：
  good     两张图各有 4 个互不交叠的节点框、标签都在框内        -> FAIL 0
  overlap  两个节点框重叠 1000 pt² 以上                        -> FAIL（节点框交叠）
  cross    一块文字横跨两个节点框 + 一个字被框边裁掉            -> FAIL（跨框 / 裁切）
  gone     只有题注关键字、没有节点框（图形被删）                -> FAIL（图形消失）
  textcut  示意图区间里两块文字叠在一起 + 红色线条穿过文字        -> FAIL（规则 6/7）
  edge     节点框右沿探出 v1.9.0 版心（1.75/1.95 cm）           -> FAIL（越出版心）

版式约定（与 check/figure_layout.py 的区间算法对应）：
  图题注 y 越小，其纵向区间为 [10, 本题注+40]；后一张图的区间从「前一张题注+10」开始，
  所以后一张图的节点框必须整体位于前一张题注下方 40 pt 以外，否则会被前一张图吞掉。
  示意图（图 9.3）只需题注关键字落在区间内即可被识别，故三张题注都要写进合成页，
  否则套件会按「在整本 PDF 中都没找到」判 FAIL。

输出全为 ASCII（控制台是 GBK），套件完整输出写入 tmp/figlayout-fixture/logs/*.txt。

用法：python check/fixtures/gate_figure_layout.py
依赖：PyMuPDF（pip install PyMuPDF）
"""
import hashlib
import io
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUITE = os.path.join(ROOT, 'check', 'figure_layout.py')
FIXDIR = os.path.join(ROOT, 'tmp', 'figlayout-fixture')

CAP1 = '图 9.1 漂移波（可能）不稳定性的正反馈回路。'
CAP2 = '图 9.2 等离子体不稳定性分类树状图。'
CAP3 = '图 9.3 瑞利–泰勒不稳定性物理图像。尖峰处重流体下沉，谷处轻流体上升。'

# 图 9.1：4 个框在 y=100..250；题注 y=300
BOX1 = [(60, 100, 200, 140), (300, 100, 440, 140), (60, 180, 200, 220), (300, 180, 440, 220)]
TXT1 = ['密度涨落', '电势涨落', '密度涨落被放大', 'E×B 漂移']
CAP1_Y = 300
# 图 9.2：4 个框在 y=420..550（在前一题注+40 之外）；题注 y=700
BOX2 = [(60, 420, 200, 460), (300, 420, 440, 460), (60, 510, 200, 550), (300, 510, 440, 550)]
TXT2 = ['等离子体', '宏观 MHD', '微观动理学', '漂移波']
CAP2_Y = 700
CAP3_Y = 780


def build_pdf(path, kind):
    import pymupdf
    doc = pymupdf.open()
    page = doc.new_page()                                  # A4：595x842 pt
    page.draw_rect(pymupdf.Rect(20, 20, 575, 800), color=None, fill=(0.98, 0.98, 0.98))
    b1, b2 = list(BOX1), list(BOX2)
    if kind == 'overlap':
        b1[1] = (150, 110, 290, 150)                       # 与 b1[0] 重叠 >1000 pt²
    elif kind == 'cross':
        b1[1] = (205, 100, 345, 140)                       # 与 b1[0] 紧邻，留 5 pt 缝
    elif kind == 'gone':
        b1, b2 = [], []
    elif kind == 'edge':
        b1[0] = (430, 100, 575, 140)                       # 右沿 575 pt > 版心右界 540 pt
        b1[1] = (150, 100, 290, 140)                       # 挪开，避免同时触发交叠判定
    for r in b1 + b2:
        page.draw_rect(pymupdf.Rect(*r), color=(0.2, 0.2, 0.2), fill=(0.9, 0.95, 1.0), width=0.8)
    for (x0, y0, x1, y1), txt in zip(b1, TXT1 * 2):
        page.insert_text((x0 + 12, y0 + 26), txt, fontname='china-s', fontsize=11)
    for (x0, y0, x1, y1), txt in zip(b2, TXT2 * 2):
        page.insert_text((x0 + 12, y0 + 26), txt, fontname='china-s', fontsize=11)
    if kind == 'cross':
        # 一块文字正好骑在两框接缝上（73% 落在左框、23% 落在右框 → 跨框）
        page.insert_text((95, 126), '把高密度区输运到低密度区', fontname='china-s', fontsize=11)
        # 一个字被框下沿裁掉（框内 45%、探出 6 pt → 裁切）
        page.insert_text((60, 146), '越界文字', fontname='china-s', fontsize=11)
    if kind == 'textcut':
        # 示意图区间（[710, 820]）内：两块文字叠在一起（规则 6）
        page.insert_text((100, 745), '重流体区域标签', fontname='china-s', fontsize=11)
        page.insert_text((112, 745), '重力箭头标签', fontname='china-s', fontsize=11)
        # 红色线条穿过上面的文字（规则 7）
        page.draw_line(pymupdf.Point(60, 739), pymupdf.Point(400, 739),
                       color=(1, 0, 0), width=1.5)
    page.insert_text((150, CAP1_Y), CAP1, fontname='china-s', fontsize=11)
    page.insert_text((150, CAP2_Y), CAP2, fontname='china-s', fontsize=11)
    page.insert_text((150, CAP3_Y), CAP3, fontname='china-s', fontsize=11)
    doc.save(path)
    doc.close()


def run(pdf, logpath):
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    p = subprocess.run([sys.executable, SUITE, '--pdf', pdf],
                       capture_output=True, env=env, cwd=ROOT)
    out = p.stdout.decode('utf-8', 'replace')
    io.open(logpath, 'w', encoding='utf-8', newline='').write(out)
    return out, p.returncode


# ---- 配色检查用例：静态扫源码，用「植入 + 复原」而不是合成 PDF ----

PAL_TARGET = os.path.join(ROOT, 'src', 'ch3_supplement.tex')
PAL_PAT = re.compile(r'(?<![A-Za-z!])fig(?:Blue|Verm|Green|Orange|Purple|Sky|Gray|Ink)')


def implant_palette():
    """把源码里第一处色板名改成内置色名（blue），返回 (path, 原始字节, sha256)。"""
    raw = io.open(PAL_TARGET, 'rb').read()
    text = raw.decode('utf-8')
    crlf = '\r\n' in text
    if crlf:
        text = text.replace('\r\n', '\n')
    m = PAL_PAT.search(text)
    if not m:
        raise RuntimeError('源码里找不到色板名 figXxx，无法植入')
    out = text[:m.start()] + 'blue' + text[m.end():]
    io.open(PAL_TARGET, 'w', encoding='utf-8', newline='').write(
        out.replace('\n', '\r\n') if crlf else out)
    return PAL_TARGET, raw, hashlib.sha256(raw).hexdigest()


def restore_palette(backup):
    path, raw, want = backup
    io.open(path, 'wb').write(raw)
    got = hashlib.sha256(io.open(path, 'rb').read()).hexdigest()
    return got == want


def main():
    try:
        import pymupdf  # noqa: F401
    except ImportError:
        print('PyMuPDF not installed; fixture skipped (pip install PyMuPDF)')
        return 2
    logs = os.path.join(FIXDIR, 'logs')
    if not os.path.isdir(logs):
        os.makedirs(logs)
    cases = [('good', 0, []),
             ('overlap', 1, ['节点框交叠']),
             ('cross', 1, ['文字被节点框边裁切', '文字同时压住两个节点框']),
             ('gone', 1, ['在整本 PDF 中都没找到']),
             ('textcut', 1, ['文字互相压住', '彩色线条/箭头压住文字']),
             ('edge', 1, ['越出版心']),
             # v1.26.3 新增：插图配色规范（静态扫源码）——把色板名改成内置色名即应报警
             ('palette', 1, ['用了写死颜色'])]
    ok = True
    try:
        for kind, want_rc, needles in cases:
            pdf = os.path.join(FIXDIR, kind + '.pdf')
            build_pdf(pdf, kind)
            backup, restored = None, True
            if kind == 'palette':
                backup = implant_palette()
            try:
                out, rc = run(pdf, os.path.join(logs, kind + '.txt'))
            finally:
                if backup is not None:
                    restored = restore_palette(backup)
            if not restored:
                print('[X ] palette 源码未按原字节恢复')
                ok = False
            hits = [n for n in needles if n in out]
            good = (rc == want_rc) and len(hits) == len(needles)
            ok = ok and good
            print('[%s] %-8s exit=%d want=%d needles=%d/%d log=tmp/figlayout-fixture/logs/%s.txt'
                  % ('OK' if good else 'X ', kind, rc, want_rc, len(hits), len(needles), kind))
    finally:
        for f in os.listdir(FIXDIR):
            if f.endswith('.pdf'):
                os.remove(os.path.join(FIXDIR, f))
    print('fixture verdict: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
