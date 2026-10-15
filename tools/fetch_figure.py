# -*- coding: utf-8 -*-
"""下载并规范化一张开放许可配图到 figures/，同时登记到 figures/CREDITS.md。

用法：
  python tools\\fetch_figure.py --url <直链> --name tokamak_jet --title "JET 内部" \
      --author "..." --license "CC BY-SA 4.0" --source "<文件页 URL>" [--width-cm 12] [--credit "说明"]

要点：
  * 运行前会清掉四个代理环境变量，并用 curl --noproxy '*'（本机代理 127.0.0.1:7897 是死的）。
  * 需要 Pillow 计算尺寸；目标 300 dpi（width_cm 对应像素 = width_cm/2.54*300）。
  * 只在成功校验图像后写盘；CREDITS.md 以追加方式维护（LF、末尾换行）。
"""
from __future__ import print_function

import argparse
import io
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, 'figures')
CREDITS = os.path.join(FIGDIR, 'CREDITS.md')

try:
    from PIL import Image
except ImportError:
    Image = None


def clean_env():
    env = dict(os.environ)
    for k in ('HTTP_PROXY', 'HTTPS_PROXY', 'http_proxy', 'https_proxy'):
        env.pop(k, None)
    return env


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--url', required=True)
    ap.add_argument('--name', required=True, help='输出主名（不含扩展名），只用 [a-z0-9_]')
    ap.add_argument('--title', required=True)
    ap.add_argument('--author', required=True)
    ap.add_argument('--license', required=True)
    ap.add_argument('--source', required=True, help='文件页 URL（署名用）')
    ap.add_argument('--width-cm', type=float, default=12.0)
    ap.add_argument('--credit', default='')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    if not os.path.isdir(FIGDIR):
        os.makedirs(FIGDIR)
    raw = os.path.join(FIGDIR, '_raw_' + a.name)
    env = clean_env()
    cmd = ['curl.exe', '-sS', '-L', '--noproxy', '*', '-A',
           'plasma-textbook-figure-fetch/1.0 (educational; contact: repo owner)', '-o', raw, a.url]
    print('curl:', ' '.join(cmd))
    if not a.dry_run:
        r = subprocess.call(cmd, env=env)
        if r != 0:
            print('[FAIL] curl exit=%d' % r)
            return 1
    if a.dry_run:
        print('[dry-run] 不下载')
        return 0

    size = os.path.getsize(raw)
    if size < 5000:
        print('[FAIL] 文件过小 %d B，可能不是图片' % size)
        return 1
    if Image is None:
        print('[WARN] 无 Pillow，跳过尺寸校验')
        w = h = 0
    else:
        try:
            im = Image.open(raw)
            w, h = im.size
        except Exception as e:
            print('[FAIL] 不是可识别图像：%s' % e)
            return 1
    need = int(round(a.width_cm / 2.54 * 300))
    dpi_eff = w / (a.width_cm / 2.54) if a.width_cm else 0
    print('原图 %dx%d px（%d B）；按 %.2f cm 排版 -> 有效 %.0f dpi；目标 %d px' %
          (w, h, size, a.width_cm, dpi_eff, need))

    # 规范化：统一输出 PNG（无损、LaTeX 直接可用）；过小则拒绝
    if Image is not None:
        out = os.path.join(FIGDIR, a.name + '.png')
        im2 = Image.open(raw)
        if im2.mode not in ('RGB', 'L'):
            im2 = im2.convert('RGB')
        im2.save(out, 'PNG', optimize=True)
        print('已写 %s（%d B）' % (out, os.path.getsize(out)))
    else:
        out = raw
    os.remove(raw)

    head = ''
    if not os.path.exists(CREDITS):
        head = ('# 配图来源与许可（figures/）\n\n'
                '本目录图片均为开放许可（CC0 / 公有领域 / CC BY / CC BY-SA），'
                '按各自许可要求随文署名。若某图许可要求「相同方式共享」，'
                '其衍生排版同样以该许可发布。\n\n'
                '| 文件 | 内容 | 作者 | 许可 | 来源页 | 说明 |\n'
                '|---|---|---|---|---|---|\n')
    with io.open(CREDITS, 'a', encoding='utf-8', newline='') as f:
        if head:
            f.write(head)
        f.write('| `%s` | %s | %s | %s | %s | %s |\n' %
                (os.path.basename(out), a.title, a.author, a.license, a.source, a.credit))
    print('已登记到', CREDITS)
    print('[OK] %s  %dx%d px' % (out, w, h))
    return 0


if __name__ == '__main__':
    sys.exit(main())
