# -*- coding: utf-8 -*-
r"""生成 docs/self-study-audit.md（判据来自 check/selfstudy_criteria.py）。

用法：
    python tools\selfstudy_audit.py            # dry-run：只打印缺口
    python tools\selfstudy_audit.py --write    # 覆盖 docs/self-study-audit.md
"""
from __future__ import print_function

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'src')
OUT = os.path.join(ROOT, 'docs', 'self-study-audit.md')
sys.path.insert(0, os.path.join(ROOT, 'check'))

import selfstudy_criteria as SC      # noqa: E402


def main():
    write = '--write' in sys.argv
    ch, tools, topics, gaps = SC.audit(SRC)
    print(u'审计：%d 章、%d 工具型、%d 专题型；缺口 %d 项'
          % (len(ch), len(tools), len(topics), len(gaps)))
    for who, name in gaps[:20]:
        print(u'  缺：%s — %s' % (who, name))
    text = SC.markdown(SRC)
    # 面向读者的口径：A–H 内联在主文件、没有独立源文件；J–L 已移入配套手册，仍按源文件审计
    text = text.replace('tmp/selfstudy_audit.py', 'tools/selfstudy_audit.py')
    text = text.replace('tmp\\selfstudy_audit.py', 'tools\\selfstudy_audit.py')
    text = text.replace('工具型附录（A–L）', '工具型附录（A–I）')
    note = (u'\n\n口径说明：A–H 的正文内联在主文件里，没有独立源文件，本表只审有独立文件的附录；'
            u'J／K／L（经典教材对照、跨教材例题、分章自测题库）已移出主册、编入配套手册，'
            u'本表仍按源文件审计，方便改稿时核对。\n')
    marker = '判据同时接受，不强行统一。'
    if marker in text:
        text = text.replace(marker, marker + note, 1)
    text = text.replace('## 二、工具型附录（A–I）',
                        '## 二、工具型附录（A–I；J–L 已移入配套手册，仍按源文件列出）')
    text = text.replace(
        '**无缺口**：12 章、4 个工具型附录、40 个专题附录全部满足各自标准。',
        '**无缺口**：12 章、40 个专题附录，以及 4 个有独立源文件的工具型附录（I–L）全部满足各自标准'
        '（A–H 内联在主文件中，不在本表）。')
    data = text.encode('utf-8')
    if not write:
        print(u'[DRY-RUN] 未写 %s（将写入 %d 字节）' % (OUT, len(data)))
        return 0
    if not os.path.isdir(os.path.dirname(OUT)):
        os.makedirs(os.path.dirname(OUT))
    io.open(OUT, 'wb').write(data)
    print(u'[OK] 已写 %s（%d 字节）' % (OUT, len(data)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
