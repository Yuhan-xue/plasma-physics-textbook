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
