# -*- coding: utf-8 -*-
"""插图来源的单一分类口径（docs/figures.md 与 docs/overview.png 共用）。

只分两类，规则写死在这里，避免两处文档各写一套数字：
  inst   = 政府机构、国家实验室、大学与研究机构（含 ITER/ESA 等国际组织）
  person = 个人摄影师、企业与开放社区

用法：
    from figure_sources import load_rows, counts
    rows = load_rows('figures/CREDITS.md')
    counts(rows)   # 例：{'inst': 45, 'person': 73, 'total': 118}（随 CREDITS.md 变）
"""
import io
import os
import re

ROW = re.compile(r'^\| `([^`]+)` \| (.*?) \| (.*?) \| (.*?) \| (\S+) \| (.*?) \|\s*$')

INST_KW = (
    # 政府 / 国家实验室 / 国际组织
    'NASA', 'ESA', 'European Space Agency', 'ITER', 'LLNL', 'Lawrence Livermore',
    'Lawrence Berkeley', 'Los Alamos', 'INL', 'Idaho National', 'EUROfusion',
    'USAF', 'U.S.', 'Air Force', 'Navy', 'Department of Energy', 'DOE', 'CERN',
    'PPPL', 'Princeton Plasma', 'NIF', 'ORNL', 'Oak Ridge', 'Sandia', 'JET',
    'UKAEA', 'CCFE', 'NOAA', 'Goddard', 'SDO', 'STEREO', 'JAEA', 'QST', 'NIFS',
    'CEA', 'General Atomics', 'IPP', 'Max Planck', 'Forschungszentrum', 'KIT',
    'ENEA', 'CNRS', 'DIII-D', 'ASDEX', 'KSTAR', 'EAST', 'Wendelstein',
    # 大学 / 研究所 / 天文台
    'Universit', 'Institut', 'Politecnico', 'College', 'School', 'Academy',
    'Técnico', 'Observatory', 'Hospital', 'Laboratory', 'Lab',
)


def load_rows(credits_path):
    rows = []
    with io.open(credits_path, encoding='utf-8') as f:
        for line in f:
            m = ROW.match(line.rstrip('\n'))
            if not m:
                continue
            fname, name, author, lic, url, where = m.groups()
            rows.append({'file': fname, 'name': name, 'author': author.strip(),
                         'lic': lic, 'url': url, 'where': where.strip()})
    return rows


def classify(author):
    a = author or ''
    return 'inst' if any(k in a for k in INST_KW) else 'person'


def counts(rows):
    cls = [classify(r['author']) for r in rows]
    return {'inst': cls.count('inst'), 'person': cls.count('person'),
            'total': len(cls)}


if __name__ == '__main__':
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(counts(load_rows(os.path.join(root, 'figures', 'CREDITS.md'))))
