# -*- coding: utf-8 -*-
"""在 Wikimedia Commons 上按概念检索图像，抓取许可/尺寸元数据并给出可下载直链。

用法：
  python tools\\commons_pick.py            # 全部概念
  python tools\\commons_pick.py tokamak    # 只跑某个键

只使用标准库（urllib），自动绕过本机死代理。
"""
from __future__ import print_function

import json
import os
import re
import sys
import urllib.parse
import urllib.request

API = 'https://commons.wikimedia.org/w/api.php'
UA = 'plasma-textbook-figure-research/1.0 (educational; local use)'
ALLOW = re.compile(r'(cc0|cc-by|cc by|public domain|pd-|gfdl)', re.I)
BAD = re.compile(r'(-nc|-nd|non-?commercial|no-?deriv)', re.I)
OK_LIC = re.compile(r'^(cc0|cc by|cc-by|public domain|pd)', re.I)

CONCEPTS = {
    'tokamak': 'tokamak interior plasma',
    'tokamak_cutaway': 'tokamak schematic cutaway diagram',
    'stellarator': 'Wendelstein 7-X stellarator',
    'iter': 'ITER tokamak construction',
    'nif': 'National Ignition Facility target chamber',
    'icf': 'inertial confinement fusion hohlraum',
    'aurora': 'aurora borealis ISS',
    'solar_corona': 'solar corona prominence',
    'solar_wind': 'solar wind heliosphere',
    'magnetosphere': 'Earth magnetosphere diagram',
    'reconnection': 'magnetic reconnection',
    'glow_discharge': 'glow discharge tube plasma',
    'plasma_globe': 'plasma globe lamp',
    'hall_thruster': 'Hall effect thruster firing',
    'langmuir_probe': 'Langmuir probe plasma',
    'plasma_jet': 'plasma jet atmospheric',
    'zpinch': 'Z-pinch plasma',
    'lightning': 'lightning strike cloud',
    'sun_prominence': 'solar prominence eruption',
    'etch_chamber': 'plasma etching semiconductor chamber',
    'magnetic_mirror': 'magnetic mirror plasma',
    'torus_field': 'toroidal magnetic field lines',
}


def api(params):
    import time
    params = dict(params)
    params['format'] = 'json'
    params['formatversion'] = '2'
    url = API + '?' + urllib.parse.urlencode(params)
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode('utf-8'))
            time.sleep(2.5)   # 礼貌限速：Commons 对匿名脚本会 429
            return data
        except Exception as e:
            last = e
            time.sleep(4.0 * (attempt + 1))
    raise last


def search(term, limit=12):
    d = api({'action': 'query', 'list': 'search', 'srsearch': term,
             'srnamespace': 6, 'srlimit': limit})
    return [x['title'] for x in d.get('query', {}).get('search', [])]


def info(titles):
    out = []
    for i in range(0, len(titles), 10):
        chunk = titles[i:i + 10]
        d = api({'action': 'query', 'titles': '|'.join(chunk), 'prop': 'imageinfo',
                 'iiprop': 'url|size|extmetadata|mime'})
        for pg in d.get('query', {}).get('pages', []):
            ii = (pg.get('imageinfo') or [{}])[0]
            if not ii:
                continue
            md = ii.get('extmetadata', {})
            g = lambda k: (md.get(k, {}) or {}).get('value', '')
            lic = re.sub(r'<[^>]+>', '', g('LicenseShortName')).strip()
            artist = re.sub(r'<[^>]+>', '', g('Artist')).strip()
            artist = re.sub(r'\s+', ' ', artist)[:60]
            out.append({
                'title': pg.get('title', ''),
                'w': ii.get('width', 0), 'h': ii.get('height', 0),
                'mime': ii.get('mime', ''), 'url': ii.get('url', ''),
                'license': lic, 'artist': artist,
                'desc': re.sub(r'<[^>]+>', '', g('ImageDescription')).strip()[:110],
            })
    return out


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
        print('用法：python -X utf8 tools/commons_pick.py [概念键 …]')
        print('候选概念：' + '、'.join(CONCEPTS))
        return
    keys = sys.argv[1:] or list(CONCEPTS)
    for k in keys:
        term = CONCEPTS[k]
        print('=' * 100)
        print('## %s   (search: %s)' % (k, term))
        try:
            titles = search(term)
        except Exception as e:
            print('  search failed: %s' % e)
            continue
        try:
            rows = info(titles)
        except Exception as e:
            print('  info failed: %s' % e)
            continue
        rows = [r for r in rows if r['mime'].startswith('image/')]
        good = [r for r in rows if r['w'] >= 1200 and OK_LIC.match(r['license'] or '')
                and not BAD.search(r['license'] or '')]
        good.sort(key=lambda r: -(r['w'] * r['h']))
        if not good:
            good = [r for r in rows if r['w'] >= 800 and OK_LIC.match(r['license'] or '')
                    and not BAD.search(r['license'] or '')]
        for r in good[:3]:
            fn = r['title'].replace('File:', '').replace(' ', '_')
            direct = 'https://commons.wikimedia.org/wiki/Special:FilePath/' + \
                urllib.parse.quote(fn) + '?width=2000'
            print('  %-58s %5dx%-5d %-18s %s' % (r['title'][5:63], r['w'], r['h'], r['license'],
                                                 r['artist'][:34]))
            print('      direct: %s' % direct)
            print('      page  : https://commons.wikimedia.org/wiki/%s' %
                  urllib.parse.quote(r['title'].replace(' ', '_')))
        if not good:
            print('  (无合格候选，前 3 个含许可信息:)')
            for r in rows[:3]:
                print('     %-52s %5dx%-5d %-16s' % (r['title'][5:57], r['w'], r['h'], r['license']))
        print()


if __name__ == '__main__':
    main()
