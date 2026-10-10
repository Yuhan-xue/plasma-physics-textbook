# -*- coding: utf-8 -*-
"""PP-W01 验收夹具：用合成日志驱动 check/full_check.py 与 check/ui_readability.py。

验收口径（reviews/2026-10-10-independent/TODO.md PP-W01）：
  1) 只含 U+000A 伪缺字的日志      -> PASS（可豁免）
  2) 真实缺字（α / Ω 等）           -> 被拒（FAIL），并点名码点与字体
  3) 一次 200 pt 溢出               -> 被拒（FAIL）
  4) 5--50 pt 溢出                  -> WARN（人工复核），不 FAIL
  5) 干净日志                       -> 全 PASS
  6) 真实 TeX 错误（! 开头）        -> 仍可检出（FAIL）

做法：备份真实 .log（编译产物，不是源文件）→ 逐个写入夹具 → 跑两套检查
（`--no-compile` / `--json`）→ 按级别断言消息与退出码 → `finally` 恢复真实日志。

注意：本夹具断言的是**退出码**，而 `full_check.py` 的退出码是整份报告的结论，
所以它隐含要求「其它判定（README 页数、版本号……）处于 FAIL 0」。
若别处还有 FAIL，这里会报「退出码与期望不符」——那说明套件整体尚未修好，不是判定失效；
各条目级别的断言（消息文本与级别）始终独立有效。

依赖：无第三方库；需要已编译过一次（src/plasma_physics_textbook_v2.log 存在）。
用法：python check/fixtures/gate_compile_log.py
"""
import io
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LOG = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.log')

CLEAN = ('This is XeTeX, Version 3.141592653\n'
         'Output written on plasma_physics_textbook_v2.pdf (222 pages).\n')

FIXTURES = [
    ('干净日志', CLEAN,
     [('LaTeX 错误（! 开头）：0', 'PASS'), ('真实缺字 0 处', 'PASS'), ('均 <5 pt', 'PASS')],
     [('缺字：0 处', 'PASS'), ('均 <5 pt', 'PASS')]),
    ('仅 U+000A 伪缺字', CLEAN +
     'Missing character: There is no \n (U+000A) in font '
     '[lmroman12-regular]:mapping=tex-text;!\n',
     [('豁免 U+000A 换行符 1 处，真实缺字 0 处', 'PASS')],
     [('仅 1 处 U+000A 换行符（可豁免）', 'PASS')]),
    ('真实缺字 α(U+03B1)', CLEAN +
     'Missing character: There is no α (U+03B1) in font [lmroman12-regular]!\n',
     [('真实缺字 1 处', 'FAIL'), ('U+03B1', 'FAIL')],
     [('缺字（非 U+000A）：1 处', 'FAIL')]),
    ('真实缺字 Ω(U+03A9)', CLEAN +
     'Missing character: There is no Ω (U+03A9) in font [lmmono10-regular]!\n',
     [('真实缺字 1 处', 'FAIL'), ('U+03A9', 'FAIL')],
     [('缺字（非 U+000A）：1 处', 'FAIL')]),
    ('200 pt 溢出', CLEAN +
     'Overfull \\hbox (200.0pt too wide) in paragraph at lines 100--102\n',
     [('最大 200.0 pt（≥50 pt 属明显排版缺陷）', 'FAIL')],
     [('≥50 pt 属明显排版缺陷', 'FAIL')]),
    ('8 pt 溢出', CLEAN +
     'Overfull \\hbox (8.0pt too wide) in paragraph at lines 200--202\n',
     [('最大 8.0 pt（5--50 pt 需人工复核）', 'WARN')],
     [('5--50 pt 需人工复核', 'WARN')]),
    ('真实 TeX 错误', CLEAN +
     '! Undefined control sequence.\nl.10 \\nosuchmacro\n',
     [('LaTeX 错误（! 开头）：1', 'FAIL')],
     []),
]


def run(script, *extra):
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    r = subprocess.run([sys.executable, os.path.join('check', script)] + list(extra),
                       cwd=ROOT, capture_output=True, env=env)
    return r.stdout.decode('utf-8', errors='replace'), r.returncode


def buckets_text(out):
    """按 '--- LEVEL (n) ---' 分组，返回 {level: [message lines]}"""
    res, cur = {}, None
    for ln in out.split('\n'):
        m = re.match(r'--- (FAIL|WARN|INFO|PASS) \(\d+\) ---', ln.strip())
        if m:
            cur = m.group(1)
            res.setdefault(cur, [])
            continue
        if cur and ln.strip().startswith('['):
            res[cur].append(ln.strip())
    return res


def buckets_json(out):
    res = {}
    for row in json.loads(out):
        res.setdefault(row['level'], []).append(row['msg'])
    return res


def main():
    if not os.path.exists(LOG):
        print('缺少编译日志 %s：先运行 python check/full_check.py' % LOG)
        return 2
    bak = LOG + '.fixture-bak'
    shutil.copy2(LOG, bak)
    ok = True
    try:
        for name, text, exp_full, exp_ui in FIXTURES:
            with io.open(LOG, 'w', encoding='utf-8', newline='') as f:
                f.write(text)
            out_f, rc_f = run('full_check.py', '--no-compile')
            out_u, rc_u = run('ui_readability.py', '--json')
            for out, rc, exp, tag, is_json in (
                    (out_f, rc_f, exp_full, 'full_check', False),
                    (out_u, rc_u, exp_ui, 'ui_readability', True)):
                buck = buckets_json(out) if is_json else buckets_text(out)
                want_fail = any(l == 'FAIL' for _, l in exp)
                if (rc == 1) != want_fail:
                    ok = False
                    print('  [%s] %-16s 退出码 %d 与期望%s不符'
                          % (tag, name, rc, '（应有 FAIL）' if want_fail else '（不应有 FAIL）'))
                for needle, lvl in exp:
                    if not any(needle in m for m in buck.get(lvl, [])):
                        ok = False
                        print('  [%s] %-16s %s 级别未命中：%s'
                              % (tag, name, lvl, needle))
            print('  %-16s 已检' % name)
    finally:
        shutil.copy2(bak, LOG)
        os.remove(bak)
        with io.open(LOG, encoding='utf-8') as f:
            n = f.read().count('Missing character')
        print('已恢复真实日志（Missing character 计数 = %d）' % n)

    print('PP-W01 夹具总判定：%s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
