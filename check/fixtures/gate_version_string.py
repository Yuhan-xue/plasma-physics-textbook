# -*- coding: utf-8 -*-
"""可证伪性夹具：验 check/full_check.py 的「版本串一致性」判定（E 类）真的会红。

背景：v1.26.0 发布时封面与 pdftitle 都停在 v1.25.0，已发布的 PDF 封面上印着上一版的号，
而六套件全绿——旧判定是「版本号在主文件里**任意位置**出现过就算一致」且只判 WARN，
主文件注释里留着一堆历史版本号（v1.8.0…v1.25.0），于是这条 WARN 常年飘着没人当回事。
新判定改成 pdftitle / 封面可见串 / 已编译 PDF 的 Title 元数据三处分别精确取值、分别判 FAIL。
本夹具就是那条判定的负例，防止它日后又被放松回去。

五个用例（前四个负例、最后一个正对照）：
  1) 封面版本串退回旧版           -> FAIL，点名「封面版本串 是 v0.0.1」
  2) pdftitle 退回旧版            -> FAIL，点名「pdftitle 是 v0.0.1」
  3) 封面版本串被删空             -> FAIL，点名「封面版本串 里读不到版本号」
  4) pdftitle 与封面同时退回旧版  -> 两条 FAIL 同时出现（两处独立判定，不是一处带过）
  5) 原样（正对照）               -> 三条全 PASS

做法：备份主文件 -> 逐个植入 -> 跑 `full_check.py --no-compile` -> 按级别断言消息与退出码
-> `finally` 按**原字节**恢复。恢复用 `newline=''` 写回，保证不把 LF 改写成 CRLF。

注意一：与 gate_compile_log.py 同理，本夹具断言的是 full_check 的**整体退出码**，
所以隐含要求「其它判定处于 FAIL 0」。若别处还有 FAIL，这里会报「退出码与期望不符」——
那说明套件整体尚未修好，不是本判定失效；各条目级别的消息断言始终独立有效。

注意二：正对照要求 README 版本历史表最新一行与主文件封面/pdftitle 一致。发版时若只改了
README 还没改主文件（正是 v1.26.0 的漏法），正对照会红——这是设计意图，不是夹具坏了。

注意三：本夹具植入期间若工作区有自动构建（IDE 文件监听 / 后台 latexmk）恰好重编，
`pdf/` 可能短暂留下带假版本号的 PDF，使正对照多出一条「PDF 元数据 Title」的 FAIL。
源码恢复后重编一次即自愈；判定消息本身不受影响。

依赖：无第三方库；需要 pdf/plasma_physics_textbook_v2.pdf 存在（供第 3 处比对取 Title）。
用法：python check/fixtures/gate_version_string.py
"""
import io
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAIN = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.tex')
STALE = 'v0.0.1'          # 合成的过期版本号，只需与 README 最新行不同

RE_PDFTITLE = r'pdftitle=\{[^}]*?v(\d+\.\d+\.\d+)\}'
RE_COVER = r'\\textcolor\{white!45\}\{v(\d+\.\d+\.\d+)\}'


def run_full_check():
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    r = subprocess.run([sys.executable, os.path.join('check', 'full_check.py'), '--no-compile'],
                       cwd=ROOT, capture_output=True, env=env)
    return r.stdout.decode('utf-8', errors='replace'), r.returncode


def buckets(out):
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


def build_cases(cur_pdftitle, cur_cover):
    """用主文件里的真实版本串构造植入用例，故发版后无需改本夹具。"""
    old_pt = 'pdftitle={等离子体物理自学教材 %s}' % cur_pdftitle
    old_cv = '\\textcolor{white!45}{%s}' % cur_cover
    new_pt = 'pdftitle={等离子体物理自学教材 %s}' % STALE
    new_cv = '\\textcolor{white!45}{%s}' % STALE
    return [
        ('封面退回旧版', [(old_cv, new_cv)],
         [('封面版本串 是 %s' % STALE, 'FAIL')]),
        ('pdftitle 退回旧版', [(old_pt, new_pt)],
         [('pdftitle 是 %s' % STALE, 'FAIL')]),
        ('封面版本串被删空', [(old_cv, '\\textcolor{white!45}{}')],
         [('封面版本串 里读不到版本号', 'FAIL')]),
        ('两处同时退回旧版', [(old_pt, new_pt), (old_cv, new_cv)],
         [('pdftitle 是 %s' % STALE, 'FAIL'), ('封面版本串 是 %s' % STALE, 'FAIL')]),
        ('原样（正对照）', [],
         [('pdftitle 与 README 版本历史表一致', 'PASS'),
          ('封面版本串 与 README 版本历史表一致', 'PASS'),
          ('PDF 元数据 Title 与 README 版本历史表一致', 'PASS')]),
    ]


def main():
    orig_bytes = io.open(MAIN, 'rb').read()
    orig = orig_bytes.decode('utf-8')
    m_pt = re.search(RE_PDFTITLE, orig)
    m_cv = re.search(RE_COVER, orig)
    if not m_pt or not m_cv:
        print('[X] 主文件里取不到 pdftitle / 封面版本串，夹具无法构造用例')
        return 2
    cur_pt, cur_cv = m_pt.group(1), m_cv.group(1)
    print('主文件当前版本串：pdftitle=v%s 封面=v%s' % (cur_pt, cur_cv))

    ok = True
    try:
        for name, edits, exp in build_cases('v' + cur_pt, 'v' + cur_cv):
            text = orig
            for old, new in edits:
                if old not in text:
                    print('  [X] %-18s 锚点未找到：%s' % (name, old[:60]))
                    ok = False
                    text = None
                    break
                text = text.replace(old, new, 1)
            if text is None:
                continue
            io.open(MAIN, 'w', encoding='utf-8', newline='').write(text)
            out, rc = run_full_check()
            buck = buckets(out)
            want_fail = any(lvl == 'FAIL' for _, lvl in exp)
            case_ok = True
            if (rc == 1) != want_fail:
                case_ok = False
                print('  [X] %-18s 退出码 %d 与期望%s不符'
                      % (name, rc, '（应有 FAIL）' if want_fail else '（不应有 FAIL）'))
            for needle, lvl in exp:
                if not any(needle in m for m in buck.get(lvl, [])):
                    case_ok = False
                    print('  [X] %-18s %s 级别未命中：%s' % (name, lvl, needle))
            if case_ok:
                print('  [OK] %-18s 已检（退出码 %d）' % (name, rc))
            ok = ok and case_ok
    finally:
        io.open(MAIN, 'w', encoding='utf-8', newline='').write(orig)

    back = io.open(MAIN, 'rb').read()
    if back != orig_bytes:
        print('[X] 主文件未按原字节恢复（长度差 %+d 字节，CRLF %+d 处）'
              % (len(back) - len(orig_bytes),
                 back.count(b'\r\n') - orig_bytes.count(b'\r\n')))
        ok = False
    else:
        print('[OK] 主文件已按原字节恢复')
    print('版本串夹具总判定：%s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
