# -*- coding: utf-8 -*-
"""可证伪性夹具：故意在源码里植入本轮的典型回退，确认 check/physics_regression.py 会 FAIL。

流程：备份 → 植入 16 处回退 → 跑套件（期望 退出码≠0 且 FAIL≥MIN_FAIL）→ `finally` 恢复 → 再跑（期望 0）。
其中四处打在 v1.8.0/v1.8.1 新增的附录 I/J/K 上（K.5 的阻尼率符号、K.7 的鞘厚量级与单位、K.4 的磁冻结结论、J 的逐章对照表章号），
三处打在 v1.10.0 新增的附录 M 上（束速量级、H 模阈值标度、欧姆/自加热倍数口径），
一处是 v1.10.0 的结构性回退（语义盒退回裸 `\begin{tip}`，正文首字会被吞、标签可能重复），
三处打在 v1.11.0 新增的附录 N 上（脱靶辐射功率量级、ELM 侵蚀倍数、表标签被改名）。
用来证明 P 类新检查确实会失败——不是「从来没红过」的装饰：K.5 那条只翻符号，绝对值一字不变。
恢复用 `newline=''` 写回，保证不把 LF 改写成 CRLF（否则整个文件变成「已修改」）。

对应验收：reviews/2026-10-10-independent/TODO.md 完成条件「复算本轮负例与正对照……
使『修复通过』不依赖旧错误仍存在」。

用法：python check/fixtures/gate_physics_negative.py
"""
import io
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAIN = os.path.join(ROOT, 'src', 'plasma_physics_textbook_v2.tex')
CH9 = os.path.join(ROOT, 'src', 'ch9_supplement.tex')
CH6 = os.path.join(ROOT, 'src', 'ch6_supplement.tex')
APPJ = os.path.join(ROOT, 'src', 'appendix_j_textbook_map.tex')
APPK = os.path.join(ROOT, 'src', 'appendix_k_problems.tex')
APPM = os.path.join(ROOT, 'src', 'appendix_m_heating.tex')
APPN = os.path.join(ROOT, 'src', 'appendix_n_boundary.tex')
SUITE = os.path.join(ROOT, 'check', 'physics_regression.py')

PATCHES = [
    (CH9, r'56.6\,\mathrm{ps}', r'5.66\,\mathrm{ps}'),                     # 例题数值回退
    (MAIN, '0.0543', '0.0600'),                                            # Buneman 标度回退
    (CH6, r'\omega=\frac{\omega_{pe}}{\cos\theta_0}', r'\omega=\omega_{pe}\cos\theta_0'),
    (MAIN, r'w^2\approx-0.1101', r'w^2\approx-0.2202'),                    # 双流复核例回退
    # 附录 K.5：只翻符号、量级不变——旧的「按绝对值查存在性」检查抓不到，只有带符号匹配能抓
    (APPK, r'=-2.00\times10^{-2}$', r'=+2.00\times10^{-2}$'),
    (APPK, r'=1.81\times10^{-3}\ \mathrm{m}', r'=1.81\times10^{-2}\ \mathrm{m}'),  # 鞘厚量级回退
    (APPK, r'\approx1.8\ \mathrm{mm}', r'\approx1.8\ \mathrm{cm}'),         # 鞘厚单位回退
    # 附录 K.4：把「v_A/c 小」重新当成磁冻结的充分条件（旧错误结论回归）
    (APPK, '不能据此断定磁冻结成立', '磁冻结图像在非相对论等离子体中非常好'),
    (APPJ, '第 10 章 & 约束与受控核聚变', '第 13 章 & 约束与受控核聚变'),   # 附录 J 幻章/漏行
    # 附录 M（v1.10.0）：三个不同层面的回退——数值量级、经验标度、结论口径
    (APPM, r'\approx3.10\times10^{6}\ \mathrm{m/s}', r'\approx3.10\times10^{5}\ \mathrm{m/s}'),
    (APPM, r'P_{\mathrm{LH}}\approx86.5\ \mathrm{MW}', r'P_{\mathrm{LH}}\approx8.65\ \mathrm{MW}'),
    (APPM, '差 22 倍', '差 220 倍'),
    # v1.10.0 结构性回退：两参数语义盒退回裸调用（会吞掉正文首字、并可能造成标签重复）
    (MAIN, r'\begin{tip}{磁应力记忆口诀}{tip-magstress-mnemonic}', r'\begin{tip}'),
    # 附录 N（v1.11.0）：三处不同层面——数值量级、结论里的倍数、表标签被改名（结构性）
    (APPN, r'\approx0.861', r'\approx0.761'),                              # 脱靶辐射份额量级回退
    (APPN, r'ELM 造成的侵蚀比稳态高约 $390$ 倍', r'ELM 造成的侵蚀比稳态高约 $39$ 倍'),
    (APPN, r'\label{tab:app-n-impurity}', r'\label{tab:app-n-impuritx}'),
]

MIN_FAIL = 16        # 实测值（16 处植入回退合并计入）：见下
# 逐条实测（tmp/probe_neg_v1100.py 与 tmp/probe_neg_v1110.py）每处补丁单独植入都至少产生 1 个 FAIL，
# 单独计共 19 个；合并植入时降为 16，因为 P01（附录 K 数值）、P08（附录 M 数值）与 P10（附录 N 数值）
# 各自把同一附录内的多处不符汇总成一行。故 16 是「合并植入」的实测下界，
# per-patch 可证伪性由上面的探针单独保证（v1.11.0 实测：16/16 全部命中，无落空）。


def run_suite():
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'          # 控制台是 GBK，直接跑会 UnicodeEncodeError
    p = subprocess.run([sys.executable, SUITE], capture_output=True, env=env)
    t = p.stdout.decode('utf-8', 'replace')
    nf = 0
    for line in t.split('\n'):
        if line.startswith('汇总：'):
            nf = int(line.split('FAIL')[1].split('／')[0].strip())
    if nf == 0 and '汇总' not in t:
        err = p.stderr.decode('utf-8', 'replace').strip().split('\n')[-1:]
        print('    （套件未输出汇总，疑似异常：%s）' % (err or ['?']))
    return p.returncode, nf, t


def write(path, text):
    """newline='' 保证不把 LF 改写成 CRLF（否则整个文件都变「已修改」）。"""
    io.open(path, 'w', encoding='utf-8', newline='').write(text)


def main():
    backup, original = {}, {}
    for path, _, _ in PATCHES:
        if path not in backup:
            original[path] = io.open(path, 'rb').read()
            backup[path] = original[path].decode('utf-8')
    ok = True
    try:
        for path, old, new in PATCHES:
            t = io.open(path, encoding='utf-8').read()
            if old not in t:
                print('[X] 夹具锚点未找到：%s -> %s' % (os.path.basename(path), old))
                ok = False
                continue
            io.open(path, 'w', encoding='utf-8', newline='').write(t.replace(old, new, 1))
        code, nf, _ = run_suite()
        print('[%s] 植入回退后：退出码=%d FAIL=%d（期望 !=0 且 >=%d）'
              % ('OK' if (code != 0 and nf >= MIN_FAIL) else 'X', code, nf, MIN_FAIL))
        ok = ok and code != 0 and nf >= MIN_FAIL
    finally:
        for path, t in backup.items():
            write(path, t)
    code, nf, _ = run_suite()
    print('[%s] 恢复后：退出码=%d FAIL=%d（期望 0 且 0）'
          % ('OK' if (code == 0 and nf == 0) else 'X', code, nf))
    ok = ok and code == 0 and nf == 0
    # 按**原字节**比对恢复结果：既不允许多出 CRLF，也不允许任何行尾/字节漂移。
    # （不能只查「有没有 CRLF」——core.autocrlf=true 的检出本身可能就是 CRLF。）
    for path in backup:
        b = io.open(path, 'rb').read()
        if b != original[path]:
            d = abs(len(b) - len(original[path]))
            print('[X] %s 未按原字节恢复（CRLF 变化 %+d 处，长度差 %d 字节）'
                  % (os.path.basename(path), b.count(b'\r\n') - original[path].count(b'\r\n'), d))
            ok = False
    print('夹具总判定：%s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
