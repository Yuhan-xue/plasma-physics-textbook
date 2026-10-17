# -*- coding: utf-8 -*-
"""附录 AQ「干涉与相位诊断」数值唯一来源（只 print，不写文件）。

用法：python tmp\\appAQ_calc.py | Out-File -Encoding utf8 tmp\\appAQ_calc.txt

DECLARED（声明输入，非本书推导）：
  工作点与 Q/AN/Z 同：n_e0 = 1e20 m^-3、T_e = 9.194 keV、B = 5.3 T、a = 2.0 m、R = 6.2 m
  密度剖面取抛物型 n_e(r) = n_e0 (1 - r^2/a^2)，中心弦线积分 = (4/3) n_e0 a
  三个工作波长：CO2 激光 10.6 um、远红外 118.8 um（2.52 THz）、毫米波 1.0 mm（300 GHz）
  色散干涉用 Nd:YAG 1064 nm + 倍频 532 nm（波长比恰为 2，共光路）
  弦布局 9 道：|y| = 0 / 0.45 / 0.90 / 1.35 / 1.80 m（间距 0.45 m）；对照用 21 道与 41 道
  镜面机械振动幅值 1 um；外差相位刷新率 1 MHz（等效噪声带宽 1 MHz）
  到达探测器的信号功率 1 uW（声明，含窗口与传输损耗）；系统噪声温度 3000 K
  源链：Gunn 110 GHz 100 mW -> 三倍频器效率 10% -> 330 GHz 10 mW；本振 1 mW
  混频器变频损耗 8 dB；标定不确定度 1%；弦位置不确定度 0.5 mm
  相对论修正按 T_e/m_e c^2 量级声明（不推导系数）
"""
from __future__ import print_function

import math

# ---- 物理常数 ----
E = 1.602176634e-19        # C
ME = 9.1093837015e-31      # kg
EPS0 = 8.8541878128e-12    # F/m
C = 2.99792458e8           # m/s
KB = 1.380649e-23          # J/K
RE = E * E / (4.0 * math.pi * EPS0 * ME * C * C)   # 经典电子半径 m

# ---- DECLARED ----
NE0 = 1.0e20               # m^-3
TE_EV = 9.194e3            # eV
A_MINOR = 2.0              # m
LAM_CO2 = 10.6e-6          # m
LAM_FIR = 118.8e-6         # m
LAM_MM = 1.0e-3            # m
LAM_YAG = 1064e-9          # m
LAM_SHG = 532e-9           # m
DX_VIB = 1.0e-6            # m
F_REFRESH = 1.0e6          # Hz
P_SIG = 1.0e-6             # W
T_SYS = 3000.0             # K
P_GUNN = 0.1               # W
ETA_TRIPLER = 0.10
F_GUNN = 110e9             # Hz
CL_DB = 8.0                # dB
CAL_ERR = 0.01
DY_CHORD = 0.5e-3          # m
DY_STEP = 0.45             # m
Y_CHORDS = [i * DY_STEP for i in range(5)]   # 0 / 0.45 / 0.90 / 1.35 / 1.80 m
ME_C2_EV = 510.99895e3     # eV


def p(tag, name, value, unit=''):
    print('[%s] %-44s %s %s' % (tag, name, ('%.6g' % value), unit))


def n_cutoff(lam):
    """截止密度 n_c = eps0 m_e omega^2 / e^2（omega = 2 pi c / lam）。"""
    w = 2.0 * math.pi * C / lam
    return EPS0 * ME * w * w / (E * E)


def n_line(y):
    """抛物剖面 n_e(r) = n_e0 (1 - r^2/a^2) 在碰撞参数 y 处的线积分。"""
    return (4.0 / 3.0) * NE0 * (A_MINOR * A_MINOR - y * y) ** 1.5 / (A_MINOR * A_MINOR)


def peel(nvals, ysv, amin):
    """洋葱剥皮反演：壳层边界取各弦碰撞参数，最外层壳延伸到 a。

    返回每个壳层的等效密度（n[0] 为中心壳层），壳层假设内为常数。
    """
    m = len(ysv)
    r = list(ysv) + [amin]
    n = [0.0] * m
    for i in range(m - 1, -1, -1):
        s = 0.0
        for j in range(i + 1, m):
            s += 2.0 * n[j] * (math.sqrt(r[j + 1] ** 2 - ysv[i] ** 2)
                               - math.sqrt(r[j] ** 2 - ysv[i] ** 2))
        n[i] = (nvals[i] - s) / (2.0 * math.sqrt(r[i + 1] ** 2 - ysv[i] ** 2))
    return n


def chord_set(step):
    ys = []
    y = 0.0
    while y < A_MINOR - 1e-9:
        ys.append(y)
        y += step
    return ys


print('=== [AQ.1] 相移、条纹与三个波长 ===')
p('AQ.1', '经典电子半径 r_e', RE, 'm')
p('AQ.1', '中心弦线积分密度 N_L(0)', n_line(0.0), 'm^-2')
p('AQ.1', '中心弦等效密度 N_L/(2a)', n_line(0.0) / (2.0 * A_MINOR), 'm^-3')
for lam, name in ((LAM_CO2, 'CO2 10.6 um'), (LAM_FIR, 'FIR 118.8 um'), (LAM_MM, 'mm 1.0 mm')):
    f = C / lam
    nc = n_cutoff(lam)
    dphi = RE * lam * n_line(0.0)
    p('AQ.1', '%s 频率' % name, f, 'Hz')
    p('AQ.1', '%s 截止密度' % name, nc, 'm^-3')
    p('AQ.1', '%s 工作密度占比 n_e0/n_c' % name, NE0 / nc)
    p('AQ.1', '%s 折射率亏损 1-N' % name, NE0 / (2.0 * nc))
    p('AQ.1', '%s 中心弦相移' % name, dphi, 'rad')
    p('AQ.1', '%s 中心弦条纹数' % name, dphi / (2.0 * math.pi), 'fringe')
    p('AQ.1', '%s 一条条纹对应的线积分密度' % name, 2.0 * math.pi / (RE * lam), 'm^-2')

print('')
print('=== [AQ.2] 相位精度、外差与振动 ===')
for lam, name in ((LAM_CO2, 'CO2 10.6 um'), (LAM_FIR, 'FIR 118.8 um'), (LAM_MM, 'mm 1.0 mm')):
    p('AQ.2', '%s 每 0.001 rad 对应的线积分密度' % name, 1e-3 / (RE * lam), 'm^-2')
    p('AQ.2', '%s 每 0.001 rad 对应的电子密度' % name,
      1e-3 / (RE * lam) / (n_line(0.0) / NE0), 'm^-3')
    p('AQ.2', '%s 1 um 振动引入的相位' % name, 2.0 * math.pi * DX_VIB / lam, 'rad')
    p('AQ.2', '%s 1 um 振动折合的条纹' % name, DX_VIB / lam, 'fringe')
    p('AQ.2', '%s 相位误差控制在 0.001 rad 的允许振幅' % name,
      1e-3 * lam / (2.0 * math.pi), 'm')
    p('AQ.2', '%s 相位误差控制在 0.001 rad 的允许振幅' % name,
      1e-3 * lam / (2.0 * math.pi) / 1e-9, 'nm')
p_if = P_SIG * 10.0 ** (-CL_DB / 10.0)
p('AQ.2', '混频器输出中频功率', p_if, 'W')
snr = p_if / (KB * T_SYS * F_REFRESH)
p('AQ.2', '中频信噪比（带宽 1 MHz）', snr)
sig_phi = 1.0 / math.sqrt(2.0 * snr)
p('AQ.2', '窄带相位噪声 sigma_phi', sig_phi, 'rad')
p('AQ.2', '相位噪声折合的线积分密度', sig_phi / (RE * LAM_FIR), 'm^-2')
p('AQ.2', '相位噪声折合的电子密度', sig_phi / (RE * LAM_FIR) / (n_line(0.0) / NE0), 'm^-3')
p('AQ.2', '解缠上限：1 us 内允许的线积分密度变化',
  math.pi / (1e-6 * RE * LAM_FIR), 'm^-2/s')
p('AQ.2', '解缠上限折合的电子密度变化率',
  math.pi / (1e-6 * RE * LAM_FIR) / (n_line(0.0) / NE0), 'm^-3/s')

print('')
print('=== [AQ.3] 色散干涉与双波长 ===')
d_fund = RE * LAM_YAG * n_line(0.0)
vib_yag = 2.0 * math.pi * DX_VIB / LAM_YAG
p('AQ.3', '基频 1064 nm 的色散相移 D', d_fund, 'rad')
p('AQ.3', '倍频 532 nm 的色散相移 2D', 2.0 * d_fund, 'rad')
p('AQ.3', '1 um 振动在基频上的相位 V', vib_yag, 'rad')
p('AQ.3', '实测基频相移 D + V', d_fund + vib_yag, 'rad')
p('AQ.3', '实测倍频相移 2D + V', 2.0 * d_fund + vib_yag, 'rad')
p('AQ.3', '组合量 (2f - f) 复原的色散相移', (2.0 * d_fund + vib_yag) - (d_fund + vib_yag), 'rad')
p('AQ.3', '组合量 (2f - f) 留下的振动残差',
  abs((2.0 * d_fund + vib_yag) - (d_fund + vib_yag) - d_fund), 'rad')
p('AQ.3', '单波长：1 um 振动折合的中心密度误差',
  vib_yag / (RE * LAM_YAG) / (n_line(0.0) / NE0), 'm^-3')
p('AQ.3', '单波长：1 um 振动折合的中心密度误差（相对工作密度）',
  vib_yag / (RE * LAM_YAG) / (n_line(0.0) / NE0) / NE0)
p('AQ.3', '色散干涉：相位噪声折合的中心密度误差',
  sig_phi / (RE * LAM_YAG) / (n_line(0.0) / NE0), 'm^-3')
p('AQ.3', '色散干涉：相位噪声折合的中心密度相对误差',
  sig_phi / (RE * LAM_YAG) / (n_line(0.0) / NE0) / NE0)
p('AQ.3', '振动被压低的倍数', vib_yag / sig_phi)
p('AQ.3', 'CO2 与 FIR 的波长比（双色干涉的色散比）', LAM_FIR / LAM_CO2)
p('AQ.3', 'FIR 与 CO2 双色独立测密度的相位噪声折合误差',
  sig_phi / (RE * LAM_CO2) / (n_line(0.0) / NE0), 'm^-3')

print('')
print('=== [AQ.4] 阿贝尔反演与弦布局 ===')
for y in Y_CHORDS:
    p('AQ.4', 'y = %.2f m 弦的线积分密度' % y, n_line(y), 'm^-2')
p('AQ.4', '弦间距', DY_STEP, 'm')
p('AQ.4', '中心道与最外道线积分之比', n_line(0.0) / n_line(1.80))
n9 = peel([n_line(y) for y in Y_CHORDS], Y_CHORDS, A_MINOR)
p('AQ.4', '9 道剥皮反演的中心壳层密度', n9[0], 'm^-3')
p('AQ.4', '9 道反演相对真值的偏差', n9[0] / NE0 - 1.0)
for step in (0.1, 0.05):
    ys = chord_set(step)
    nf = peel([n_line(y) for y in ys], ys, A_MINOR)
    p('AQ.4', '%d 道剥皮反演的中心壳层密度' % len(ys), nf[0], 'm^-3')
    p('AQ.4', '%d 道反演相对真值的偏差' % len(ys), nf[0] / NE0 - 1.0)
ys41 = chord_set(0.05)
n41 = peel([n_line(y) for y in ys41], ys41, A_MINOR)
n9_outer = peel([n_line(y) * (1.01 if abs(y - 1.80) < 1e-9 else 1.0) for y in Y_CHORDS],
                Y_CHORDS, A_MINOR)
p('AQ.4', '9 道：最外道偏 1% 时中心密度的变化', n9_outer[0] - n9[0], 'm^-3')
p('AQ.4', '9 道：最外道偏 1% 折合中心密度的相对变化', (n9_outer[0] - n9[0]) / n9[0])
n9_center = peel([n_line(y) * (1.01 if y == 0.0 else 1.0) for y in Y_CHORDS],
                 Y_CHORDS, A_MINOR)
p('AQ.4', '9 道：中心道偏 1% 时中心密度的变化', n9_center[0] - n9[0], 'm^-3')
n41_outer = peel([n_line(y) * (1.01 if abs(y - ys41[-1]) < 1e-9 else 1.0) for y in ys41],
                 ys41, A_MINOR)
p('AQ.4', '41 道：最外道偏 1% 折合中心密度的相对变化', (n41_outer[0] - n41[0]) / n41[0])
p('AQ.4', '弦位置 0.5 mm 在 y=1.35 m 处引起的线积分变化',
  abs(n_line(1.35 + DY_CHORD) - n_line(1.35)), 'm^-2')
p('AQ.4', '弦位置 0.5 mm 相对该道线积分的份额',
  abs(n_line(1.35 + DY_CHORD) - n_line(1.35)) / n_line(1.35))
p('AQ.4', '弦位置 0.5 mm 在最外道引起的相对变化',
  abs(n_line(1.80 + DY_CHORD) - n_line(1.80)) / n_line(1.80))
p('AQ.4', '最外层壳层厚度（9 道）', A_MINOR - 1.80, 'm')

print('')
print('=== [AQ.5] 源与探测器链 ===')
p('AQ.5', 'Gunn 源功率', P_GUNN, 'W')
p('AQ.5', '三倍频后功率', P_GUNN * ETA_TRIPLER, 'W')
f_mm = F_GUNN * 3.0
p('AQ.5', '三倍频后频率', f_mm, 'Hz')
lam_mm_src = C / f_mm
p('AQ.5', '对应波长', lam_mm_src, 'm')
p('AQ.5', '该频率的截止密度', n_cutoff(lam_mm_src), 'm^-3')
p('AQ.5', '截止密度相对工作密度的余量', n_cutoff(lam_mm_src) / NE0)
p('AQ.5', '传输 3 dB 后到达等离子体的功率', P_GUNN * ETA_TRIPLER * 0.5, 'W')
p('AQ.5', '本振功率', 1e-3, 'W')
p('AQ.5', '混频器变频损耗', CL_DB, 'dB')
p('AQ.5', '系统噪声温度', T_SYS, 'K')
p('AQ.5', '热噪声功率 kT B（1 MHz）', KB * T_SYS * F_REFRESH, 'W')
p('AQ.5', '接收功率与本振功率之比', P_SIG / 1e-3)
p('AQ.5', '相位刷新率对应的时间分辨', 1.0 / F_REFRESH, 's')
p('AQ.5', '本振与信号的功率比（dB）', 10.0 * math.log10(1e-3 / P_SIG), 'dB')

print('')
print('=== [AQ.6] 相对论修正与误差预算 ===')
rel = TE_EV / ME_C2_EV
p('AQ.6', 'T_e / m_e c^2', rel)
p('AQ.6', '相对论修正折合的密度偏置', rel * NE0, 'm^-3')
p('AQ.6', '毫米波下线性近似的相对偏差', 
  (1.0 - math.sqrt(1.0 - NE0 / n_cutoff(LAM_MM))) / (NE0 / (2.0 * n_cutoff(LAM_MM))) - 1.0)
p('AQ.6', 'CO2 下线性近似的相对偏差',
  (1.0 - math.sqrt(1.0 - NE0 / n_cutoff(LAM_CO2))) / (NE0 / (2.0 * n_cutoff(LAM_CO2))) - 1.0)
items = [
    ('相位噪声', sig_phi / (RE * LAM_FIR) / (n_line(0.0) / NE0) / NE0),
    ('标定', CAL_ERR),
    ('弦位置', abs(n_line(1.80 + DY_CHORD) - n_line(1.80)) / n_line(1.80)),
    ('相对论修正', rel),
]
tot = 0.0
for name, val in items:
    p('AQ.6', '误差项：%s' % name, val)
    tot += val * val
p('AQ.6', '合成相对误差（RSS）', math.sqrt(tot))
p('AQ.6', '合成后的中心密度绝对误差', math.sqrt(tot) * NE0, 'm^-3')
p('AQ.6', '只算相位噪声时的相对误差', sig_phi / (RE * LAM_FIR) / (n_line(0.0) / NE0) / NE0)

print('')
print('=== [AQ.7] 练习解答 ===')
p('AQ.7', '练习 1：118.8 um、0.001 rad 对应的电子密度',
  1e-3 / (RE * LAM_FIR) / (n_line(0.0) / NE0), 'm^-3')
p('AQ.7', '练习 1：相对份额', (1e-3 / (RE * LAM_FIR) / (n_line(0.0) / NE0)) / NE0)
p('AQ.7', '练习 2：色散干涉复原的 D', (2.0 * d_fund + vib_yag) - (d_fund + vib_yag), 'rad')
p('AQ.7', '练习 2：留下的振动相位 V', 2.0 * (d_fund + vib_yag) - (2.0 * d_fund + vib_yag), 'rad')
p('AQ.7', '练习 3：SNR 与相位噪声', sig_phi, 'rad')
p('AQ.7', '练习 4：0.01 条纹对应的允许振幅（10.6 um）', 0.01 * LAM_CO2, 'm')
p('AQ.7', '练习 5：n_e/n_c = 0.1 时线性近似的相对偏差',
  (1.0 - math.sqrt(0.9)) / 0.05 - 1.0)
p('AQ.7', '练习 6：最外道偏 1% 折合中心密度的相对变化', (n9_outer[0] - n9[0]) / n9[0])
p('AQ.7', '练习 6：中心道偏 1% 时中心密度的变化', n9_center[0] - n9[0], 'm^-3')
p('AQ.7', '练习 6：中心道偏 1% 折合中心密度的相对变化', (n9_center[0] - n9[0]) / n9[0])
p('AQ.7', '练习 6：9 道反演的弦数不足偏差', n9[0] / NE0 - 1.0)

print('')
print('=== [AQ.9] 百分数对照（供正文表格直接引用） ===')
p('AQ.9', '相位噪声折合的相对密度误差（百分数）',
  sig_phi / (RE * LAM_FIR) / (n_line(0.0) / NE0) / NE0 * 100.0, '%')
p('AQ.9', '色散干涉相对误差（百分数）',
  sig_phi / (RE * LAM_YAG) / (n_line(0.0) / NE0) / NE0 * 100.0, '%')
p('AQ.9', '9 道反演相对真值的偏差（百分数）', n9[0] / NE0 * 100.0 - 100.0, '%')
p('AQ.9', '20 道反演相对真值的偏差（百分数）',
  peel([n_line(y) for y in chord_set(0.1)], chord_set(0.1), A_MINOR)[0] / NE0 * 100.0 - 100.0, '%')
p('AQ.9', '40 道反演相对真值的偏差（百分数）', n41[0] / NE0 * 100.0 - 100.0, '%')
p('AQ.9', '最外道偏 1% 折合中心密度变化（百分数）', (n9_outer[0] - n9[0]) / n9[0] * 100.0, '%')
p('AQ.9', '中心道偏 1% 折合中心密度变化（百分数）', (n9_center[0] - n9[0]) / n9[0] * 100.0, '%')
p('AQ.9', '弦位置 0.5 mm 在 y=1.35 m 的份额（百分数）',
  abs(n_line(1.35 + DY_CHORD) - n_line(1.35)) / n_line(1.35) * 100.0, '%')
p('AQ.9', '弦位置 0.5 mm 在最外道的份额（百分数）',
  abs(n_line(1.80 + DY_CHORD) - n_line(1.80)) / n_line(1.80) * 100.0, '%')
p('AQ.9', '相对论修正幅度（百分数）', rel * 100.0, '%')
p('AQ.9', '毫米波线性近似偏差（百分数）',
  ((1.0 - math.sqrt(1.0 - NE0 / n_cutoff(LAM_MM))) / (NE0 / (2.0 * n_cutoff(LAM_MM))) - 1.0) * 100.0, '%')
p('AQ.9', '合成相对误差（百分数）', math.sqrt(tot) * 100.0, '%')
p('AQ.9', '标定误差项（百分数）', CAL_ERR * 100.0, '%')
p('AQ.9', '相位噪声误差项（百分数）',
  sig_phi / (RE * LAM_FIR) / (n_line(0.0) / NE0) / NE0 * 100.0, '%')
p('AQ.9', '弦位置误差项（百分数）',
  abs(n_line(1.80 + DY_CHORD) - n_line(1.80)) / n_line(1.80) * 100.0, '%')
p('AQ.9', 'n_e/n_c = 0.1 时线性近似偏差（百分数）',
  ((1.0 - math.sqrt(0.9)) / 0.05 - 1.0) * 100.0, '%')
p('AQ.9', '工作密度占比（CO2，百分数）', NE0 / n_cutoff(LAM_CO2) * 100.0, '%')
p('AQ.9', '工作密度占比（FIR，百分数）', NE0 / n_cutoff(LAM_FIR) * 100.0, '%')
p('AQ.9', '工作密度占比（毫米波，百分数）', NE0 / n_cutoff(LAM_MM) * 100.0, '%')
p('AQ.9', '折射率亏损（CO2，百分数）', NE0 / (2.0 * n_cutoff(LAM_CO2)) * 100.0, '%')
p('AQ.9', '折射率亏损（FIR，百分数）', NE0 / (2.0 * n_cutoff(LAM_FIR)) * 100.0, '%')
p('AQ.9', '折射率亏损（毫米波，百分数）', NE0 / (2.0 * n_cutoff(LAM_MM)) * 100.0, '%')

print('')
print('=== [AQ.8] 文献级对照（声明） ===')
p('AQ.8', 'ITER 密度上限量级（声明）', 1.2e20, 'm^-3')
p('AQ.8', '色散干涉常用基频（声明）', LAM_YAG, 'm')
p('AQ.8', '远红外干涉仪常用波长（声明）', LAM_FIR, 'm')
p('AQ.8', '弦数典型值（声明）', 9, '道')

# [AQ.10] 练习补充派生量（自足计算，不依赖上方变量名）
import math as _m
_P_IF = 1.58489e-7
_T_SYS = 3000.0
_B_IF = 1.0e6
_snr = _P_IF / (1.380649e-23 * _T_SYS * _B_IF)
print('=== [AQ.10] 练习补充派生量 ===')
p('AQ.10', '0.01 条纹对应的相位', 0.01 * 2 * _m.pi, 'rad')
p('AQ.10', '0.01 条纹允许的 CO2 振动振幅', 0.01 * 2 * _m.pi / 5.92753e5, 'm')
p('AQ.10', '0.01 条纹允许的 CO2 振动振幅（微米）', 0.01 * 2 * _m.pi / 5.92753e5 * 1e6, 'um')
p('AQ.10', 'x/4 粗略估计（x = 0.1）', 0.1 / 4.0, '')
p('AQ.10', '中频信噪比（dB）', 10 * _m.log10(_snr), 'dB')

# [AQ.11] 剖面局域密度与源频率（练习与 AQ.5 正文引用）
print('=== [AQ.11] 局域密度与源频率 ===')
p('AQ.11', 'y = 1.35 m 处局域密度', 1.0e20 * (1.0 - (1.35 / 2.0) ** 2), 'm^-3')
p('AQ.11', 'y = 1.80 m 处局域密度', 1.0e20 * (1.0 - (1.80 / 2.0) ** 2), 'm^-3')
p('AQ.11', 'Gunn 源频率', 110.0, 'GHz')
p('AQ.11', '三倍频输出频率', 330.0, 'GHz')
