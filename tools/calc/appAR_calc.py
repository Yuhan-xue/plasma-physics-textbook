# -*- coding: utf-8 -*-
"""附录 AR「电源与储能：从脉冲电源到电网接口」的唯一数值来源（只 print，不落盘）。

用法：python tmp/appAR_calc.py
工作点沿用 D1（与附录 M/S/U/V/W/X/Y/Z/AA/.../AQ 同一组）：
  n_e=1e20 m^-3, T_e=9.194 keV, B=5.3 T, a=2.0 m, R=6.2 m, V=832.2 m^3,
  Z_eff=2, P_fus=556.5 MW, P_alpha=111.3 MW, W_th=371.4 MJ, tau_CQ=76.9074 ms
声明输入（DECLARED）：
  eta_th=0.40（热-电转换）, eta_conv=0.97（电源变换效率）
  P_CD=50.0 MW（加热与电流驱动电功率）, P_cryo=13.4157 MW（= 附录 AM 的 1.31333+12.1024）
  P_plant=8.0 MW（泵/加料/诊断/控制等厂用电）
  六脉波纹波系数 13.4%、十二脉波 4.2%（文献级）
  飞轮：m=20000 kg, r=2.0 m, n=3000 rpm, rho_steel=7850 kg/m^3, sigma_steel=1.0 GPa,
        rho_cf=1600 kg/m^3, sigma_cf=3.0 GPa
  电容：E=2.0 MJ, V=20.0 kV（另一例 C=5 mF）
  电感储能：L=10 H, I=40 kA
  十二脉波换相角 mu=10 deg，触发角 alpha=15 deg
  电网：S_grid=1000 MVA, H_inertia=5 s, 扰动 100 MW
"""

M = 1.0

def p(tag, name, val, unit='', extra=''):
    if isinstance(val, float):
        s = ('%.6g' % val)
    else:
        s = str(val)
    print('[%s] %-34s = %s %s%s' % (tag, name, s, unit, ('  ' + extra) if extra else ''))

def main():
    # ---------- 公用 ----------
    P_fus = 556.5e6
    P_alpha = 111.3e6
    P_th = P_fus
    eta_th = 0.40
    eta_conv = 0.97
    P_CD = 50.0e6
    P_cryo = 13.4157e6
    P_plant = 8.0e6
    p('AR.1', 'P_fus（MW 口径）', P_fus / 1e6, 'MW')
    p('AR.1', 'P_CD（加热与电流驱动）', P_CD, 'W', '= %.6g MW' % (P_CD / 1e6))
    p('AR.1', 'P_cryo（低温，附录 AM 口径）', P_cryo, 'W', '= %.6g MW' % (P_cryo / 1e6))
    p('AR.1', 'P_plant（其他厂用电）', P_plant, 'W', '= %.6g MW' % (P_plant / 1e6))

    print('==== [AR.1] 三本能量账：脉冲爬升、稳态环流、瞬态转移 ====')
    p('AR.1', 'P_fus', P_fus, 'W')
    p('AR.1', 'P_alpha', P_alpha, 'W')
    p('AR.1', 'P_th（全部入包层）', P_th, 'W')
    P_gross = eta_th * P_th
    p('AR.1', 'P_gross = eta_th*P_th', P_gross, 'W', '= %.6g MW' % (P_gross / 1e6))
    P_deliv = P_CD + P_cryo + P_plant
    p('AR.1', 'P_deliv（送到负载）', P_deliv, 'W', '= %.6g MW' % (P_deliv / 1e6))
    P_rec = P_deliv / eta_conv
    p('AR.1', 'P_rec = P_deliv/eta_conv', P_rec, 'W', '= %.6g MW' % (P_rec / 1e6))
    P_closs = P_rec - P_deliv
    p('AR.1', '变换损耗', P_closs, 'W', '= %.6g MW' % (P_closs / 1e6))
    f_rec = P_rec / P_gross
    p('AR.1', '环流份额 f_rec', f_rec, '', '= %.6g %%' % (100 * f_rec))
    P_net = P_gross - P_rec
    p('AR.1', '净电 P_net', P_net, 'W', '= %.6g MW' % (P_net / 1e6))
    eta_net = P_net / P_th
    p('AR.1', '电厂净效率 eta_net', eta_net, '', '= %.6g %%' % (100 * eta_net))
    p('AR.1', '每 1 MW 净电对应的聚变功率', P_th / P_net, 'MW/MW')
    # 脉冲爬升（中心螺线管）
    L_CS = 10.0
    I_pk = 40.0e3
    E_CS = 0.5 * L_CS * I_pk ** 2
    p('AR.1', 'CS 磁储能 0.5*L*I^2', E_CS, 'J', '= %.6g GJ' % (E_CS / 1e9))
    t_ramp = 10.0
    P_ramp = E_CS / t_ramp
    p('AR.1', '10 s 爬升平均功率', P_ramp, 'W', '= %.6g MW' % (P_ramp / 1e6))
    dI = 2 * I_pk / t_ramp
    p('AR.1', '爬升速率 dI/dt', dI, 'A/s')
    V_CS = L_CS * dI
    p('AR.1', 'CS 端电压 L*dI/dt', V_CS, 'V', '= %.6g kV' % (V_CS / 1e3))
    t_ramp2 = 100.0
    p('AR.1', '100 s 爬升平均功率', E_CS / t_ramp2, 'W', '= %.6g MW' % (E_CS / t_ramp2 / 1e6))
    p('AR.1', '爬升功率/净电（10 s）', P_ramp / P_net, '')
    # 瞬态转移（破裂缓解）
    E_mit = 2.0e6
    t_mit = 10.0e-3
    P_mit = E_mit / t_mit
    p('AR.1', '破裂缓解脉冲功率（2 MJ/10 ms）', P_mit, 'W', '= %.6g MW' % (P_mit / 1e6))
    p('AR.1', '瞬态脉冲/稳态环流', P_mit / P_rec, '')

    print('==== [AR.2] 变流器与纹波：十二脉波为什么是标准 ====')
    ripple6 = 0.134
    ripple12 = 0.042
    p('AR.2', '六脉波纹波系数（文献）', ripple6, '', '= %.6g %%' % (100 * ripple6))
    p('AR.2', '十二脉波纹波系数（文献）', ripple12, '', '= %.6g %%' % (100 * ripple12))
    p('AR.2', '纹波比（六/十二）', ripple6 / ripple12, '')
    alpha = 15.0
    import math
    ar = math.radians(alpha)
    p('AR.2', '触发角 alpha', alpha, 'deg')
    p('AR.2', 'cos(alpha)', math.cos(ar), '')
    p('AR.2', '位移因数 tan(alpha)', math.tan(ar), '')
    Q_conv = math.tan(ar) * P_gross
    p('AR.2', '无功需求 Q = P*tan(alpha)', Q_conv, 'var', '= %.6g MVAr' % (Q_conv / 1e6))
    mu = 10.0
    p('AR.2', '换相角 mu', mu, 'deg')
    p('AR.2', 'cos(alpha+mu/2)', math.cos(math.radians(alpha + mu / 2)), '')
    p('AR.2', '功率因数下降（相对 cos alpha）', 1 - math.cos(math.radians(alpha + mu / 2)) / math.cos(ar), '', '= %.6g %%' % (100 * (1 - math.cos(math.radians(alpha + mu / 2)) / math.cos(ar))))
    for n in (5, 7, 11, 13):
        p('AR.2', '谐波 %d 次相对幅值 1/n' % n, 1.0 / n, '')
    p('AR.2', '11 与 13 次合成', 1.0 / 11 + 1.0 / 13, '')
    U_d0 = 24.0e3
    p('AR.2', '理想空载直流电压 U_d0', U_d0, 'V', '= %.6g kV' % (U_d0 / 1e3))
    p('AR.2', '带触发角输出 U_d0*cos(alpha)', U_d0 * math.cos(ar), 'V', '= %.6g kV' % (U_d0 * math.cos(ar) / 1e3))
    I_dc = 40.0e3
    p('AR.2', '直流电流 I_dc', I_dc, 'A', '= %.6g kA' % (I_dc / 1e3))
    P_dc = U_d0 * math.cos(ar) * I_dc
    p('AR.2', '直流功率', P_dc, 'W', '= %.6g MW' % (P_dc / 1e6))
    p('AR.2', '交流侧视在功率（cos(alpha+mu/2) 口径）', P_dc / math.cos(math.radians(alpha + mu / 2)), 'VA', '= %.6g MVA' % (P_dc / math.cos(math.radians(alpha + mu / 2)) / 1e6))
    p('AR.2', '纹波电压（十二脉波，U_d0 口径）', U_d0 * ripple12, 'V', '= %.6g V' % (U_d0 * ripple12))

    print('==== [AR.3] 储能：飞轮、电容与它们能撑多久 ====')
    m_fw = 20000.0
    r_fw = 2.0
    rpm = 3000.0
    omega = 2 * math.pi * rpm / 60.0
    p('AR.3', '飞轮角速度 omega', omega, 'rad/s', '= %.6g rpm' % rpm)
    I_fw = 0.5 * m_fw * r_fw ** 2
    p('AR.3', '实心圆盘转动惯量 I', I_fw, 'kg m^2')
    E_fw = 0.5 * I_fw * omega ** 2
    p('AR.3', '飞轮储能 0.5*I*omega^2', E_fw, 'J', '= %.6g GJ' % (E_fw / 1e9))
    v_tip = omega * r_fw
    p('AR.3', '边缘线速度', v_tip, 'm/s')
    rho_s = 7850.0
    sigma_s = 1.0e9
    p('AR.3', '钢的环向应力 rho*v^2', rho_s * v_tip ** 2, 'Pa', '= %.6g GPa' % (rho_s * v_tip ** 2 / 1e9))
    v_max_s = math.sqrt(sigma_s / rho_s)
    p('AR.3', '钢允许线速度 sqrt(sigma/rho)', v_max_s, 'm/s')
    p('AR.3', '对应允许转速（r=2 m）', v_max_s / r_fw * 60 / (2 * math.pi), 'rpm')
    rho_cf = 1600.0
    sigma_cf = 3.0e9
    v_max_cf = math.sqrt(sigma_cf / rho_cf)
    p('AR.3', '碳纤维允许线速度', v_max_cf, 'm/s')
    p('AR.3', '碳纤维/钢 允许线速度比', v_max_cf / v_max_s, '')
    p('AR.3', '同质量同半径碳纤维飞轮储能', 0.5 * I_fw * (v_max_cf / r_fw) ** 2, 'J', '= %.6g GJ' % (0.5 * I_fw * (v_max_cf / r_fw) ** 2 / 1e9))
    P_fw = 100.0e6
    p('AR.3', '100 MW 放电持续时间', E_fw / P_fw, 's')
    # 电容
    E_cap = 2.0e6
    V_cap = 20.0e3
    C_cap = 2 * E_cap / V_cap ** 2
    p('AR.3', '电容储能所需 C = 2E/V^2', C_cap, 'F', '= %.6g mF' % (C_cap * 1e3))
    C_given = 5.0e-3
    p('AR.3', 'C=5 mF 时电压 sqrt(2E/C)', math.sqrt(2 * E_cap / C_given), 'V', '= %.6g kV' % (math.sqrt(2 * E_cap / C_given) / 1e3))
    p('AR.3', '2 MJ 在 10 ms 释放的功率', E_cap / 10e-3, 'W', '= %.6g MW' % (E_cap / 10e-3 / 1e6))
    p('AR.3', '20 kV 下释放 2 MJ 的平均电流', E_cap / 10e-3 / V_cap, 'A', '= %.6g kA' % (E_cap / 10e-3 / V_cap / 1e3))
    p('AR.3', '峰值电流（线性放电、10 ms、C 口径 I=CV/t）', C_cap * V_cap / 10e-3, 'A', '= %.6g MA' % (C_cap * V_cap / 10e-3 / 1e6))
    p('AR.3', '电容储能密度（按 2 MJ / 50 L）', E_cap / 50e-3, 'J/m^3', '= %.6g MJ/m^3' % (E_cap / 50e-3 / 1e6))
    # 电感储能
    L_es = 10.0
    I_es = 40.0e3
    E_es = 0.5 * L_es * I_es ** 2
    p('AR.3', '电感储能 0.5*L*I^2', E_es, 'J', '= %.6g GJ' % (E_es / 1e9))
    p('AR.3', '电感储能断开时间常数（R=1 mOhm）', L_es / 1e-3, 's')
    # 电池口径
    p('AR.3', '1 MWh', 3.6e9, 'J')
    p('AR.3', '飞轮等效 MWh', E_fw / 3.6e9, 'MWh')
    p('AR.3', '厂用电 8 MW 撑 1 h 所需能量', 8e6 * 3600, 'J', '= %.6g GJ' % (8e6 * 3600 / 1e9))

    print('==== [AR.4] 环流功率：电厂净效率的分母 ====')
    p('AR.4', 'CD 份额', P_CD / P_rec, '', '= %.6g %%' % (100 * P_CD / P_rec))
    p('AR.4', '低温份额', P_cryo / P_rec, '', '= %.6g %%' % (100 * P_cryo / P_rec))
    p('AR.4', '厂用电份额', P_plant / P_rec, '', '= %.6g %%' % (100 * P_plant / P_rec))
    p('AR.4', '变换损耗份额', P_closs / P_rec, '', '= %.6g %%' % (100 * P_closs / P_rec))
    dP = 5.0e6
    p('AR.4', 'CD 省 5 MW 后的 P_rec', P_rec - dP, 'W', '= %.6g MW' % ((P_rec - dP) / 1e6))
    p('AR.4', '对应净效率增量', (P_gross - (P_rec - dP)) / P_th - eta_net, '', '= %.6g 个百分点' % (100 * ((P_gross - (P_rec - dP)) / P_th - eta_net)))
    p('AR.4', '若 CD 电效率提高 10% 省下的功率', 0.10 * P_CD, 'W', '= %.6g MW' % (0.10 * P_CD / 1e6))
    p('AR.4', '环流份额降到 25% 所需 P_rec', 0.25 * P_gross, 'W', '= %.6g MW' % (0.25 * P_gross / 1e6))
    p('AR.4', '需要削减的环流功率', P_rec - 0.25 * P_gross, 'W', '= %.6g MW' % ((P_rec - 0.25 * P_gross) / 1e6))
    p('AR.4', '工程增益 Q_eng = P_net/P_rec', P_net / P_rec, '')
    p('AR.4', '科学增益 Q_pl = P_fus/P_CD', P_fus / P_CD, '')

    print('==== [AR.5] 电网接口：并网、频率响应与无功 ====')
    S_grid = 1000.0e6
    H_in = 5.0
    dP_g = 100.0e6
    rocof = dP_g / (2 * H_in * S_grid)
    p('AR.5', '电网惯性常数 H', H_in, 's')
    p('AR.5', '短路容量 S_grid', S_grid, 'VA', '= %.6g MVA' % (S_grid / 1e6))
    p('AR.5', '扰动 dP', dP_g, 'W', '= %.6g MW' % (dP_g / 1e6))
    p('AR.5', '频率变化率 df/dt', rocof, 'Hz/s')
    p('AR.5', '允许 0.2 Hz 跌落的时间', 0.2 / rocof, 's')
    p('AR.5', '厂用负荷占电网短路容量比', P_deliv / S_grid, '', '= %.6g %%' % (100 * P_deliv / S_grid))
    p('AR.5', '无功补偿容量（按 Q_conv）', Q_conv, 'var', '= %.6g MVAr' % (Q_conv / 1e6))
    p('AR.5', '补偿后功率因数', math.cos(ar), '')
    f0 = 50.0
    for n in (11, 13):
        Lf = 1.0 / ((2 * math.pi * n * f0) ** 2 * 100e-6)
        p('AR.5', '%d 次谐波滤波器电感（C=100 uF）' % n, Lf, 'H', '= %.6g mH' % (Lf * 1e3))
    p('AR.5', '净电送出占短路容量比', P_net / S_grid, '', '= %.6g %%' % (100 * P_net / S_grid))
    p('AR.5', '黑启动所需能量（厂用电 8 MW，1 h）', 8e6 * 3600, 'J', '= %.6g GJ' % (8e6 * 3600 / 1e9))

    print('==== [AR.6] 保护与故障：毫秒级动作 ====')
    L_p = 5.0e-6
    dI_p = 1.0e5
    t_fast = 1.0e-3
    p('AR.6', '等离子体电感 L_p', L_p, 'H', '= %.6g uH' % (L_p * 1e6))
    p('AR.6', '快控需要在 1 ms 内改变电流', dI_p, 'A', '= %.6g kA' % (dI_p / 1e3))
    V_fast = L_p * dI_p / t_fast
    p('AR.6', '快控所需电压 L*dI/dt', V_fast, 'V', '= %.6g kV' % (V_fast / 1e3))
    p('AR.6', '快控所需功率（1 ms 平均）', V_fast * dI_p / 2, 'W', '= %.6g MW' % (V_fast * dI_p / 2 / 1e6))
    R_coil = 1.0e-3
    p('AR.6', '线圈失电衰减时间常数 L/R', L_CS / R_coil, 's', '= %.6g h' % (L_CS / R_coil / 3600))
    p('AR.6', '失电后 1 s 残余电流份额 exp(-t/tau)', math.exp(-1.0 / (L_CS / R_coil)), '')
    p('AR.6', '保护动作时间（SPI 阀门）', 1.0e-3, 's', '= 1 ms')
    p('AR.6', '在热沉积时间尺度 tau_CQ=76.9074 ms 内的动作余量', 76.9074e-3 / 1.0e-3, '')
    p('AR.6', '破裂缓解能量 2 MJ / 环流功率', 2.0e6 / P_rec, 's')
    p('AR.6', '破裂缓解等效环流时间（ms 口径）', 2.0e6 / P_rec * 1e3, 'ms')
    p('AR.6', '直流母线 40 kA 短路时 24 kV 侧功率', 40e3 * 24e3, 'W', '= %.6g MW' % (40e3 * 24e3 / 1e6))

    print('==== [AR.7] 练习解答 ====')
    # 1) 电厂净效率对 eta_th 的敏感度
    p('AR.7', '练 1：eta_th=0.45 时 P_gross', 0.45 * P_th, 'W', '= %.6g MW' % (0.45 * P_th / 1e6))
    p('AR.7', '练 1：对应 f_rec', P_rec / (0.45 * P_th), '', '= %.6g %%' % (100 * P_rec / (0.45 * P_th)))
    p('AR.7', '练 1：对应 eta_net', (0.45 * P_th - P_rec) / P_th, '', '= %.6g %%' % (100 * (0.45 * P_th - P_rec) / P_th))
    p('AR.7', '练 1：净电相对提升倍数', (0.45 * P_th - P_rec) / P_net, '')
    p('AR.7', '练 1：eta_th=0.45 时的净电', 0.45 * P_th - P_rec, 'W', '= %.6g MW' % ((0.45 * P_th - P_rec) / 1e6))
    # 2) 飞轮 4000 rpm
    om2 = 2 * math.pi * 4000 / 60.0
    p('AR.7', '练 2：4000 rpm 飞轮角速度', om2, 'rad/s')
    p('AR.7', '练 2：4000 rpm 飞轮储能', 0.5 * I_fw * om2 ** 2, 'J', '= %.6g GJ' % (0.5 * I_fw * om2 ** 2 / 1e9))
    p('AR.7', '练 2：相对 3000 rpm 提升倍数', (om2 / omega) ** 2, '')
    p('AR.7', '练 2：4000 rpm 边缘线速度', om2 * r_fw, 'm/s')
    p('AR.7', '练 2：对应钢内应力', rho_s * (om2 * r_fw) ** 2, 'Pa', '= %.6g GPa' % (rho_s * (om2 * r_fw) ** 2 / 1e9))
    p('AR.7', '练 2：是否超过钢强度（1 GPa）', 'yes' if rho_s * (om2 * r_fw) ** 2 > 1.0e9 else 'no')
    # 3) 十二脉波改二十四脉波
    p('AR.7', '练 3：24 脉波纹波（文献级 2.1%）', 0.021, '', '= 2.1 %%')
    p('AR.7', '练 3：相对 12 脉波改善倍数', ripple12 / 0.021, '')
    p('AR.7', '练 3：最低谐波次数', 24 - 1, '')
    # 4) CS 爬升 20 s
    p('AR.7', '练 4：20 s 爬升平均功率', E_CS / 20.0, 'W', '= %.6g MW' % (E_CS / 20.0 / 1e6))
    p('AR.7', '练 4：dI/dt', 2 * I_pk / 20.0, 'A/s')
    p('AR.7', '练 4：端电压', L_CS * 2 * I_pk / 20.0, 'V', '= %.6g kV' % (L_CS * 2 * I_pk / 20.0 / 1e3))
    p('AR.7', '练 4：爬升功率/净电', (E_CS / 20.0) / P_net, '')
    # 5) 储能替代：用 8 GJ 飞轮撑 100 MW
    p('AR.7', '练 5：8 GJ 撑 100 MW 时间', 8e9 / 100e6, 's')
    p('AR.7', '练 5：相对 2 GJ 飞轮倍数', 8e9 / E_fw, '')
    p('AR.7', '练 5：碳纤维飞轮所需质量（E=8 GJ，v_max）', 2 * 8e9 / v_max_cf ** 2, 'kg', '= %.6g t' % (2 * 8e9 / v_max_cf ** 2 / 1e3))
    p('AR.7', '练 5：钢飞轮所需质量（v_max_s）', 2 * 8e9 / v_max_s ** 2, 'kg', '= %.6g t' % (2 * 8e9 / v_max_s ** 2 / 1e3))
    # 6) 无功补偿
    p('AR.7', '练 6：alpha=25 deg 时 tan', math.tan(math.radians(25.0)), '')
    p('AR.7', '练 6：对应无功', math.tan(math.radians(25.0)) * P_gross, 'var', '= %.6g MVAr' % (math.tan(math.radians(25.0)) * P_gross / 1e6))
    p('AR.7', '练 6：相对 alpha=15 deg 增加倍数', math.tan(math.radians(25.0)) / math.tan(ar), '')
    p('AR.7', '练 6：补偿到 0.98 所需容量', P_gross * math.tan(math.acos(0.98)), 'var', '= %.6g MVAr' % (P_gross * math.tan(math.acos(0.98)) / 1e6))

main()
