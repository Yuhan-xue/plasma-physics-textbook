# -*- coding: utf-8 -*-
"""附录 AP「维护与遥操作」数值唯一来源（只 print，不写文件）。

用法：python tmp\\appAP_calc.py | Out-File -Encoding utf8 tmp\\appAP_calc.txt

DECLARED（声明输入，非本书推导）：
  工作点与 Z/AJ/AK/AO 同：P_fus = 556.5 MW
  电厂净电功率 P_e = 200 MW、电价 0.4 元/kWh、年 8760 h
  计划维护 60 d/年、非计划 20 d/年
  包层总质量 1e6 kg（与 AK 的构件质量同口径）、扇形段 16 段
  包层更换周期 5 年、电厂寿命 40 年
  停堆 1 天后剂量率（短寿命分量）1e4 Sv/h、T_half = 2.5 d
  长寿命分量（Co-60 主导）10 mSv/h、T_half = 5.2711 年（与 AK 同）
  远程操作阈值 1 Sv/h；人工进入阈值 10 uSv/h
  现场屏蔽因子 1e-4；器件耐总剂量 1e5 Gy（1 Gy ~ 1 Sv 近似）
  机械臂臂长 5 m、定位精度 0.2 mm/m、力反馈 0.1 N、主从延迟 100 ms
  切割/焊接装配间隙 0.5 mm、插拔力 100 N
  热室工位 2、单段周转 30 d、cask 2 个（一次 1 段）
  人工剂量限值 20 mSv/年、单次检修 0.5 mSv、在役检查 100 次/年、远程比例 95%
"""
from __future__ import print_function

import math

LN2 = math.log(2.0)

# ---- DECLARED ----
P_FUS = 556.5e6          # W
P_E = 200e6              # W
PRICE = 0.4              # 元/kWh
H_YEAR = 8760.0          # h
D_PLAN = 60.0            # d
D_UNPLAN = 20.0          # d
M_BLANKET = 1.0e6        # kg
N_SECTOR = 16
T_BLANKET = 5.0          # 年
T_PLANT = 40.0           # 年
D0_SHORT = 1.0e4         # Sv/h（停堆 1 天后）
T_SHORT = 2.5            # d
D_LONG = 1.0e-2          # Sv/h = 10 mSv/h
T_LONG = 5.2711          # 年
D_REMOTE = 1.0            # Sv/h
D_HANDS = 1.0e-5          # Sv/h = 10 uSv/h
F_SHIELD = 1.0e-4
DOSE_LIMIT_DEV = 1.0e5   # Gy
ARM_LEN = 5.0            # m
ACC_PER_M = 0.2e-3       # m/m
F_FB = 0.1               # N
LATENCY = 0.1            # s
GAP = 0.5e-3             # m
F_PLUG = 100.0           # N
N_HOTCELL = 2
T_TURN = 30.0            # d
N_CASK = 2
DOSE_LIMIT_HUMAN = 20e-3  # Sv/年
DOSE_INSPECT = 0.5e-3     # Sv/次
N_INSPECT = 100          # 次/年
F_REMOTE = 0.95


def p(tag, name, value, unit=''):
    print('[%s] %-44s %s %s' % (tag, name, ('%.6g' % value), unit))


print('=== [AP.1] 可用因子与停机代价 ===')
down_h = (D_PLAN + D_UNPLAN) * 24.0
avail = 1.0 - down_h / H_YEAR
p('AP.1', '每年停机小时', down_h, 'h')
p('AP.1', '可用因子', avail)
lost_mwh = P_E / 1e6 * down_h          # MWh
p('AP.1', '停机的发电量损失', lost_mwh, 'MWh')
p('AP.1', '停机的售电损失', lost_mwh * 1e3 * PRICE, '元')
year_mwh = P_E / 1e6 * H_YEAR
p('AP.1', '全年满发发电量', year_mwh, 'MWh')
p('AP.1', '相对全年发电量的损失份额', 1.0 - avail)
p('AP.1', '每停机一天的损失', P_E / 1e6 * 24 * 1e3 * PRICE, '元/天')

print('')
print('=== [AP.2] 包层更换周期 ===')
m_sec = M_BLANKET / N_SECTOR
p('AP.2', '单个扇形段质量', m_sec, 'kg')
p('AP.2', '单个扇形段质量', m_sec / 1e3, 't')
n_change = T_PLANT / T_BLANKET
p('AP.2', '全寿期更换次数', n_change, '次')
p('AP.2', '全寿期更换段数', n_change * N_SECTOR, '段')
p('AP.2', '平均每年更换段数', n_change * N_SECTOR / T_PLANT, '段/年')
p('AP.2', '平均每年更换质量', n_change * N_SECTOR / T_PLANT * m_sec / 1e3, 't/年')
p('AP.2', '每段平均可用工期（按每年非计划+计划停机）', down_h / (n_change * N_SECTOR / T_PLANT), 'h/段')

print('')
print('=== [AP.3] 剂量率衰减与冷却窗口 ===')
t_remote = T_SHORT * math.log(D0_SHORT / (D_REMOTE - D_LONG), 2.0)
p('AP.3', '降到远程操作阈值所需冷却时间', t_remote, 'd')
t_01 = T_SHORT * math.log(D0_SHORT / (0.1 - D_LONG), 2.0)
p('AP.3', '降到 0.1 Sv/h 所需冷却时间', t_01, 'd')
p('AP.3', '一年后的剂量率', D_LONG + D0_SHORT * 2 ** (-365.0 / T_SHORT), 'Sv/h')
p('AP.3', '十年后的剂量率', D_LONG * 2 ** (-10.0 / T_LONG), 'Sv/h')
t_hands = T_LONG * math.log(D_LONG / D_HANDS, 2.0)
p('AP.3', '降到人工阈值所需冷却时间', t_hands, '年')
p('AP.3', '人工阈值相当于', D_HANDS * 1e6, 'uSv/h')
# 从远程阈值起、工期 30 天内的累积剂量（未屏蔽）
t0, t1 = t_remote, t_remote + T_TURN
short_int = D0_SHORT * (T_SHORT / LN2) * (2 ** (-t0 / T_SHORT) - 2 ** (-t1 / T_SHORT))
acc = D_LONG * (T_TURN * 24.0) + short_int
p('AP.3', '30 天工期内的累积剂量（未屏蔽）', acc, 'Sv')
p('AP.3', '屏蔽后累积剂量', acc * F_SHIELD, 'Sv')
p('AP.3', '屏蔽后相对器件耐总剂量的份额', acc * F_SHIELD / DOSE_LIMIT_DEV)
p('AP.3', '屏蔽后累积到器件上限所需工期', DOSE_LIMIT_DEV / (acc * F_SHIELD / (T_TURN * 24.0)), 'h')

print('')
print('=== [AP.4] 遥操作精度与力反馈 ===')
err = ARM_LEN * ACC_PER_M
p('AP.4', '臂端定位误差', err, 'm')
p('AP.4', '相对装配间隙的余量', GAP / err)
p('AP.4', '力反馈相对插拔力的分辨率', F_FB / F_PLUG)
p('AP.4', '主从延迟对应的最大操作速度（按 1 cm 目标）', 0.01 / LATENCY, 'm/s')
p('AP.4', '带宽上限对应的周期', 1.0 / 10.0, 's')

print('')
print('=== [AP.5] 热室与 cask 流转 ===')
need = n_change * N_SECTOR / T_PLANT
cap = N_CASK * 365.0 / T_TURN
p('AP.5', '年更换需求', need, '段/年')
p('AP.5', 'cask 年转运能力', cap, '段/年')
p('AP.5', '能力余量', cap / need)
p('AP.5', '单工位年周转能力', 365.0 / T_TURN, '段/年')
p('AP.5', '热室同时占用的工位需求', need * T_TURN / 365.0)
p('AP.5', '每次转运的屏蔽质量（按 1:1 屏蔽比，声明）', m_sec, 'kg')

print('')
print('=== [AP.6] 人因与在役检查 ===')
hands_ok = DOSE_LIMIT_HUMAN / DOSE_INSPECT
p('AP.6', '按单次剂量允许的人工次数', hands_ok, '次/年')
n_hands = N_INSPECT * (1.0 - F_REMOTE)
p('AP.6', '在役检查中的人工次数', n_hands, '次/年')
p('AP.6', '人工年剂量', n_hands * DOSE_INSPECT * 1e3, 'mSv/年')
p('AP.6', '占限值份额', n_hands * DOSE_INSPECT / DOSE_LIMIT_HUMAN)
p('AP.6', '远程检修相对人工的次数倍数', F_REMOTE / (1.0 - F_REMOTE))

print('')
print('=== [AP.7] 练习解答 ===')
# 1) 计划维护 45 d、非计划 20 d -> 可用因子
a1 = 1.0 - (45 + 20) * 24.0 / H_YEAR
p('AP.7-1', '维护 45+20 d 的可用因子', a1)
p('AP.7-1', '相对 60+20 d 方案提高的份额', a1 - avail)
# 2) 包层 3 年一换 -> 段/年
p('AP.7-2', '包层 3 年一换的年更换段数', T_PLANT / 3.0 * N_SECTOR / T_PLANT, '段/年')
p('AP.7-2', '相对 5 年方案的倍数', (T_PLANT / 3.0) / n_change)
# 3) 短寿命 T_half = 1 d -> 到 1 Sv/h 时间
p('AP.7-3', 'T_half=1 d 时的冷却时间', 1.0 * math.log(D0_SHORT / (D_REMOTE - D_LONG), 2.0), 'd')
# 4) 屏蔽因子 1e-5 -> 累积
p('AP.7-4', '屏蔽因子 1e-5 时 30 天累积剂量', acc * 1e-5, 'Sv')
p('AP.7-4', '相对器件上限份额', acc * 1e-5 / DOSE_LIMIT_DEV)
# 5) 限值 10 mSv -> 人工次数
p('AP.7-5', '限值 10 mSv 时允许人工次数', 10e-3 / DOSE_INSPECT, '次/年')
# 6) cask 4 个 -> 年能力
p('AP.7-6', 'cask 4 个的年转运能力', 4 * 365.0 / T_TURN, '段/年')
p('AP.7-6', '相对需求的余量', 4 * 365.0 / T_TURN / need)

print('')
print('=== [AP.8] 派生比值（供交叉核对） ===')
p('AP.8', '更换一段的相对成本（停机小时/段）', down_h / need, 'h/段')
p('AP.8', 'P_fus / P_e（声明值的比值）', P_FUS / P_E)
p('AP.8', '每年停机损失与 P_e 的比值（按满发小时）', (1.0 - avail))
p('AP.8', '长寿命分量占停堆 1 天后剂量率的份额', D_LONG / (D_LONG + D0_SHORT))
p('AP.8', '远程阈值相对长寿命分量的倍数', D_REMOTE / D_LONG)
p('AP.8', '人工阈值相对远程阈值的倍数', D_REMOTE / D_HANDS)


# ---- [AP.8b] 派生比值 ----
p('AP.8b', '热室能力余量百分数', (cap / need) * 100.0, '%')
p('AP.8b', '人员剂量余量百分数',
  (1.0 - n_hands * DOSE_INSPECT / DOSE_LIMIT_HUMAN) * 100.0, '%')


# ---- [AP.8c] 百分比与练习派生量 ----
p('AP.8c', '占限值余量百分数',
  100.0 * (1.0 - n_hands * DOSE_INSPECT / DOSE_LIMIT_HUMAN), '%')
p('AP.8c', '练习 1 压缩的计划维护天数', 60.0 - 45.0, 'd')
p('AP.8c', '练习 1 可用因子提高的百分点',
  100.0 * (a1 - (1.0 - (60 + 20) * 24.0 / H_YEAR)), '%')
p('AP.8c', '练习 2 三年一换的年更换质量',
  (M_BLANKET / 1000.0) / 3.0, 't')
