# -*- coding: utf-8 -*-
"""附录 AK《聚变安全与环境：从氚的克数到废物的分级》数值唯一来源。

只 print，不写任何文件。正文与门禁 P43 只能引用本脚本打印出来的数字。
生成 canon：python tmp\\appAK_calc.py | Out-File -Encoding utf8 tmp\\appAK_calc.txt
"""
import math

# ---- CODATA 2018 常数 ----
e = 1.602176634e-19          # C
u = 1.66053906660e-27        # kg
NA = 6.02214076e23           # 1/mol
MeV = 1e6 * e                # J
yr = 3.15576e7               # s（儒略年 365.25 d）
d = 86400.0

print('=' * 72)
print('附录 AK 数值唯一来源：tmp/appAK_calc.py（CODATA 2018）')
print('=' * 72)

# ---- DECLARED：工作点与工程输入（与附录 Z / AJ 同点）----
P_fus = 556.5e6              # W
P_n = 445.2e6                # W
P_alpha = 111.3e6            # W
q_n = 0.90944e6              # W/m^2（中子壁负载，附录 Z）
dpa_rate = 4.5472            # dpa/年（附录 Z）
life = 10.996                # 年（第一壁寿命，附录 Z）
A_wall = 489.532             # m^2（第一壁面积，附录 X 口径）
E_n = 14.1 * MeV             # J（D-T 中子能量）
T_half_T = 12.32 * yr        # s（氚半衰期）
m_T = 3.0160492 * u          # kg（氚原子质量）
m_burn = 31.191              # kg/年（氚燃耗，附录 Z/AJ）
m_flow = 162.72              # kg/年（加料氚流，附录 AJ）
retention = 0.01             # 氚滞留份额（附录 AJ）
inv_limit = 2.0              # kg（氚库存上限，附录 Z/AJ）
A_spec_steel = 1.0e6         # Bq/kg（寿期末活化钢比活度，文献量级·声明）
clearance = 1.0e4            # Bq/kg（钢的豁免/解控水平，文献量级·声明）
T_half_Co = 5.2711 * yr      # s（Co-60 半衰期）
decay_heat = 0.05            # W/kg（寿期末活化钢衰变热，文献量级·声明）
rho_steel = 7800.0           # kg/m^3
M_waste = 1.0e6              # kg（堆内构件待处置质量，声明）
j0 = 1.0e9                   # Bq/(m^2 s)（无屏障氚渗透通量，文献量级·声明）
PRF = 100.0                  # 渗透降低因子（带屏障，声明）
d_conc = 1.8e-11             # Sv/Bq（HTO 吸入剂量系数，ICRP 量级·声明）
d_conc_HT = 1.8e-15          # Sv/Bq（HT 气体，ICRP 量级·声明）
lim_pub = 1.0e-3             # Sv/年（公众剂量限值）
lim_occ = 20.0e-3            # Sv/年（职业剂量限值）
rel_lim_g = 1.0              # g/年（氚年释放管理目标，文献量级·声明）
D_per_g = 20.0e-6            # Sv（每释放 1 g 氚对应的公众剂量，文献量级·声明）
V_glove = 1000.0             # m^3（手套箱体积，声明）
A_glove0 = 1.0e6             # Bq/m^3（初始氚浓度，声明）
A_glove1 = 1.0e3             # Bq/m^3（目标浓度，声明）
k_clean = 1.0                # 1/h（除氚系统一级速率常数，声明）
leak_frac = 1.0e-3           # 事件期间安全壳泄漏份额（声明）
DF_detrit = 1.0e-2           # 排放前除氚剩余份额（声明）
rel_frac = 0.10              # 事故中堆内氚释放份额（声明）
print('[DECLARED] P_fus=%.4g MW  P_n=%.4g MW  q_n=%.4g MW/m^2  dpa=%.4g/年  寿命=%.4g 年'
      % (P_fus / 1e6, P_n / 1e6, q_n / 1e6, dpa_rate, life))
print('[DECLARED] A_wall=%.6g m^2  E_n=14.1 MeV  T_1/2(T)=12.32 年  m_T=3.0160492 u' % A_wall)
print('[DECLARED] 燃耗=%.5g kg/年  滞留=%.4g  库存上限=%.4g kg' % (m_burn, retention, inv_limit))
print('[DECLARED] 活化钢比活度=%.4g Bq/kg  豁免=%.4g Bq/kg  Co-60 T_1/2=5.2711 年  衰变热=%.4g W/kg'
      % (A_spec_steel, clearance, decay_heat))
print('[DECLARED] 渗透 j0=%.4g Bq/(m^2 s)  PRF=%.4g  HTO 剂量系数=%.4g Sv/Bq  HT=%.4g Sv/Bq'
      % (j0, PRF, d_conc, d_conc_HT))
print('[DECLARED] 公众限值=%.4g mSv/年  职业限值=%.4g mSv/年  年释放目标=%.4g g/年  每克释放剂量=%.4g uSv'
      % (lim_pub * 1e3, lim_occ * 1e3, rel_lim_g, D_per_g * 1e6))
print('[DECLARED] 手套箱=%.4g m^3  %.4g -> %.4g Bq/m^3  k=%.4g /h' % (V_glove, A_glove0, A_glove1, k_clean))
print('[DECLARED] 安全壳泄漏=%.4g  除氚后剩余=%.4g  事故释放份额=%.4g' % (leak_frac, DF_detrit, rel_frac))
print()

# ---- [AK.1] 第四本账：把安全写成数字 ----
print('[AK.1] 第四本账：把安全写成数字')
phi = q_n / E_n
print('  中子壁负载 %.6g W/m^2  ->  中子通量 phi = %.6g n/(m^2 s)' % (q_n, phi))
phi_y = phi * yr
print('  年注量 = %.6g n/(m^2 年)' % phi_y)
sigma_d = dpa_rate / phi_y
print('  位移截面标定 sigma_d = %.6g m^2 = %.6g b（与附录 Z 的 %.5g dpa/年 同口径）'
      % (sigma_d, sigma_d * 1e28, dpa_rate))
dpa_life = dpa_rate * life
print('  寿期累计位移损伤 = %.6g dpa' % dpa_life)
P_n_frac = P_n / P_fus
print('  中子功率份额 = %.6g（%.6g%% of P_fus）' % (P_n_frac, P_n_frac * 100))
print('  三条账（Z/AJ）之外第四本账的两个自变量：氚活度与活化比活度')
print()

# ---- [AK.2] 氚的活度与剂量 ----
print('[AK.2] 氚的活度与剂量：克数怎么变成西弗')
lam_T = math.log(2) / T_half_T
print('  衰变常数 lambda = %.6g /s' % lam_T)
N_1g = 1e-3 / m_T
print('  1 g 氚的原子数 = %.6g' % N_1g)
A_1g = lam_T * N_1g
print('  1 g 氚活度 = %.6g Bq = %.6g TBq（比活度）' % (A_1g, A_1g / 1e12))
print('  0.1 g -> %.6g Bq；1 mg -> %.6g Bq' % (A_1g * 0.1, A_1g * 1e-3))
A_inv = A_1g * inv_limit * 1e3
print('  库存上限 %.4g kg 对应活度 = %.6g Bq' % (inv_limit, A_inv))
D_1g = A_1g * d_conc
print('  1 g 全部氧化吸入的待积有效剂量 = %.6g Sv' % D_1g)
A_1mSv = lim_pub / d_conc
print('  公众限值 %.4g mSv 对应摄入 = %.6g Bq = %.6g g = %.6g ug' % (lim_pub * 1e3, A_1mSv, A_1mSv / A_1g, A_1mSv / A_1g * 1e6))
A_occ = lim_occ / d_conc
print('  职业限值 %.4g mSv 对应摄入 = %.6g Bq = %.6g g' % (lim_occ * 1e3, A_occ, A_occ / A_1g))
print('  HT 气体同活度的剂量 = %.6g Sv（是 HTO 的 %.6g 分之一）'
      % (A_1g * d_conc_HT, d_conc / d_conc_HT))
A_burn = A_1g * m_burn * 1e3
print('  年燃耗 %.5g kg = %.6g Bq/年' % (m_burn, A_burn))
D_1pct = A_burn * 0.01 * d_conc
print('  年燃耗的 1%% 全部氧化吸入 = %.6g Sv（一个人吸入）' % D_1pct)
print('  剂量阶梯（表 AK.2 口径）：1 mg -> %.6g Sv；0.1 g -> %.6g Sv；'
      '2 kg -> %.6g Sv；年燃耗 100%% -> %.6g Sv'
      % (A_1g * 1e-3 * d_conc, A_1g * 0.1 * d_conc, A_inv * d_conc, A_burn * d_conc))
print('  练习 4 输入：库存 1 kg 的活度 = %.6g Bq' % (A_1g * 1e3))
print()

# ---- [AK.3] 氚的包容：渗透、除氚与释放限值 ----
print('[AK.3] 氚的包容：渗透、除氚与释放限值')
j = j0 / PRF
print('  带屏障渗透通量 j = j0/PRF = %.6g Bq/(m^2 s)' % j)
perm_rate = j * A_wall
print('  全壁渗透率 = %.6g Bq/s = %.6g Bq/天 = %.6g TBq/年'
      % (perm_rate, perm_rate * d, perm_rate * yr / 1e12))
perm_g_day = perm_rate * d / A_1g
perm_kg_y = perm_rate * yr / A_1g / 1e3
print('  折合质量 = %.6g g/天 = %.6g kg/年' % (perm_g_day, perm_kg_y))
print('  与燃耗 %.5g kg/年 之比 = %.6g' % (m_burn, perm_kg_y / m_burn))
A_release = A_1g * rel_lim_g
print('  年释放管理目标 %.4g g/年 = %.6g Bq/年' % (rel_lim_g, A_release))
print('  渗透量 / 管理目标 = %.6g 倍（必须靠除氚系统补上）' % (perm_rate * yr / A_release))
DF_need = perm_rate * yr / A_release
print('  需要的总除氚因子 DF = %.6g' % DF_need)
print('  另一条线（附录 AJ 口径）：加料氚流 %.6g kg/年 的 %.4g%% 滞留 = %.6g kg/年；'
      '库存上限 %.4g kg -> %.6g 年' % (m_flow, retention * 100, m_flow * retention,
                                    inv_limit, inv_limit / (m_flow * retention)))
t_clean = math.log(A_glove0 / A_glove1) / k_clean
print('  手套箱 %.4g m^3 从 %.4g 降到 %.4g Bq/m^3：%.6g 小时' % (V_glove, A_glove0, A_glove1, t_clean))
print('  一级去除：浓度比 %.6g 需要 ln = %.6g 个时间常数' % (A_glove0 / A_glove1, math.log(A_glove0 / A_glove1)))
print()

# ---- [AK.4] 中子活化 ----
print('[AK.4] 中子活化：壁材料怎么变成放射性废物')
print('  中子通量 phi = %.6g n/(m^2 s)、年注量 = %.6g n/(m^2 年)' % (phi, phi_y))
print('  位移截面 sigma_d = %.6g b（由附录 Z 的 %.5g dpa/年 标定）' % (sigma_d * 1e28, dpa_rate))
print('  寿期（%.4g 年）累计 = %.6g dpa' % (life, dpa_life))
He_rate = 10.0            # appm/dpa（声明：Fe 的 He 产额量级）
print('  若 He/dpa = %.4g appm/dpa，则寿期累计 He = %.6g appm（声明值）' % (He_rate, He_rate * dpa_life))
V_waste = M_waste / rho_steel
print('  待处置堆内构件 %.4g kg -> 体积 %.6g m^3' % (M_waste, V_waste))
A_waste = A_spec_steel * M_waste
print('  比活度 %.4g Bq/kg -> 总活度 %.6g Bq = %.6g TBq' % (A_spec_steel, A_waste, A_waste / 1e12))
print('  是豁免水平 %.4g Bq/kg 的 %.6g 倍' % (clearance, A_spec_steel / clearance))
print('  单位质量动能对照：1 g 氚 = %.6g Bq；1 kg 活化钢 = %.6g Bq' % (A_1g, A_spec_steel))
print()

# ---- [AK.5] 衰变热与废物分级 ----
print('[AK.5] 衰变热与废物分级：冷却多久才能豁免')
lam_Co = math.log(2) / T_half_Co
print('  Co-60 衰变常数 = %.6g /s' % lam_Co)
t_clear = math.log(A_spec_steel / clearance) / lam_Co
print('  从 %.4g 降到 %.4g Bq/kg 需要 %.6g 年' % (A_spec_steel, clearance, t_clear / yr))
A_50 = A_spec_steel * math.exp(-lam_Co * 50 * yr)
print('  冷却 50 年后的比活度 = %.6g Bq/kg（豁免水平的 %.6g 倍）' % (A_50, A_50 / clearance))
A_100 = A_spec_steel * math.exp(-lam_Co * 100 * yr)
print('  冷却 100 年后的比活度 = %.6g Bq/kg' % A_100)
P_decay0 = decay_heat * M_waste
print('  停机瞬时衰变热 = %.4g W/kg x %.4g kg = %.6g W = %.6g kW'
      % (decay_heat, M_waste, P_decay0, P_decay0 / 1e3))
P_50 = P_decay0 * math.exp(-lam_Co * 50 * yr)
print('  50 年后衰变热 = %.6g W（衰减 %.6g 倍）' % (P_50, P_decay0 / P_50))
print('  废物分级只看两个量：比活度与衰变热；两者都由 Co-60 类核素主导（声明）')
print()

# ---- [AK.6] 事故与纵深防御 ----
print('[AK.6] 事故与纵深防御：LOCA 情景的剂量账')
f_site = D_per_g / A_1g
print('  场址剂量因子 f = %.4g uSv / %.4g g = %.6g Sv/Bq' % (D_per_g * 1e6, rel_lim_g, f_site))
A_core = A_inv
print('  堆内氚库存 = %.6g Bq（%.4g kg = %.4g g）' % (A_core, inv_limit, inv_limit * 1e3))
A_rel = A_core * rel_frac
print('  事故释放份额 %.4g -> 进安全壳 %.6g Bq = %.6g g' % (rel_frac, A_rel, A_rel / A_1g))
D_no = A_rel * f_site
print('  无包容（直接排到边界）剂量 = %.6g Sv = %.6g mSv' % (D_no, D_no * 1e3))
A_env = A_rel * leak_frac * DF_detrit
print('  有包容（事件期泄漏 %.4g、除氚后剩余 %.4g）-> 环境释放 %.6g Bq = %.6g g'
      % (leak_frac, DF_detrit, A_env, A_env / A_1g))
D_yes = A_env * f_site
print('  有包容剂量 = %.6g Sv = %.6g uSv' % (D_yes, D_yes * 1e6))
print('  包容带来的剂量降低倍数 = %.6g' % (D_no / D_yes))
print('  与公众限值 %.4g mSv 之比 = %.6g（有包容）' % (lim_pub * 1e3, D_yes / lim_pub))
print('  事故余量 = %.6g 倍' % (lim_pub / D_yes))
print('  三个应回答的问题：库存多少、释放多少、到人多少（乘法定律）')
print()

# ---- [AK.7] 六道练习解答 ----
print('[AK.7] 六道练习解答')
A_10mg = A_1g * 0.01
print('  1) 10 mg 氚：活度 = %.6g Bq = %.6g GBq；全部吸入剂量 = %.6g Sv'
      % (A_10mg, A_10mg / 1e9, A_10mg * d_conc))
A_HT = lim_pub / d_conc_HT
print('  2) 用 HT 剂量系数：限值对应 %.6g Bq = %.6g g（比 HTO 宽 %.6g 倍）'
      % (A_HT, A_HT / A_1g, d_conc / d_conc_HT))
j2 = 1.0e9 / 1000.0
perm2 = j2 * A_wall * d / A_1g
print('  3) PRF=1000：j=%.6g Bq/(m^2 s)，%.6g g/天 = %.6g kg/年' % (j2, perm2, perm2 * 365.25 / 1e3))
A_core2 = A_1g * 1e3
A_env2 = A_core2 * 0.10 * leak_frac * DF_detrit
print('  4) 库存 1 kg、释放 10%%、同样包容：环境释放 %.6g Bq，剂量 %.6g uSv'
      % (A_env2, A_env2 * f_site * 1e6))
t_clear2 = math.log(A_spec_steel / (clearance * 10)) / lam_Co
print('  5) 比活度降到 10 倍豁免（%.4g Bq/kg）需要 %.6g 年' % (clearance * 10, t_clear2 / yr))
He2 = 5.0 * dpa_life
print('  6) He/dpa=%.4g appm/dpa 时寿期累计 He = %.6g appm（%.4g 声明值的一半）' % (5.0, He2, He_rate))
print()
print('以上为可以引用的全部数值（正文与门禁 P43 只引用本输出）。')
