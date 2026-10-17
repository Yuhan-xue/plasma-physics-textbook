# -*- coding: utf-8 -*-
"""附录 AV《氚循环与氚工厂：从增殖到衡算》数值唯一来源（只 print，不写文件）。

口径：
  - 工作点沿用全书 D1：P_fus = 556.5 MW、P_alpha = 111.3 MW、V = 832.2 m^3、
    R = 6.2 m、a = 2.0 m、n_e = 1e20 m^-3、T_e = 9.194 keV；第一壁面积取附录 X 的 489.532 m^2；
  - 核数据：D-T 总释能 17.586 MeV（中子 14.068、alpha 3.518）；
    6Li(n,alpha)T 释能 4.78 MeV；7Li(n,n'alpha)T 阈值 2.47 MeV；9Be(n,2n) 阈值 1.9 MeV；
  - 氚半衰期 12.32 年、平均 beta 能量 5.7 keV、比活度按 T 原子 3.0160492 u 计；
  - 声明输入（设计取值，门禁白名单）：包层覆盖率/厚度/密度、屏蔽衰减长度、
    渗透率有效值 Phi_eff、提取效率、燃烧率、分离因子、库存与衡算精度、氚水排放限值。

运行：python tmp\\appAV_calc.py > tmp\\appAV_calc.txt（PowerShell 用 Out-File -Encoding utf8）
"""
from __future__ import print_function

import math

# ---------------- 常数 ----------------
E_CHG = 1.602176634e-19      # C
U_KG = 1.66053906660e-27     # kg
N_A = 6.02214076e23          # 1/mol
SEC_YEAR = 3.15576e7         # s/a（儒略年）
SEC_DAY = 86400.0

# ---------------- 质量（u） ----------------
M_D = 2.01410177812
M_T = 3.01604928199
M_LI6 = 6.0151228874
M_LI7 = 7.0160034366
M_BE9 = 9.012183065
M_HE4 = 4.00260325413

# ---------------- 声明输入 ----------------
P_FUS = 556.5e6              # W
P_ALPHA = 111.3e6            # W
VOL = 832.2                  # m^3
A_WALL = 489.532             # m^2（4 pi^2 R a，与附录 X 一致）
NE0 = 1.0e20                 # m^-3
E_FUS_MEV = 17.586
E_N_MEV = 14.068
E_ALPHA_MEV = 3.518
T_HALF_Y = 12.32             # a
E_BETA_KEV = 5.7             # 平均 beta 能量（keV）
LI6_ABUND = 0.0759           # 天然锂中 6Li 丰度
LI6_ENRICH = 0.60            # 包层锂富集度（声明值）
LI7_ENRICH = 1.0 - LI6_ENRICH
TBR_DESIGN = 1.05            # 设计氚增殖比（声明值）
COVERAGE = 0.80              # 包层覆盖率（声明值）
BLANKET_THICK = 0.5          # m（声明值）
PACKING = 0.60               # 陶瓷堆积因子（声明值）
RHO_LI4SIO4 = 2400.0         # kg/m^3（声明值）
SHIELD_LAMBDA = 0.075        # m（钢/水混合有效衰减长度，声明值）
SHIELD_REDUCTION = 1.0e6     # 需要降低的倍数（声明值）
PEAK_FACTOR = 1.3            # 中子壁负载峰值因子（声明值）
PHI_EFF = 1.0e-10            # mol m^-1 s^-1 Pa^-0.5（渗透率有效值，声明值）
PERM_THICK = 2.0e-3          # m（声明值）
P_T = 10.0                   # Pa（氚分压，声明值）
A_HOT = 500.0                # m^2（热表面，声明值）
PRF = 100.0                  # 渗透降低因子（声明值）
ETA_EXT = 0.90               # 提取效率（声明值）
BURN_FRAC = 0.02             # 燃烧率（声明值）
ALPHA_DT = 1.5               # 低温蒸馏 D/T 分离因子（声明值）
ALPHA_DT2 = 1.3              # 对照分离因子（声明值）
X_FEED = 0.50                # 蒸馏进料氚份额（声明值）
X_PROD = 0.90                # 产品氚份额（声明值）
INV_KG = 10.0                # 氚库存（声明值）
ACC_ERR = 0.01               # 衡算不确定度（声明值）
ACC_ERR2 = 0.005             # 目标衡算不确定度（声明值）
HTO_LIMIT = 1.0e4            # Bq/L（饮用水氚限值，声明值）
CO_DEP_AJ = 603.647          # g（附录 AJ 共沉积氚）

print('==== 声明输入 ====')
print('P_fus = %.1f MW, P_alpha = %.1f MW, V = %.1f m^3, A_wall = %.3f m^2'
      % (P_FUS / 1e6, P_ALPHA / 1e6, VOL, A_WALL))
print('E_fus = %.3f MeV (n %.3f + alpha %.3f), T_half = %.2f a, E_beta = %.1f keV'
      % (E_FUS_MEV, E_N_MEV, E_ALPHA_MEV, T_HALF_Y, E_BETA_KEV))
print('6Li 丰度 %.4g（天然）/ 富集度 %.4g；TBR 设计 %.4g' % (LI6_ABUND, LI6_ENRICH, TBR_DESIGN))
print('包层：覆盖 %.4g、厚 %.4g m、堆积 %.4g、密度 %.4g kg/m^3；屏蔽 lambda = %.4g m'
      % (COVERAGE, BLANKET_THICK, PACKING, RHO_LI4SIO4, SHIELD_LAMBDA))
print('渗透：Phi_eff = %.6g mol/(m s Pa^0.5)、d = %.4g m、p_T = %.4g Pa、A_hot = %.4g m^2、PRF = %.4g'
      % (PHI_EFF, PERM_THICK, P_T, A_HOT, PRF))
print('提取效率 %.4g、燃烧率 %.4g、分离因子 %.4g / %.4g、库存 %.4g kg、衡算 %.4g / %.4g'
      % (ETA_EXT, BURN_FRAC, ALPHA_DT, ALPHA_DT2, INV_KG, ACC_ERR, ACC_ERR2))
print()

# ================= AV.1 三本账：反应率、消耗与比耗 =================
print('==== AV.1 三本账：一个反应、两份燃料、三笔账 ====')
E_FUS_J = E_FUS_MEV * 1e6 * E_CHG
E_N_J = E_N_MEV * 1e6 * E_CHG
E_ALPHA_J = E_ALPHA_MEV * 1e6 * E_CHG
R_N = P_FUS / E_FUS_J                      # 中子/氚产生率（1/s）
frac_n = E_N_MEV / E_FUS_MEV
frac_alpha = E_ALPHA_MEV / E_FUS_MEV
print('每次 D-T 反应释能 E_fus = %.6g J -> 反应率 R = P_fus/E_fus = %.6g s^-1'
      % (E_FUS_J, R_N))
print('  中子功率份额 = %.5f%%（中子 %.6g J）；alpha 功率份额 = %.5f%%'
      % (100 * frac_n, E_N_J, 100 * frac_alpha))
print('  中子功率 P_n = %.6g MW；alpha 功率 = P_fus - P_n = %.6g MW'
      % (P_FUS * frac_n / 1e6, (P_FUS - P_FUS * frac_n) / 1e6))
m_T_g = M_T * U_KG * 1e3                  # 每个氚原子质量（g）
m_D_g = M_D * U_KG * 1e3
t_burn_g_s = R_N * m_T_g                  # g/s
print('氚原子质量 = %.6g g；燃烧掉的氚 = %.6g g/s' % (m_T_g, t_burn_g_s))
print('  -> %.6g g/day；%.6g kg/a' % (t_burn_g_s * SEC_DAY, t_burn_g_s * SEC_YEAR / 1e3))
print('  -> 同法氘消耗 %.6g g/day；每 MW 聚变功率耗氚 %.6g g/(MW day)'
      % (R_N * m_D_g * SEC_DAY, t_burn_g_s * SEC_DAY / (P_FUS / 1e6)))
# 反推工作点所需的 <sigma v>（自洽性检查）
n_d = 0.5 * NE0
sv_need = R_N / (n_d * n_d * VOL)
print('自洽性检查：由 R = n_D n_T <sigma v> V 反推 <sigma v> = %.6g m^3/s（n_D = n_T = %.6g m^-3）'
      % (sv_need, n_d))
print('  每秒钟发生的反应数 = %.6g；一年 %.6g 次' % (R_N, R_N * SEC_YEAR))
print()

# ================= AV.2 增殖：TBR 与锂 =================
print('==== AV.2 增殖：TBR 是「每烧一个氚要生回几个」 ====')
t_bred_g_s = TBR_DESIGN * R_N * m_T_g
print('TBR 设计 %.4g -> 产生率 = TBR x R = %.6g 氚/s = %.6g g/s = %.6g g/day（燃烧 %.6g g/day）'
      % (TBR_DESIGN, TBR_DESIGN * R_N, t_bred_g_s, t_bred_g_s * SEC_DAY, t_burn_g_s * SEC_DAY))
print('6Li(n,alpha)T 放能 %.4g MeV（热中子截面极大）；7Li(n,n\'alpha)T 阈值 %.4g MeV；9Be(n,2n) 阈值 %.4g MeV'
      % (4.78, 2.47, 1.9))
li6_need_g_s = TBR_DESIGN * R_N * M_LI6 * U_KG * 1e3
print('6Li 消耗 = %.6g g/s = %.6g kg/a（每个氚耗一个 6Li）'
      % (li6_need_g_s, li6_need_g_s * SEC_YEAR / 1e3))
li_need = li6_need_g_s / LI6_ENRICH
li_nat = li6_need_g_s / LI6_ABUND
print('富集度 %.4g 的锂消耗 = %.6g kg/a；天然锂 = %.6g kg/a'
      % (LI6_ENRICH, li_need * SEC_YEAR / 1e3, li_nat * SEC_YEAR / 1e3))
m_li4sio4 = 4 * M_LI6 + 28.0855 + 4 * 15.999   # Li4SiO4 摩尔质量（g/mol）
print('Li4SiO4 摩尔质量 = %.6g g/mol；按富集度 %.4g 折算陶瓷消耗 = %.6g kg/a'
      % (m_li4sio4, LI6_ENRICH, (li_need / 4.0) * m_li4sio4 * SEC_YEAR / 1e3))
# 7Li 与 Be 的倍增作用
print('7Li(n,n\'alpha)T 阈值 %.4g MeV、9Be(n,2n) 阈值 %.4g MeV：前者不耗中子，后者把 1 个中子变 2 个'
      % (2.47, 1.9))
print('若全靠 6Li：TBR 最多 = 1（每中子一个氚）；要靠 7Li 与 Be 才能把 TBR 抬到 %.4g' % TBR_DESIGN)
print()

# ================= AV.3 包层与中子壁负载 =================
print('==== AV.3 包层与中子壁负载：0.9 MW/m^2 是什么概念 ====')
P_N = P_FUS * frac_n
nwl = P_N / A_WALL
nwl_peak = nwl * PEAK_FACTOR
print('P_n = %.6g MW；A_wall = %.3f m^2 -> 平均中子壁负载 = %.6g W/m^2 -> %.6f MW/m^2'
      % (P_N / 1e6, A_WALL, nwl, nwl / 1e6))
print('  峰值因子 %.4g -> 峰值 %.6g W/m^2 -> %.6f MW/m^2' % (PEAK_FACTOR, nwl_peak, nwl_peak / 1e6))
vol_blanket = A_WALL * COVERAGE * BLANKET_THICK
m_blanket = vol_blanket * PACKING * RHO_LI4SIO4
print('包层体积 = A_wall x 覆盖 x 厚度 = %.6g m^3；Li4SiO4 装量 = %.6g kg = %.6g t'
      % (vol_blanket, m_blanket, m_blanket / 1e3))
li6_in_blanket = m_blanket * 1e3 / m_li4sio4 * 4 * LI6_ENRICH * M_LI6
print('  其中 6Li = %.6g kg（按 %.4g 富集）' % (li6_in_blanket / 1e3, LI6_ENRICH))
print('  按 %.6g kg/a 的 6Li 消耗，装量可支撑 %.6g a（不计在线补充）'
      % (li6_need_g_s * SEC_YEAR / 1e3, li6_in_blanket / 1e3 / (li6_need_g_s * SEC_YEAR / 1e3)))
print('屏蔽：要把通量降 %.0e 倍需 %.6g m（lambda = %.4g m，e^-d/lambda；ln(1e6) = %.6g）'
      % (SHIELD_REDUCTION, SHIELD_LAMBDA * math.log(SHIELD_REDUCTION), SHIELD_LAMBDA,
         math.log(SHIELD_REDUCTION)))
print('  每米钢水混合屏蔽的衰减倍数 = e^(1/lambda) = %.6g' % (math.exp(1.0 / SHIELD_LAMBDA)))
print()

# ================= AV.4 提取与同位素分离 =================
print('==== AV.4 提取与同位素分离：燃烧率决定了工厂多大 ====')
t_inject = t_burn_g_s / BURN_FRAC
print('燃烧率 %.4g -> 注入氚 = 燃烧/燃烧率 = %.6g g/s = %.6g kg/day'
      % (BURN_FRAC, t_inject, t_inject * SEC_DAY / 1e3))
print('  年通过量 = %.6g t/a；未燃份额 %.4g 必须回收；通过量是燃烧量的 %.6g 倍'
      % (t_inject * SEC_YEAR / 1e6, 1 - BURN_FRAC, 1.0 / BURN_FRAC))
tbr_min_ext = 1.0 / ETA_EXT
print('提取效率 %.4g -> 若全部靠增殖补回，TBR 至少要 1/eta = %.6g' % (ETA_EXT, tbr_min_ext))
# 低温蒸馏理论板数
def n_stage(x_feed, x_prod, alpha):
    return math.log((x_prod / (1 - x_prod)) / (x_feed / (1 - x_feed))) / math.log(alpha)
n1 = n_stage(X_FEED, X_PROD, ALPHA_DT)
n2 = n_stage(X_FEED, X_PROD, ALPHA_DT2)
print('低温蒸馏：进料 x = %.4g -> 产品 x = %.4g，alpha = %.4g -> 理论板数 %.6g'
      % (X_FEED, X_PROD, ALPHA_DT, n1))
print('  若 alpha = %.4g -> %.6g 板（分离因子差 %.4g 倍，板数差 %.6g 倍）'
      % (ALPHA_DT2, n2, ALPHA_DT / ALPHA_DT2, n2 / n1))
x99 = 0.99
print('  要做到 x = %.4g -> %.6g 板（alpha = %.4g）' % (x99, n_stage(X_FEED, x99, ALPHA_DT), ALPHA_DT))
# 吹扫气与渗透窗面积
n_T_gas = P_T * 1.0 / (8.314 * 500.0)      # mol/m^3（理想气体 500 K）
print('吹扫气：p_T = %.4g Pa、500 K -> 氚分子密度 %.6g mol/m^3' % (P_T, n_T_gas))
J = PHI_EFF * math.sqrt(P_T) / PERM_THICK
print('渗透窗通量 J = Phi sqrt(p)/d = %.6g mol/(m^2 s) = %.6g g/(m^2 day)'
      % (J, J * M_T * SEC_DAY))
a_window = t_inject / (J * M_T)
print('  要提取 %.6g g/s 需渗透窗面积 %.6g m^2（理想混合、无回流）' % (t_inject, a_window))
J_hi = PHI_EFF * math.sqrt(1.0e4) / 1.0e-4
a_window_hi = t_inject / (J_hi * M_T)
print('  提高到 p_T = 1e4 Pa、壁厚 0.1 mm：通量 %.6g g/(m^2 day) -> 只需 %.6g m^2'
      % (J_hi * M_T * SEC_DAY, a_window_hi))
print('  为提高回收率可两级串联（各 %.4g -> 总 %.6g），但面积不降：面积由通量决定'
      % (0.9, 1 - (1 - 0.9) ** 2))
print('  反过来：10 m^2 的窗（p_T = 10 Pa、2 mm）每天只能提出 %.6g g，占需求的 %.6g%%'
      % (J * M_T * SEC_DAY * 10.0, 100 * J * M_T * SEC_DAY * 10.0 / (t_inject * SEC_DAY)))
print()

# ================= AV.5 渗透与滞留 =================
print('==== AV.5 渗透与滞留：库存不是一个数，是一条链 ====')
J_hot = PHI_EFF * math.sqrt(P_T) / PERM_THICK
perm_year = J_hot * A_HOT * M_T * SEC_YEAR / PRF
print('热表面 %.0f m^2、Phi_eff = %.4g、PRF = %.0f -> 渗透 %.6g g/a（无壁垒为 %.6g g/a）'
      % (A_HOT, PHI_EFF, PRF, perm_year, perm_year * PRF))
print('  占燃烧量 %.6g%%；占注入量 %.6g%%'
      % (100 * perm_year / (t_burn_g_s * SEC_YEAR), 100 * perm_year / (t_inject * SEC_YEAR)))
lam = math.log(2.0) / (T_HALF_Y * SEC_YEAR)
inv_g = INV_KG * 1e3
n_atoms = inv_g / m_T_g
act = lam * n_atoms
decay_g_day = inv_g * (1.0 - math.exp(-lam * SEC_DAY))
p_decay = act * E_BETA_KEV * 1e3 * E_CHG
print('半衰期 %.2f a -> lambda = %.6g s^-1' % (T_HALF_Y, lam))
print('库存 %.4g kg -> 原子数 %.6g、活度 %.6g Bq（%.6g TBq）' % (INV_KG, n_atoms, act, act / 1e12))
print('  衰变热 = A x %.1f keV = %.6g W = %.6g kW' % (E_BETA_KEV, p_decay, p_decay / 1e3))
print('  每天衰变掉 %.6g g（不加料就每天少这么多）' % decay_g_day)
print('  衰变掉的氚 = 燃烧量的 %.6g%%' % (100 * decay_g_day / (t_burn_g_s * SEC_DAY)))
print('共沉积滞留（附录 AJ）%.3f g 相当于库存的 %.6g%%' % (CO_DEP_AJ, 100 * CO_DEP_AJ / inv_g))
print()

# ================= AV.6 氚自持与倍增时间 =================
print('==== AV.6 自持：TBR 的及格线在哪里 ====')
tbr_min = (t_burn_g_s * SEC_DAY + decay_g_day) / (t_burn_g_s * SEC_DAY)
print('最小自持 TBR =（燃烧 + 衰变）/燃烧 = %.6g' % tbr_min)
print('  含提取效率 %.4g：TBR_min = %.6g / eta = %.6g' % (ETA_EXT, tbr_min, tbr_min / ETA_EXT))
for tbr in (1.02, 1.05, 1.08, 1.10):
    net = tbr * t_burn_g_s * SEC_DAY - t_burn_g_s * SEC_DAY - decay_g_day
    if net > 0:
        print('  TBR = %.4g -> 净增 %.6g g/day -> 库存 %.4g kg 的倍增时间 %.6g a'
              % (tbr, net, INV_KG, inv_g / net / 365.25))
    else:
        print('  TBR = %.4g -> 净增 %.6g g/day（负值：库存永远补不回来）' % (tbr, net))
print()

# ================= AV.7 衡算与安全 =================
print('==== AV.7 衡算与安全：0.5% 的精度意味着什么 ====')
print('衡算不确定度 %.4g -> 库存 %.4g kg 上 %.6g g' % (ACC_ERR, INV_KG, inv_g * ACC_ERR))
print('  相当于 %.6g 天的燃烧量；目标 %.4g -> %.6g g（%.6g 天）'
      % (inv_g * ACC_ERR / (t_burn_g_s * SEC_DAY), ACC_ERR2, inv_g * ACC_ERR2,
         inv_g * ACC_ERR2 / (t_burn_g_s * SEC_DAY)))
acc_need = (t_burn_g_s * SEC_DAY) / inv_g
print('  若要求衡算误差 <= 1 天燃烧量，精度必须优于 %.6g（%.6g%%）' % (acc_need, 100 * acc_need))
act_per_g_T = act / inv_g
act_per_g_hto = act_per_g_T * 3.0 / 20.0
print('氚比活度 = %.6g Bq/g（T）；氚水 HTO 比活度 = %.6g Bq/g（HTO）'
      % (act_per_g_T, act_per_g_hto))
print('按饮用水限值 %.0e Bq/L，1 g HTO 需 %.6g m^3 水稀释；1 g 氚（折成 HTO）需 %.6g m^3'
      % (HTO_LIMIT, act_per_g_hto / HTO_LIMIT / 1e3, act_per_g_T / HTO_LIMIT / 1e3))
print('库存 %.4g kg 相当于 %.6g 天的燃烧量（即「库存天数」）'
      % (INV_KG, inv_g / (t_burn_g_s * SEC_DAY)))
print('库存 %.4g kg 的总活度若全部以 HTO 形式释放，需要 %.6g m^3 水稀释'
      % (INV_KG, (act) / HTO_LIMIT / 1e3))
print()

# ================= AV.8 来源、分级与版权说明（无数值） =================
print('==== AV.8 来源、分级与版权说明（本节只列来源与边界，不复算） ====')
print()

# ================= AV.9 练习与解答 =================
print('==== AV.9 估算练习与解答 ====')
# 练习 1
p_new = 800.0e6
r_new = p_new / E_FUS_J
ans1_day = r_new * m_T_g * SEC_DAY
ans1_year = r_new * m_T_g * SEC_YEAR / 1e3
ans1_spec = ans1_day / (p_new / 1e6)
print('[练习 1] P_fus = 800 MW -> 氚 %.6g g/day、%.6g kg/a、比耗 %.6g g/(MW day)'
      % (ans1_day, ans1_year, ans1_spec))
print('  解答 1：R = P/E = %.6g s^-1 -> %.6g g/day、%.6g kg/a、%.6g g/(MW day)'
      % (r_new, ans1_day, ans1_year, ans1_spec))
# 练习 2
inv2 = 8.0e3
tbr2 = 1.08
net2 = tbr2 * t_burn_g_s * SEC_DAY - t_burn_g_s * SEC_DAY - decay_g_day * (inv2 / inv_g)
ans2 = inv2 / net2 / 365.25
print('[练习 2] TBR = 1.08、库存 8 kg -> 净增 %.6g g/day、倍增时间 %.6g a' % (net2, ans2))
print('  解答 2：净增 = (TBR-1) x 燃烧 - 衰变 = %.6g g/day -> %.6g a'
      % (net2, ans2))
# 练习 3
eta3 = 0.85
ans3 = tbr_min / eta3
print('[练习 3] 提取效率 85%% -> TBR_min = %.6g / %.4g = %.6g' % (tbr_min, eta3, ans3))
print('  解答 3：TBR_min = (1 + 衰变/燃烧)/eta = %.6g' % ans3)
# 练习 4
prf4 = 200.0
ans4 = J_hot * A_HOT * M_T * SEC_YEAR / prf4
print('[练习 4] PRF = 200 -> 年渗透 %.6g g/a，占燃烧量 %.6g%%'
      % (ans4, 100 * ans4 / (t_burn_g_s * SEC_YEAR)))
print('  解答 4：%.6g g/a（壁垒加一倍，渗透减半）' % ans4)
# 练习 5
ans5a = n_stage(0.5, 0.95, 1.3)
ans5b = n_stage(0.5, 0.95, 1.5)
print('[练习 5] alpha = 1.3、x 0.5 -> 0.95：%.6g 板；alpha = 1.5 时 %.6g 板（差 %.6g 倍）'
      % (ans5a, ans5b, ans5a / ans5b))
print('  解答 5：N = ln[(x/(1-x))_P/(x/(1-x))_F]/ln(alpha) = %.6g 板 / %.6g 板' % (ans5a, ans5b))
# 练习 6
ans6a = inv_g * ACC_ERR2
ans6b = ans6a / (t_burn_g_s * SEC_DAY)
print('[练习 6] 衡算 0.5%% 于 10 kg -> %.6g g，相当于 %.6g 天燃烧量' % (ans6a, ans6b))
print('  解答 6：%.6g g、%.6g 天（要压到 1 天需 %.6g%% ，见 AV.7）' % (ans6a, ans6b, 100 * acc_need))
print()
print('==== 附：派生量补充（供正文引用前先在此打印）====')
print('TBR 设计 %.4g 下的年产生量 = %.6g kg/a；净增（含衰变）= %.6g g/day'
      % (TBR_DESIGN, t_bred_g_s * SEC_YEAR / 1e3,
         t_bred_g_s * SEC_DAY - t_burn_g_s * SEC_DAY - decay_g_day))
print('库存 10 kg 的循环时间（按注入量）= %.6g day' % (inv_g / (t_inject * SEC_DAY)))
print('循环一圈衰变掉的份额 = %.6g%%（= 循环时间 x lambda）'
      % (100 * (1 - math.exp(-lam * (inv_g / (t_inject * SEC_DAY)) * SEC_DAY))))
print('氚水分子量 20.0275 g/mol；1 g HTO 含氚 %.6g g' % (3 * M_T / (2 * 1.007825 + 16.0 + 3 * M_T)))
