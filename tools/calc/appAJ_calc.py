# -*- coding: utf-8 -*-
"""附录 AJ（加料与排灰）的唯一数值来源。

只 print，不写文件。用法：
    python tmp\\appAJ_calc.py | Out-File -Encoding utf8 tmp\\appAJ_calc.txt

三级数值：
  现算级 —— 本脚本打印（正文与门禁只允许引用本文件的数）
  声明级 —— 脚本顶部输入或明确标注为量级/文献值的常数
  文献级 —— 只作定性对照，不参与任何数值断言

工作点与附录 M/S/U/V/W/X/Y/Z/AD/AI 同一点：n_e=1e20 m^-3、T_e=9.194 keV、B=5.3 T、
a=2.0 m、R=6.2 m、V=832.2 m^3、P_fus=556.5 MW、P_alpha=111.3 MW。
氚质量取附录 Z 的 m_T=3.0160492 u；库存上限取附录 Z 的 2.0 kg。
"""
import math

# ---------- CODATA 2018 ----------
c = 2.99792458e8
e = 1.602176634e-19
mp = 1.67262192369e-27
u = 1.66053906660e-27
kB = 1.380649e-23
MeV = 1.602176634e-13
year = 3.15576e7           # s，儒略年
day = 86400.0

print('=' * 72)
print('附录 AJ 加料与排灰：数值来源（tmp/appAJ_calc.py）')
print('=' * 72)

# ---------- DECLARED 输入 ----------
n_e = 1.0e20              # m^-3   芯部电子密度（与全书 D1 工作点同点）
T_e_eV = 9194.0           # eV     电子温度
B0 = 5.3                  # T      环向磁场
a = 2.0                   # m      小半径
R0 = 6.2                  # m      大半径
V = 832.2                 # m^3    等离子体体积
P_fus = 556.5e6           # W      聚变功率
P_alpha = 111.3e6         # W      alpha 功率
E_fus = 17.6              # MeV    每次 D-T 反应释放的总能量（声明）
tau_E = 3.7               # s      能量约束时间（声明，与附录 AD/AE 同点）
tau_p = 2.5               # s      粒子约束时间（声明）
R_rec = 0.95              # -      再循环系数（声明）
eps_gas = 0.10            # -      气体加料效率（声明为量级值）
eps_pel = 0.60            # -      弹丸加料效率（声明为量级值）
r_pel = 1.5e-3            # m      弹丸半径（声明）
rho_ice = 255.0           # kg/m^3 DT 冰密度（与附录 AI 同点）
v_pel = 300.0             # m/s    弹丸速度（声明为量级值）
tau_abl = 4.0e-4          # s      工作点弹丸烧蚀时间（声明，NGS 量级值）
nHe_target = 0.05         # -      氦灰份额控制目标 n_He/n_e（声明）
f_tauHe = 5.0             # -      声明 tau_He*/tau_E
P_div = 5.0               # Pa     偏滤器中性气体压力（声明）
S_pump = 100.0            # m^3/s  抽气速率（声明）
T_gas = 300.0             # K      抽气气体温度（声明）
He_enrich = 5.0           # -      偏滤器处氦富集因子（声明）
f_ret = 0.01              # -      氚滞留份额（声明）
inv_T_lim = 2.0           # kg     氚库存上限（与附录 Z 的 2.0 kg 同口径）
T_half = 12.32            # 年     氚半衰期（声明，与附录 Z 同）
m_T = 3.0160492 * u       # kg     氚原子质量（与附录 Z 同）
m_i = 2.5 * u             # kg     D-T 平均离子质量（与附录 AD 同口径）
I_p = 11.0844e6           # A      等离子体电流（与附录 AG 的 q95=3 工作点同点）

print('\n[DECLARED] 输入')
print('  n_e = %.6g m^-3, T_e = %.6g eV, B = %.6g T, a = %.6g m, R = %.6g m, V = %.6g m^3'
      % (n_e, T_e_eV, B0, a, R0, V))
print('  P_fus = %.6g MW, P_alpha = %.6g MW, E_fus = %.6g MeV, tau_E = %.6g s, tau_p = %.6g s'
      % (P_fus / 1e6, P_alpha / 1e6, E_fus, tau_E, tau_p))
print('  R_rec = %.6g, eps_gas = %.6g, eps_pel = %.6g' % (R_rec, eps_gas, eps_pel))
print('  弹丸 r = %.6g mm, rho = %.6g kg/m^3, v = %.6g m/s, tau_abl = %.6g s'
      % (r_pel * 1e3, rho_ice, v_pel, tau_abl))
print('  n_He/n_e 目标 = %.6g, tau_He*/tau_E = %.6g, 富集因子 = %.6g'
      % (nHe_target, f_tauHe, He_enrich))
print('  抽气 P_div = %.6g Pa, S_pump = %.6g m^3/s, T_gas = %.6g K' % (P_div, S_pump, T_gas))
print('  氚滞留份额 = %.6g, 库存上限 = %.6g kg, 半衰期 = %.6g 年, m_T = %.6g kg'
      % (f_ret, inv_T_lim, T_half, m_T))
print('  m_i = %.6g kg (= 2.5 u), I_p = %.6g MA' % (m_i, I_p / 1e6))

# ---------------- AJ.1 燃料要多少：粒子平衡的三项 ----------------
print('\n[AJ.1] 燃料要多少：粒子平衡的三项')
N_e = n_e * V
R_fus = P_fus / (E_fus * MeV)
burn_ion = 2.0 * R_fus
burn_T = R_fus * m_T * year
loss_trans = N_e / tau_p
loss_total = loss_trans + burn_ion
S_ext = (1.0 - R_rec) * loss_trans + burn_ion
f_burn_net = burn_ion / S_ext
f_burn_inj_gas = burn_ion / (S_ext / eps_gas)
f_burn_inj_pel = burn_ion / (S_ext / eps_pel)
print('  电子总数 N_e = n_e V = %.6g' % N_e)
print('  聚变反应率 R_fus = P_fus/E_fus = %.6g s^-1' % R_fus)
print('  燃烧消耗率（D 与 T 各一份）= %.6g s^-1' % burn_ion)
print('  氚燃耗 = %.6g kg/年（附录 Z 同口径 31.19 kg/年）' % burn_T)
print('  输运损失 N_e/tau_p = %.6g s^-1' % loss_trans)
print('  总损失 = %.6g s^-1' % loss_total)
print('  外部加料需求 S_ext = (1-R)N/tau_p + 燃烧 = %.6g s^-1' % S_ext)
print('  加料量 / 输运损失 = %.6g' % (S_ext / loss_trans))
print('  加料量 / 燃烧量 = %.6g' % (S_ext / burn_ion))
print('  净加料里被烧掉的份额 = %.6g（%.6g%%）' % (f_burn_net, f_burn_net * 100))
print('  注进去的粒子被烧掉的份额：气体加料 %.6g%%、弹丸加料 %.6g%%'
      % (f_burn_inj_gas * 100, f_burn_inj_pel * 100))

# ---------------- AJ.2 气体加料还是弹丸加料 ----------------
print('\n[AJ.2] 气体加料还是弹丸加料：效率差在哪')
S_inj_gas = S_ext / eps_gas
mol_gas = S_inj_gas / 2.0
Q_gas = mol_gas * kB * T_gas
m_dot_pel = (S_ext / eps_pel) / 2.0 * m_i
print('  气体加料注入率 = S_ext/eps_gas = %.6g s^-1' % S_inj_gas)
print('  D2 分子流率 = %.6g s^-1' % mol_gas)
print('  气体流量 = %.6g Pa m^3/s = %.6g mbar L/s' % (Q_gas, Q_gas * 10.0))
print('  弹丸加料注入率 = S_ext/eps_pel = %.6g s^-1' % (S_ext / eps_pel))
print('  弹丸质量流率 = %.6g kg/s = %.6g mg/s' % (m_dot_pel, m_dot_pel * 1e6))
print('  = %.6g g/day = %.6g kg/年' % (m_dot_pel * 1e3 * day, m_dot_pel * year))
print('  气体注入的粒子数是弹丸的 %.6g 倍' % (eps_pel / eps_gas))

# ---------------- AJ.3 弹丸烧蚀与穿透 ----------------
print('\n[AJ.3] 弹丸烧蚀与穿透：中性气体屏蔽标度')
m_pel = 4.0 / 3.0 * math.pi * r_pel ** 3 * rho_ice
N_at = m_pel / m_i
f_pel = m_dot_pel / m_pel
L_pen = v_pel * tau_abl
m_T_pel = 0.6 * m_pel
r2 = 2.0e-3
r3 = 3.0e-3
m_pel_r2 = 4.0 / 3.0 * math.pi * r2 ** 3 * rho_ice
m_pel_r3 = 4.0 / 3.0 * math.pi * r3 ** 3 * rho_ice
tau_r2 = tau_abl * (r2 / r_pel) ** (5.0 / 3.0)
tau_r3 = tau_abl * (r3 / r_pel) ** (5.0 / 3.0)
L_r2 = v_pel * tau_r2
L_r3 = v_pel * tau_r3
L_r3_fast = 1000.0 * tau_r3
v_need_15 = v_pel * (0.5 * a) / L_pen
v_need_3 = v_pel * (0.5 * a) / L_r3
print('  弹丸质量 = %.6g kg = %.6g mg（含氚 %.6g mg）' % (m_pel, m_pel * 1e6, m_T_pel * 1e6))
print('  弹丸原子数 = %.6g（附录 AD 的 1e21 原子 = 4.18155 mg 同口径）' % N_at)
print('  弹丸频率 = %.6g Hz（周期 %.6g s）' % (f_pel, 1.0 / f_pel))
print('  穿透深度 L_pen = v tau_abl = %.6g m = %.6g cm，占小半径 %.6g%%'
      % (L_pen, L_pen * 100, L_pen / a * 100))
print('  r = 2 mm: 质量 %.6g mg, tau_abl %.6g s, L_pen %.6g m' % (m_pel_r2 * 1e6, tau_r2, L_r2))
print('  r = 3 mm: 质量 %.6g mg, tau_abl %.6g s, L_pen %.6g m（占 a 的 %.6g%%）'
      % (m_pel_r3 * 1e6, tau_r3, L_r3, L_r3 / a * 100))
print('  r = 3 mm + v = 1000 m/s: L_pen = %.6g m（L/a = %.6g）' % (L_r3_fast, L_r3_fast / a))
print('  要打到 a/2 = 1 m：r=1.5 mm 需 v >= %.6g m/s；r=3 mm 需 v >= %.6g m/s'
      % (v_need_15, v_need_3))

# ---------------- AJ.4 氦灰：稀释与排出条件 ----------------
print('\n[AJ.4] 氦灰：稀释与排出条件')
tau_He = f_tauHe * tau_E
n_He = R_fus * tau_He / V
f_He = n_He / n_e
n_fuel = n_e - 2.0 * n_He
f_fuel = n_fuel / n_e
P_dilute = P_fus * f_fuel ** 2
tau_He_req = nHe_target * n_e * V / R_fus
ratio_req = tau_He_req / tau_E
n_He_div = He_enrich * nHe_target * n_e
print('  氦产生率 = R_fus = %.6g s^-1' % R_fus)
print('  tau_He*/tau_E = %.6g -> tau_He* = %.6g s' % (f_tauHe, tau_He))
print('  稳态氦密度 n_He = %.6g m^-3，份额 %.6g%%' % (n_He, f_He * 100))
print('  燃料离子密度 = n_e - 2 n_He = %.6g m^-3（%.6g%% n_e）' % (n_fuel, f_fuel * 100))
print('  稀释后聚变功率 = %.6g MW（原来的 %.6g%%）' % (P_dilute / 1e6, P_dilute / P_fus * 100))
print('  要把份额压到 %.6g%%：tau_He* <= %.6g s，即 tau_He*/tau_E <= %.6g'
      % (nHe_target * 100, tau_He_req, ratio_req))
print('  偏滤器处氦密度（富集 %.6g 倍）= %.6g m^-3' % (He_enrich, n_He_div))

# ---------------- AJ.5 抽气与氚存量 ----------------
print('\n[AJ.5] 抽气与氚存量')
Q_pump = P_div * S_pump / (kB * T_gas)
ratio_pump = Q_pump / R_fus
frac_He_need = R_fus / Q_pump
T_fuel = S_ext / 2.0 * m_T
T_year = T_fuel * year
ret_year = f_ret * T_year
t_limit = inv_T_lim / ret_year
decay_year = math.log(2.0) / T_half * inv_T_lim
TBR_req = (burn_T + ret_year + decay_year) / burn_T
print('  抽气粒子排出率 = P S/(k T) = %.6g s^-1' % Q_pump)
print('  是氦产生率的 %.6g 倍（抽气从来不是瓶颈）' % ratio_pump)
print('  抽出的气体里只需 %.6g%% 是氦，就够把灰排掉' % (frac_He_need * 100))
print('  外部加料里的氚流 = %.6g kg/s = %.6g kg/day = %.6g kg/年'
      % (T_fuel, T_fuel * day, T_year))
print('  滞留份额 %.6g%% -> 每年滞留 %.6g kg' % (f_ret * 100, ret_year))
print('  库存上限 %.6g kg -> %.6g 年（%.6g 个月）达到上限'
      % (inv_T_lim, t_limit, t_limit * 12))
print('  衰变损失 = %.6g kg/年（库存 %.6g kg）' % (decay_year, inv_T_lim))
print('  含滞留与衰变的自持下限 TBR_req = %.6g（附录 Z 只算衰变是 1.0036）' % TBR_req)
print('  滞留 / 燃耗 = %.6g%%' % (ret_year / burn_T * 100))

# ---------------- AJ.6 密度控制回路 ----------------
print('\n[AJ.6] 密度控制回路')
n_G = I_p / 1e6 / (math.pi * a ** 2) * 1e20
ratio_G = n_e / n_G
t_ramp = N_e / S_ext
print('  Greenwald 密度 n_G = %.6g m^-3（I_p = %.6g MA）' % (n_G, I_p / 1e6))
print('  n_e/n_G = %.6g（附录 AG 同点 1.1337）' % ratio_G)
print('  弹丸周期 = %.6g s，粒子约束时间 tau_p = %.6g s，控制带宽受谁限制很清楚'
      % (1.0 / f_pel, tau_p))
print('  若加料全部用于升密度：N_e/S_ext = %.6g s（密度翻倍要这么久）' % t_ramp)

# ---------------- AJ.7 练习解答 ----------------
print('\n[AJ.7] 练习解答')
# 1) tau_p 加倍
tau_p2 = 5.0
S_ext2 = (1.0 - R_rec) * N_e / tau_p2 + burn_ion
print('  1) tau_p = 5 s: S_ext = %.6g s^-1，是原来的 %.6g 倍（降 %.6g 倍）'
      % (S_ext2, S_ext2 / S_ext, S_ext / S_ext2))
# 2) r = 2 mm 弹丸
print('  2) r = 2 mm: 质量 %.6g mg，tau_abl = %.6g ms，L_pen = %.6g m，频率 %.6g Hz'
      % (m_pel_r2 * 1e6, tau_r2 * 1e3, L_r2, m_dot_pel / m_pel_r2))
# 3) n_He/n_e = 3%
tau_He_3 = 0.03 * n_e * V / R_fus
print('  3) 份额 3%%: tau_He* = %.6g s，tau_He*/tau_E = %.6g' % (tau_He_3, tau_He_3 / tau_E))
# 4) 5% 灰的稀释
n_fuel5 = n_e - 2.0 * nHe_target * n_e
print('  4) 份额 5%%: 燃料密度 = %.6g n_e，聚变功率 = %.6g MW'
      % (n_fuel5 / n_e, P_fus * (n_fuel5 / n_e) ** 2 / 1e6))
# 5) 抽气
Q_pump2 = 10.0 * 100.0 / (kB * T_gas)
print('  5) P_div = 10 Pa, S = 100 m^3/s: 排出率 = %.6g s^-1，是氦产生率的 %.6g 倍'
      % (Q_pump2, Q_pump2 / R_fus))
# 6) 滞留
print('  6) 1%% 滞留: 年滞留 %.6g kg，达 2 kg 上限需 %.6g 年（%.6g 个月）'
      % (ret_year, t_limit, t_limit * 12))

# ---- [AJ.8] 派生比值（正文里写出的比值必须在这里能找到） ----
print('\n[AJ.8] 派生比值')
print('  输运损失 / 燃烧消耗率 = %.6g' % (loss_trans / burn_ion))
print('  气体注入量 / 燃烧消耗率 = %.6g' % (S_inj_gas / burn_ion))
print('  弹丸 r 1.5->2 mm 的质量比 = %.6g' % (m_pel_r2 / m_pel))
print('  弹丸 r 1.5->2 mm 的穿透深度增幅 = %.6g%%' % ((L_r2 / L_pen - 1.0) * 100.0))
print('  灰 5%% 时的燃料份额 = %.6g，功率份额 = %.6g'
      % (n_fuel5 / n_e, (n_fuel5 / n_e) ** 2))
print('  净加料 / 燃耗 = %.6g（氚循环要处理的路过量）' % (S_ext / burn_ion))

print('\n' + '=' * 72)
print('以上为本附录正文与门禁 P42 可以引用的全部数值')
print('=' * 72)
