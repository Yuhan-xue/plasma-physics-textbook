# -*- coding: utf-8 -*-
"""v1.13.0 附录 P（低温等离子体应用）的数值权威来源。

本脚本列出附录 P 里每一个数字的算法，正文与 check_P12 断言都必须与它一致。
运行：python tmp/appP_calc.py > tmp/appP_calc.txt
"""
import io
import math

E = 1.602176634e-19          # C
ME = 9.1093837015e-31        # kg
U = 1.66053906660e-27        # kg
EPS0 = 8.8541878128e-12      # F/m
KB = 1.380649e-23            # J/K
TORR = 133.322               # Pa
LOSCHMIDT = 2.6867805e25     # m^-3 (273.15 K, 101325 Pa)
SCCM = LOSCHMIDT * 1.6666667e-8   # 分子/s

M_AR = 39.948 * U
M_SI = 28.0855 * U
RHO_SI = 2329.0              # kg/m^3
N_SI = RHO_SI / M_SI

out = []
def p(s=''):
    out.append(s)

def debye(te_ev, ne):
    return math.sqrt(EPS0 * te_ev * E / (ne * E * E))

def cs(te_ev, m):
    return math.sqrt(te_ev * E / m)

def eth(te_ev, m):
    return math.sqrt(8 * te_ev * E / (math.pi * m))

def paschen_min(A, B, gamma):
    """V_b = B*pd/(ln(A pd) - ln(ln(1+1/gamma))); 返回 (Vmin, pd_min)"""
    k = math.log(math.log(1 + 1.0 / gamma))
    pdm = math.exp(k + 1) / A
    return B * pdm, pdm, k

def paschen(A, B, gamma, pd):
    k = math.log(math.log(1 + 1.0 / gamma))
    return B * pd / (math.log(A * pd) - k)

def child_sheath(te_ev, ne, vs):
    """Child 律鞘层厚度 s = (2/3) lam_D (2 V_s/T_e)^{3/4}"""
    return (2.0 / 3.0) * debye(te_ev, ne) * (2.0 * vs / te_ev) ** 0.75

p('=== 附录 P 数值权威来源（tmp/appP_calc.py）===')
p('常量：e=%.9e C  m_e=%.10e kg  u=%.11e kg  eps0=%.10e F/m  k_B=%.7e J/K' % (E, ME, U, EPS0, KB))
p('1 Torr=%.3f Pa  洛施密特=%.7e m^-3  1 sccm=%.5e 分子/s' % (TORR, LOSCHMIDT, SCCM))
p('m_Ar=%.6e kg  m_Si=%.6e kg  n_Si=%.6e m^-3' % (M_AR, M_SI, N_SI))
p()

p('--- P.1 参数地图：德拜长度、玻姆速度、离子通量与离子流密度 ---')
CONF = [
    ('集成电路刻蚀(CCP/ICP)', 3.0, 1e17),
    ('溅射与PECVD', 2.0, 1e16),
    ('等离子体医学(射流)', 1.5, 1e19),
    ('大气压DBD微放电通道', 2.0, 1e21),
    ('聚变边界(对照,附录N)', 15.0, 2e19),
]
p('%-26s %6s %10s %9s %11s %11s %11s' % ('工况', 'Te/eV', 'ne/m^-3', 'lam_D', 'c_s(Ar)', 'Gamma_i', 'j_i'))
for name, te, ne in CONF:
    ld = debye(te, ne)
    c = cs(te, M_AR)
    g = ne * c
    j = g * E
    p('%-26s %6.1f %10.2e %8.2f um %8.0f m/s %11.3e %8.2f A/m2' % (name, te, ne, ld * 1e6, c, g, j))
p('  (j_i 换算：1 A/m2 = 0.1 mA/cm2)')
for name, te, ne in CONF:
    p('  %s: j_i=%.3f mA/cm2' % (name, ne * cs(te, M_AR) * E * 0.1))
te0, ne0 = 3.0, 1e17
p('刻蚀基准工况 Te=%.0f eV, ne=%.0e m^-3: lam_D=%.4e m = %.2f um; c_s=%.1f m/s; Gamma=%.4e m^-2 s^-1; j_i=%.3f mA/cm2'
  % (te0, ne0, debye(te0, ne0), debye(te0, ne0) * 1e6, cs(te0, M_AR), ne0 * cs(te0, M_AR), ne0 * cs(te0, M_AR) * E * 0.1))
p('非平衡比：T_e/T_gas = 3 eV / %.4f eV = %.1f 倍' % (KB * 300 / E, 3.0 / (KB * 300 / E)))
vth = eth(2.0, ME)
p('电子热运动速率(2 eV)=%.4e m/s，热流 Gamma_e=n v_th/4 (1e17)=%.4e m^-2 s^-1，j_e=%.1f A/m2 = %.0f mA/cm2'
  % (vth, 1e17 * vth / 4, 1e17 * vth / 4 * E, 1e17 * vth / 4 * E * 0.1))
p('悬浮电位 e phi_f/kTe = 0.5 ln(m_i/(2 pi m_e)) = %.4f -> phi_f(2 eV)=%.2f V'
  % (0.5 * math.log(M_AR / (2 * math.pi * ME)), 0.5 * math.log(M_AR / (2 * math.pi * ME)) * 2.0))
p()

p('--- P.2 汤森雪崩与帕申定律 (A=15 /(cm Torr), B=365 V/(cm Torr), gamma=0.01) ---')
A, B, GAM = 15.0, 365.0, 0.01
vmin, pdmin, kk = paschen_min(A, B, GAM)
p('ln(ln(1+1/gamma))=%.4f ; 帕申最小值 V_b,min=%.2f V 出现在 pd=%.4f Torr*cm' % (kk, vmin, pdmin))
for pd in (0.5, 0.8363, 1.0, 2.0, 5.0, 10.0, 76.0):
    p('  pd=%7.4f Torr*cm -> V_b=%.2f V' % (pd, paschen(A, B, GAM, pd)))
alpha = A * 760 * math.exp(-B * 760 / 3.0e4)
p('大气压空气经验击穿场 E_br=3.0e6 V/m (30 kV/cm)：1 mm 间隙需 %.0f V' % (3.0e6 * 1e-3))
N_GAS = 101325.0 / (KB * 300)
p('1 atm 300 K 数密度 N=%.4e m^-3 ; 约化场 E/N=%.4e V*m^2 = %.1f Td (1 Td=1e-21 V m^2)'
  % (N_GAS, 3.0e6 / N_GAS, 3.0e6 / N_GAS * 1e21))
p('E/p=%.1f V/(cm*Torr)' % (3.0e4 / 760.0))
p('汤森第一系数 alpha=A p exp(-B p/E) 在 1 atm、30 kV/cm 处 = %.4f cm^-1（对应 alpha*d=%.2f，需 >=ln(1/gamma)=%.2f 才自持）'
  % (alpha, alpha * 0.1, math.log(1 / GAM)))
p('同一帕申公式外推到 76 Torr*cm 给出 %.0f V，是经验击穿电压 3000 V 的 %.2f 倍 -> 大气压必须用经验场强，不能外推'
  % (paschen(A, B, GAM, 76.0), paschen(A, B, GAM, 76.0) / 3000.0))
p()

p('--- P.3 介质阻挡放电 (DBD) 电荷与功率账 ---')
EPS_R, D_DIELECTRIC, AREA = 10.0, 1.0e-3, 0.01
C_D = EPS0 * EPS_R * AREA / D_DIELECTRIC
p('介质电容 C_d = eps0 eps_r A/d = %.6e F = %.4f nF' % (C_D, C_D * 1e9))
for UA in (4000.0, 5000.0):
    UB = 3000.0
    Q = C_D * (UA - UB)
    POW = 4 * 1.0e4 * C_D * UB * (UA - UB)
    per_cycle = POW / 1.0e4
    p('U_a=%.0f V, U_b=%.0f V: 传输电荷/半周期=%.4e C=%.4f uC  P=4 f C_d U_b (U_a-U_b)=%.1f W=%.3f W/cm2  每周期能量=%.2f mJ'
      % (UA, UB, Q, Q * 1e6, POW, POW / (AREA * 1e4), per_cycle * 1e3))
UA, UB, FIL_PER_CM2, T_NS, R_FIL = 4000.0, 3000.0, 10, 10.0, 1.0e-4
Q = C_D * (UA - UB)
POW = 4 * 1.0e4 * C_D * UB * (UA - UB)
N_FIL = FIL_PER_CM2 * AREA * 1e4          # 每个半周期的丝数
q_fil = Q / N_FIL
e_fil = (POW / 1.0e4) / (2 * N_FIL)
I_fil = q_fil / (T_NS * 1e-9)
A_fil = math.pi * R_FIL ** 2
p('取 U_a=%.0f V：每半周期丝数=%.0f（%.0f 根/cm2 x %.0f cm2），单丝电荷=%.4e C=%.3f nC'
  % (UA, N_FIL, FIL_PER_CM2, AREA * 1e4, q_fil, q_fil * 1e9))
p('  单丝能量=%.3e J=%.1f uJ ; 若持续 %.0f ns，峰值电流=%.4f A，丝内电流密度=%.3e A/m2=%.0f A/cm2'
  % (e_fil, e_fil * 1e6, T_NS, I_fil, I_fil / A_fil, I_fil / A_fil * 1e-4))
V_DIS = AREA * D_DIELECTRIC
v_fil = 4.0 / 3.0 * math.pi * R_FIL ** 3
fill_space = N_FIL * v_fil / V_DIS
fill_time = T_NS * 1e-9 * 1.0e4
p('  放电体积=%.3e m3; 单丝体积=%.3e m3; 空间填充因子=%.3e; 时间占空比=%.3e; 时空填充因子=%.3e'
  % (V_DIS, v_fil, fill_space, fill_time, fill_space * fill_time))
p('  -> 能量只落在极少时空的丝里，重粒子气体温度仍在 300-400 K（非平衡的几何原因）')
p()

p('--- P.4 鞘层与离子轰击：刻蚀方向性 ---')
for VS in (100.0, 200.0, 400.0):
    s = child_sheath(3.0, 1e17, VS)
    dth = math.sqrt(3.0 / (2.0 * VS))
    p('V_s=%5.0f V: 鞘层厚度 s=%.4e m=%.3f mm ; 离子入射角展宽 dtheta=sqrt(T_e/2E_i)=%.4f rad=%.2f deg'
      % (VS, s, s * 1e3, dth, math.degrees(dth)))
v100 = math.sqrt(2 * E * 100.0 / M_AR)
p('100 eV 氩离子速度=%.4e m/s' % v100)
P_PA = 1.3332
N_AR = P_PA / (KB * 300)
SIG_CX = 1.0e-18
p('10 mTorr=%.4f Pa: n_Ar=%.4e m^-3; 取 sigma_CX=1.0e-18 m2, 平均自由程 lam=%.4e m=%.2f mm'
  % (P_PA, N_AR, 1.0 / (N_AR * SIG_CX), 1.0 / (N_AR * SIG_CX) * 1e3))
p('  lam/s(100 V)=%.2f -> 鞘层基本无碰撞，方向性得以保留' % ((1.0 / (N_AR * SIG_CX)) / child_sheath(3.0, 1e17, 100.0)))
GAMMA_ETCH = 1e17 * cs(3.0, M_AR)
for Y in (1.0, 0.1):
    rate = GAMMA_ETCH * Y / N_SI
    p('离子通量 Gamma=%.4e m^-2 s^-1, 产额 Y=%.1f: 刻蚀速率=%.4e m/s=%.2f nm/s=%.1f nm/min'
      % (GAMMA_ETCH, Y, rate, rate * 1e9, rate * 60 * 1e9))
r12 = GAMMA_ETCH * 1.2 / N_SI
p('练习 2：Y=1.2（反应离子刻蚀 Si/Cl）-> %.4e m/s = %.1f nm/min' % (r12, r12 * 60 * 1e9))
p()

p('--- P.5 表面反应：PECVD 与溅射沉积 ---')
FLOW, DISS, WAFER_R = 20.0, 0.10, 0.10
A_WAFER = math.pi * WAFER_R ** 2
n_sih4 = FLOW * SCCM
n_si_dep = n_sih4 * DISS
flux_dep = n_si_dep / A_WAFER
rate = flux_dep / N_SI
p('%.0f sccm SiH4 = %.4e 分子/s; 解离率 %.0f%% -> 沉积 Si 原子 %.4e /s' % (FLOW, n_sih4, DISS * 100, n_si_dep))
p('200 mm 晶圆面积=%.6f m2; 沉积通量=%.4e m^-2 s^-1; 生长速率=%.4e m/s=%.3f nm/s=%.1f nm/min'
  % (A_WAFER, flux_dep, rate, rate * 1e9, rate * 60 * 1e9))
t_1um = 1e-6 / rate
n_needed = 1e-6 * A_WAFER * N_SI
t_ideal = n_needed / n_si_dep
p('  长 1 um 需要 %.1f s = %.2f min; 若 100%% 利用只需 %.1f s -> 硅烷利用率=%.1f%%'
  % (t_1um, t_1um / 60, t_ideal, t_ideal / t_1um * 100))
J_TARGET = 100.0
Y_SP = 0.5
flux_ion = J_TARGET / E
rate_sp = flux_ion * Y_SP / N_SI
p('溅射：靶面离子流密度 %.0f A/m2 -> 离子通量 %.4e m^-2 s^-1; Y=%.1f -> 沉积速率 %.4e m/s=%.1f nm/min'
  % (J_TARGET, flux_ion, Y_SP, rate_sp, rate_sp * 60 * 1e9))
J_CM2, A_TARGET, HOURS = 0.010, 300.0, 100.0
I_TOT = J_CM2 * A_TARGET
Q_TOT = I_TOT * HOURS * 3600.0
n_ion = Q_TOT / E
n_erode = n_ion * Y_SP
m_erode = n_erode * 63.546 / 6.02214076e23
p('靶侵蚀：%.3f A/cm2 x %.0f cm2 = %.2f A，%.0f h -> 电荷 %.4e C，离子 %.4e 个，Y=%.1f -> 铜 %.1f g'
  % (J_CM2, A_TARGET, I_TOT, HOURS, Q_TOT, n_ion, Y_SP, m_erode))
p()

p('--- P.6 等离子体医学与农业：剂量学 ---')
POW_JET, T_SEC, A_TREAT = 1.0, 60.0, 1.0e-4
En = POW_JET * T_SEC
p('射流 %.1f W x %.0f s = %.0f J，处理 %.0f cm2 -> 剂量 %.4e J/m2 = %.3f J/cm2'
  % (POW_JET, T_SEC, En, A_TREAT * 1e4, En / A_TREAT, En / (A_TREAT * 1e4)))
CONC, VOL = 100e-6, 1e-3
n_h2o2 = CONC * VOL * 6.02214076e23
e_per = En / n_h2o2
p('液相 H2O2 浓度 %.0f uM，体积 %.0f mL -> %.4e 个分子；每个分子分摊能量 %.4e J=%.1f eV'
  % (CONC * 1e6, VOL * 1e3, n_h2o2, e_per, e_per / E))
p('  与 O-O 键能 5.12 eV 之比 = %.4f%% -> 99.9%% 以上的电能变成了热和紫外，化学转化只是零头'
  % (5.12 / (e_per / E) * 100))
E_NO_EV = 4.35
M_N = 14.0067
kg_n_per_mol_no = M_N / 1000.0
J_per_mol = E_NO_EV * E * 6.02214076e23
p('固氮：每个 NO 分子 %.2f eV -> %.4e J/mol = %.2f MJ/mol；1 kg N 对应 %.1f mol NO -> %.1f MJ/kg N'
  % (E_NO_EV, J_per_mol, J_per_mol / 1e6, 1000.0 / kg_n_per_mol_no, J_per_mol / 1e6 * 1000.0 / kg_n_per_mol_no))
p()

p('--- P.7 环境与催化：能量效率硬账 ---')
G_KWH = 100.0
M_O3 = 47.9982
n_o3 = G_KWH / M_O3 * 6.02214076e23
e_o3 = 3.6e6 / n_o3
p('臭氧产率 %.0f g/kWh -> 1 kWh 产 %.4e 个 O3，每个分子分摊 %.4e J=%.2f eV' % (G_KWH, n_o3, e_o3, e_o3 / E))
p('  O2 解离能 5.12 eV -> 理论极限效率 %.1f%%' % (5.12 / (e_o3 / E) * 100))
SIE = 300.0
PPM = 1000e-6
n_tol = PPM * LOSCHMIDT * 1e-3
e_tol = SIE / n_tol
p('VOC：比能量输入 SIE=%.0f J/L，甲苯 %.0f ppmv -> 每升 %.4e 个分子，每个分子 %.4e J=%.1f eV'
  % (SIE, PPM * 1e6, n_tol, e_tol, e_tol / E))
p('  （同 SIE 下浓度降到 100 ppmv，每分子代价升到 %.1f eV —— 能量代价与浓度成反比）' % (e_tol / E * 10))
e_co2 = 2.9
p('CO2 解离 CO2 -> CO + 1/2 O2 吸热 %.1f eV/分子；等离子体实际能耗常见 10-100 eV/分子 -> 热力学效率仅 %.0f%%-%.0f%%'
  % (e_co2, e_co2 / 100 * 100, e_co2 / 10 * 100))
p()

p('--- P.8 六道估算练习答案 ---')
p('练 1 帕申最小：V_b,min=%.1f V, pd=%.4f Torr*cm; 大气压 1 mm 外推 %.0f V vs 经验 3000 V（%.2f 倍）'
  % (vmin, pdmin, paschen(A, B, GAM, 76.0), paschen(A, B, GAM, 76.0) / 3000.0))
p('练 2 刻蚀速率：Y=1.2 -> %.1f nm/min；Y=0.1（纯溅射）-> %.1f nm/min，选择比 = Y 之比 = %.0f'
  % (r12 * 60 * 1e9, GAMMA_ETCH * 0.1 / N_SI * 60 * 1e9, 1.2 / 0.1))
p('练 3 鞘层：100 V -> s=%.3f mm, dtheta=%.2f deg; 400 V -> s=%.3f mm, dtheta=%.2f deg; lam=%.2f mm 故 s<lam'
  % (child_sheath(3.0, 1e17, 100.0) * 1e3, math.degrees(math.sqrt(3.0 / 200.0)),
     child_sheath(3.0, 1e17, 400.0) * 1e3, math.degrees(math.sqrt(3.0 / 800.0)), 1.0 / (N_AR * SIG_CX) * 1e3))
p('练 4 PECVD：1 um 需 %.1f min，硅烷利用率 %.1f%%' % (t_1um / 60, t_ideal / t_1um * 100))
p('练 5 臭氧：%.0f g/kWh -> %.2f eV/分子，极限效率 %.1f%%，产 1 kg 需 %.1f kWh'
  % (G_KWH, e_o3 / E, 5.12 / (e_o3 / E) * 100, 1000.0 / G_KWH))
p('练 6 剂量与固氮：射流 1 W 处理 2 cm2 60 s -> %.1f J/cm2；NO %.2f eV/分子 -> %.1f MJ/kg N'
  % (En / (2.0 * A_TREAT * 1e4), E_NO_EV, J_per_mol / 1e6 * 1000.0 / kg_n_per_mol_no))
p()
p('=== 完：以上数字是正文与门禁的唯一来源 ===')

io.open('tmp/appP_calc.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('\n'.join(out))
