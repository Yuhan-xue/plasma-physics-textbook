# -*- coding: utf-8 -*-
"""附录 M 数值的独立复算（第二版，修正常数与速度约定）。

与 tmp/appM_calc.py 的差别（这两处都是自查抓出的真错）：
 (1) 氘质量改用原子质量单位 m_u 换算 m_D = 2.013553 u，
     旧脚本用 A * m_p（把 2.0136 u 当成 2 m_p），使所有氘的 v_b / f_ci / B_res 偏大约 0.7%；
 (2) 热速度统一采用本书正文约定 v_th = sqrt(2 T / m)（附录 L 第 1 章同此约定），
     旧脚本用 sqrt(T/m)，使 LHCD 的 k_par / lambda_par 差 sqrt(2) 倍。
其余公式与旧脚本一致，结果逐条对照在本文件末尾打印。
"""
import numpy as np

q = 1.602176634e-19          # C
m_e = 9.1093837015e-31       # kg
m_u = 1.66053906660e-27      # kg
m_p = 1.67262192369e-27      # kg
m_D = 2.013553212745 * m_u   # kg, 氘原子（束是中性原子，电离后近似取原子质量）
eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6
keV = 1.0e3 * q
MeV = 1.0e6 * q

out = []
def p(label, value, unit=""):
    out.append(f"{label:52s} = {value:.6g} {unit}")

# ---------- M.2 欧姆加热 ----------
def eta_spitzer(Te_eV, lnL=17.0):
    return 5.2e-5 * lnL / Te_eV ** 1.5

eta_1keV = eta_spitzer(1.0e3)
p("M.2 eta(1 keV, lnL=17)", eta_1keV, "ohm m")
p("M.2 P_ohm(j=1 MA/m2)", eta_1keV * (1e6) ** 2, "W/m3")
p("M.2 tau_R(a=2 m)", mu0 * 2.0 ** 2 / eta_1keV, "s")
P_fus = 3.10181e6  # 附录 L 第 10 章
P_alpha = P_fus / 5.0
p("M.2 P_alpha(附录 L)", P_alpha, "W/m3")
p("M.2 P_alpha/P_ohm", P_alpha / (eta_1keV * 1e12), "")
p("M.2 eta(10 keV)/eta(1 keV)", eta_spitzer(1.0e4) / eta_1keV, "")

# ---------- M.3 中性束 ----------
def v_of_E(E_J, m):
    return np.sqrt(2.0 * E_J / m)

vb_100 = v_of_E(100 * keV, m_D)
vb_1MeV = v_of_E(1 * MeV, m_D)
p("M.3 v_b(100 keV D)", vb_100, "m/s")
p("M.3 v_b(1 MeV D)", vb_1MeV, "m/s")
sigma_eff = 1e-20
n_e = 1e20
p("M.3 lambda_att(sigma=1e-20, n=1e20)", 1.0 / (n_e * sigma_eff), "m")
# E_c = 14.8 * T_e[keV] * A_b * (sum n_i Z_i^2/A_i / n_e)^(2/3)  keV
def E_c(Te_keV, A_b, frac):
    return 14.8 * Te_keV * A_b * frac ** (2.0 / 3.0)
p("M.3 E_c(Te=10 keV, D, frac=0.5)", E_c(10.0, 2.0, 0.5), "keV")
p("M.3 E_c(Te=1 keV, D, frac=0.5)", E_c(1.0, 2.0, 0.5), "keV")
p("M.3 E_c/E_b(100 keV)", E_c(10.0, 2.0, 0.5) / 100.0, "")

# ---------- M.4 射频 ----------
def f_ci(B, m):
    return q * B / (2 * np.pi * m)
def f_ce(B):
    return q * B / (2 * np.pi * m_e)
def B_res(f, m):
    return 2 * np.pi * f * m / q

B0 = 5.3
fci = f_ci(B0, m_D)
fce = f_ce(B0)
p("M.4 f_ci(D, B=5.3 T)", fci, "Hz")
p("M.4 f_ci(D, B=5.3 T)/1e6", fci / 1e6, "MHz")
p("M.4 f_ce(B=5.3 T)", fce, "Hz")
p("M.4 f_ce/B", fce / B0 / 1e9, "GHz/T")
p("M.4 f_LH=sqrt(f_ci f_ce)", np.sqrt(fci * fce), "Hz")
p("M.4 B_res(ICRH 40 MHz, D)", B_res(40e6, m_D), "T")
p("M.4 B_res(ICRH 40 MHz, H)", B_res(40e6, m_p), "T")
p("M.4 B_res(ICRH 30 MHz, D)", B_res(30e6, m_D), "T")
p("M.4 B_res(ICRH 50 MHz, D)", B_res(50e6, m_D), "T")
Te = 10 * keV
v_te = np.sqrt(2 * Te / m_e)
p("M.4 v_te(10 keV)=sqrt(2T/m)", v_te, "m/s")
f_LH = 2.45e9
w_LH = 2 * np.pi * f_LH
for mult in (2.0, 3.0):
    vph = mult * v_te
    kpar = w_LH / vph
    p(f"M.4 k_par(v_phi={mult:g} v_te)", kpar, "1/m")
    p(f"M.4 lambda_par(v_phi={mult:g} v_te)", 2 * np.pi / kpar, "m")

# ---------- M.5 运行极限 ----------
def n_G(Ip_MA, a):
    return Ip_MA / (np.pi * a ** 2)
p("M.5 n_G(Ip=15 MA, a=2 m)", n_G(15, 2.0), "1e20 m^-3")
def beta_N(beta_pct, a, B, Ip_MA):
    return beta_pct * a * B / Ip_MA
p("M.5 beta_N(beta=2.5%, a=2, B=5.3, I=15)", beta_N(2.5, 2.0, 5.3, 15), "% m T/MA")
def P_LH(n20, B, S):
    return 0.0488 * n20 ** 0.717 * B ** 0.803 * S ** 0.941
p("M.5 P_LH Martin(n20=1, B=5.3, S=683)", P_LH(1.0, 5.3, 683.0), "MW")
p("M.5 P_LH Martin(n20=0.5, B=2.5, S=100)", P_LH(0.5, 2.5, 100.0), "MW")
R0 = 6.2
for name, gam in (("LHCD", 0.3), ("NBI", 0.3), ("ECCD", 0.05)):
    p(f"M.5 P_RF(1 MA, gamma={gam:g}, {name})", 1.0 * 1.0 * R0 / gam, "MW")

# ---------- M.6 练习答案 ----------
p("M.6 [1] eta(10 keV)", eta_spitzer(1.0e4), "ohm m")
p("M.6 [1] P_ohm(j=1 MA/m2, 10 keV)", eta_spitzer(1.0e4) * 1e12, "W/m3")
p("M.6 [2] v_b(1 MeV D)", vb_1MeV, "m/s")
p("M.6 [3] B_res(170 GHz, e)", B_res(170e9, m_e), "T")
p("M.6 [4] v_ph(3 v_te, 10 keV)", 3 * v_te, "m/s")
p("M.6 [4] lambda_par(2.45 GHz, v_ph=3v_te)", 2 * np.pi * 3 * v_te / w_LH, "m")
n_G_ITER = n_G(15, 2.0) * 1e20
p("M.6 [5] f_pe(n=n_G(15MA,2m))", q * np.sqrt(n_G_ITER / (eps0 * m_e)) / (2 * np.pi), "Hz")
p("M.6 [6] P_alpha 6.2036e5 (由附录 L)", P_alpha, "W/m3")
p("M.6 [6] P_ext(=0.2 P_alpha, 若 P_loss=1.2 P_alpha)", 0.2 * P_alpha, "W/m3")

print("\n".join(out))
