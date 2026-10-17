# -*- coding: utf-8 -*-
"""附录 AN 偏振与磁光诊断：唯一数值来源（只 print，不写文件）。

重生成：python tmp\\appAN_calc.py | Out-File -Encoding utf8 tmp\\appAN_calc.txt

口径约定：
  * 法拉第旋转：R/L 模的相位差 = (omega/c) dn L；偏振面的**旋转角**是它的一半。
  * 科顿-穆顿：X/O 模的相位差本身即观测量。
  * 精确解用「剖面对应的等效均匀密度」 n_e_eff = Int(ne dl)/L，与解析线积分同口径。
"""
import math

# ---------- 常数 ----------
c = 2.99792458e8
e = 1.602176634e-19
me = 9.1093837015e-31
eps0 = 8.8541878128e-12
h = 6.62607015e-34
a0 = 5.29177210903e-11
u = 1.66053906660e-27
mD = 2.013553212745 * u

# ---------- 声明 ----------
ne0 = 1.0e20
Te = 9194.0
Bt = 5.3
Bp = 0.5
R0 = 6.2
a = 2.0
lam_CO2 = 10.6e-6
lam_FIR = 118.8e-6
Eb = 100.0e3
Tbeam = 10.0
f_PEM = 50.0e3
dphi_res = 0.01
N_chord = 9
lam_Dalpha = 656.28e-9

L_center = 2 * a
ne_L_int = ne0 * 4 * a / 3.0      # 抛物线剖面的中心弦 Int(ne dl)
ne_eff = ne_L_int / L_center      # 等效均匀密度（与线积分同口径）

print('== 声明输入 ==')
print('ne0=%.6g m^-3, Te=%.6g eV, Bt=%.6g T, Bp=%.6g T, R0=%.6g m, a=%.6g m' % (ne0, Te, Bt, Bp, R0, a))
print('lam_CO2=%.6g m, lam_FIR=%.6g m, Eb=%.6g eV, Tbeam=%.6g eV, f_PEM=%.6g Hz' % (lam_CO2, lam_FIR, Eb, Tbeam, f_PEM))
print('dphi_res=%.6g deg, N_chord=%d, lam_Dalpha=%.6g m' % (dphi_res, N_chord, lam_Dalpha))
print('L_center=%.6g m, Int(ne dl)=%.6g m^-2, n_e_eff=%.6g m^-3' % (L_center, ne_L_int, ne_eff))

# ---------- [AN.1] 特征频率与本征模 ----------
wpe = math.sqrt(ne0 * e * e / (eps0 * me))
wpe_eff = math.sqrt(ne_eff * e * e / (eps0 * me))
wce_t = e * Bt / me
wce_p = e * Bp / me
w_CO2 = 2 * math.pi * c / lam_CO2
w_FIR = 2 * math.pi * c / lam_FIR
print()
print('== [AN.1] 特征频率与本征模 ==')
print('omega_pe(芯部) = %.6g rad/s' % wpe)
print('omega_pe(等效) = %.6g rad/s' % wpe_eff)
print('omega_ce(Bt=5.3T) = %.6g rad/s' % wce_t)
print('omega_ce(Bp=0.5T) = %.6g rad/s' % wce_p)
print('omega(CO2) = %.6g rad/s ; omega(FIR) = %.6g rad/s' % (w_CO2, w_FIR))
for tag, w in (('CO2', w_CO2), ('FIR', w_FIR)):
    X = (wpe_eff / w) ** 2
    Yt = wce_t / w
    Yp = wce_p / w
    print('%s: X_eff = %.6g ; Y(Bt) = %.6g ; Y(Bp) = %.6g' % (tag, X, Yt, Yp))


def n_par(X, Y, sign):
    return math.sqrt(1.0 - X / (1.0 - sign * Y))


def n_perp(X, Y):
    return math.sqrt(1.0 - X), math.sqrt(1.0 - X * (1.0 - X) / (1.0 - X - Y * Y))


X_FIR = (wpe_eff / w_FIR) ** 2
Y_FIR_t, Y_FIR_p = wce_t / w_FIR, wce_p / w_FIR
X_CO2 = (wpe_eff / w_CO2) ** 2
Y_CO2_t, Y_CO2_p = wce_t / w_CO2, wce_p / w_CO2
print('FIR 本征模（B=Bt）: n_R = %.10f, n_L = %.10f' % (n_par(X_FIR, Y_FIR_t, +1), n_par(X_FIR, Y_FIR_t, -1)))
print('CO2 本征模（B=Bt）: n_R = %.12f, n_L = %.12f' % (n_par(X_CO2, Y_CO2_t, +1), n_par(X_CO2, Y_CO2_t, -1)))

# ---------- [AN.2] 法拉第旋转 ----------
coef_FR = e ** 3 / (8 * math.pi ** 2 * c ** 3 * eps0 * me ** 2)   # 旋转角系数（相位差的 1/2）
Int_neB = ne_L_int * Bt
print()
print('== [AN.2] 法拉第旋转 ==')
print('旋转角系数 e^3/(8 pi^2 c^3 eps0 me^2) = %.6g' % coef_FR)
print('与文献 2.62e-13 的比 = %.6g' % (coef_FR / 2.62e-13))
print('Int(ne B_par dl) = %.6g T m^-2' % Int_neB)
for tag, lam, w, Y in (('CO2', lam_CO2, w_CO2, Y_CO2_t), ('FIR', lam_FIR, w_FIR, Y_FIR_t)):
    X = (wpe_eff / w) ** 2
    dn = n_par(X, Y, +1) - n_par(X, Y, -1)
    phase = (w / c) * abs(dn) * L_center            # R/L 相位差
    rot = 0.5 * phase                                # 偏振面旋转角
    rot_an = coef_FR * lam ** 2 * Int_neB
    print('%s: |dn| = %.6g ; 相位差 = %.6g rad = %.6g deg ; 旋转角 = %.6g rad = %.6g deg' % (
        tag, abs(dn), phase, math.degrees(phase), rot, math.degrees(rot)))
    print('     解析旋转角 = %.6g rad = %.6g deg ; 精确/解析 = %.6g' % (
        rot_an, math.degrees(rot_an), rot / rot_an))
print('FIR/CO2 旋转角比 = %.6g' % ((lam_FIR / lam_CO2) ** 2))
rot_CO2 = coef_FR * lam_CO2 ** 2 * Int_neB
rot_FIR = coef_FR * lam_FIR ** 2 * Int_neB
print('FIR 旋转角 = %.6g rad = %.6g deg = %.6g 圈（2pi 模糊）' % (rot_FIR, math.degrees(rot_FIR), rot_FIR / (2 * math.pi)))
print('CO2 旋转角 = %.6g rad = %.6g deg' % (rot_CO2, math.degrees(rot_CO2)))
print('由 CO2 旋转角反推的磁场（已知密度剖面）= %.6g T' % (rot_CO2 / (coef_FR * lam_CO2 ** 2 * ne_L_int)))
print('FIR 每 360deg 对应的 Int(ne B dl) = %.6g T m^-2' % (2 * math.pi / (coef_FR * lam_FIR ** 2)))

# ---------- [AN.3] 科顿-穆顿双折射 ----------
coef_CM = e ** 4 / (16 * math.pi ** 3 * c ** 4 * eps0 * me ** 3)
Int_neB2 = ne_L_int * Bp ** 2
print()
print('== [AN.3] 科顿-穆顿双折射 ==')
print('相位系数 e^4/(16 pi^3 c^4 eps0 me^3) = %.6g（SI：lam m, ne m^-3, B T）' % coef_CM)
print('等效于（B 用高斯）= %.6g' % (coef_CM / 1e8))
print('与文献 2.45e-19（B 用高斯）的比 = %.6g' % (coef_CM / 1e8 / 2.45e-19))
print('Int(ne B_perp^2 dl)（垂直弦, B_perp=Bp）= %.6g T^2 m^-2' % Int_neB2)
for tag, lam, w, Y in (('FIR', lam_FIR, w_FIR, Y_FIR_p), ('CO2', lam_CO2, w_CO2, Y_CO2_p)):
    X = (wpe_eff / w) ** 2
    nO, nX = n_perp(X, Y)
    dn = abs(nX - nO)
    phase = (w / c) * dn * L_center
    phase_an = coef_CM * lam ** 3 * Int_neB2
    print('%s: |dn| = %.6g ; 相位 = %.6g rad = %.6g deg ; 解析 = %.6g rad ; 精确/解析 = %.6g' % (
        tag, dn, phase, math.degrees(phase), phase_an, phase / phase_an))
phi_CM_FIR = coef_CM * lam_FIR ** 3 * Int_neB2
phi_CM_CO2 = coef_CM * lam_CO2 ** 3 * Int_neB2
print('FIR 的 CM = %.6g rad = %.6g deg' % (phi_CM_FIR, math.degrees(phi_CM_FIR)))
print('CO2 的 CM = %.6g rad = %.6g deg' % (phi_CM_CO2, math.degrees(phi_CM_CO2)))
print('FIR: CM/FR 比 = %.6g ; CO2: CM/FR 比 = %.6g' % (phi_CM_FIR / rot_FIR, phi_CM_CO2 / rot_CO2))
print('CM 的 FIR/CO2 比 = %.6g（lam^3）' % ((lam_FIR / lam_CO2) ** 3))
print('切向弦（B_perp=Bt）的 CM 放大 = %.6g 倍；相位 = %.6g deg' % (
    (Bt / Bp) ** 2, math.degrees(phi_CM_FIR * (Bt / Bp) ** 2)))

# ---------- [AN.4] MSE ----------
vb = math.sqrt(2 * Eb * e / mD)
Em = vb * Bt
dE_stark = 9.0 * e * a0 * Em
dlam_stark = lam_Dalpha ** 2 * dE_stark / (h * c)
gamma = math.degrees(math.atan(Bp / Bt))
q_cyl = a * Bt / (R0 * Bp)
dlam_dopp = lam_Dalpha * math.sqrt(2 * Tbeam / (mD * c ** 2 / e))
print()
print('== [AN.4] MSE 运动斯塔克效应 ==')
print('v_beam(100 keV D) = %.6g m/s = %.6g c' % (vb, vb / c))
print('E_motional = %.6g V/m = %.6g MV/m' % (Em, Em / 1e6))
print('最大 Stark 分裂 = %.6g J = %.6g eV = %.6g nm（D_alpha）' % (dE_stark, dE_stark / e, dlam_stark / 1e-9))
print('束内多普勒宽度 = %.6g nm ; Stark/多普勒 = %.6g' % (dlam_dopp / 1e-9, dlam_stark / dlam_dopp))
print('投掷角 gamma = %.6g deg ; q_cyl = %.6g' % (gamma, q_cyl))
print('gamma 变 1 deg 所需 Bp 变化 = %.6g%%' % (100 * (math.tan(math.radians(gamma + 1)) / math.tan(math.radians(gamma)) - 1)))
print('Bp 增 10%% -> gamma = %.6g deg（+%.6g deg）' % (
    math.degrees(math.atan(1.1 * Bp / Bt)), math.degrees(math.atan(1.1 * Bp / Bt)) - gamma))

# ---------- [AN.5] 工程口径 ----------
dphi = math.radians(dphi_res)
neB_res = dphi / (coef_FR * lam_CO2 ** 2)
print()
print('== [AN.5] 工程口径 ==')
print('相位分辨率 %.6g deg = %.6g rad -> Int(ne B dl) 分辨率 = %.6g T m^-2（相对 %.6g%%）' % (
    dphi_res, dphi, neB_res, 100 * neB_res / Int_neB))
print('等效平均磁场分辨率 = %.6g T（%.6g mT）' % (neB_res / ne_L_int, 1e3 * neB_res / ne_L_int))
print('PEM 周期 = %.6g us ; 1 ms 内调制周期数 = %.6g' % (1e6 / f_PEM, 1e-3 * f_PEM))
print('9 弦弦距 = %.6g m ; 中心弦 L = %.6g m ; 边缘弦(r=0.9a) L = %.6g m' % (
    2 * a / (N_chord - 1), L_center, 2 * math.sqrt(a ** 2 - (0.9 * a) ** 2)))
print('边缘/中心弦的 Int(ne dl) 比 = %.6g' % (math.sqrt(1 - 0.9 ** 2)))
print('偏振消光比 1e-3 对应的相位误差 = %.6g deg' % math.degrees(math.sqrt(1e-3)))

# ---------- [AN.6] 反演与误差传播 ----------
print()
print('== [AN.6] 反演与误差传播 ==')
print('干涉仪 0.1%% 密度误差 + FR 0.42%% 联合 -> B 的相对误差量级 = 0.42%%')
print('q 的相对灵敏度 dq/q / dgamma = %.6g /rad = %.6g /deg（q=%.6g）' % (
    1 / math.tan(math.radians(gamma)), (1 / math.tan(math.radians(gamma))) * math.pi / 180, q_cyl))
print('1 deg 投掷角误差 -> q 相对误差 = %.6g%%' % (100 * (1 / math.tan(math.radians(gamma))) * math.pi / 180))
print('q=3 所需 Bp = %.6g T（相对当前 %.6g 倍）' % (a * Bt / (R0 * 3.0), (a * Bt / (R0 * 3.0)) / Bp))

# ---------- [AN.7] 练习解答 ----------
print()
print('== [AN.7] 练习解答 ==')
print('1) FIR/CO2 旋转角比 = %.6g ; FIR 旋转角 = %.6g deg = %.6g 圈' % (
    (lam_FIR / lam_CO2) ** 2, math.degrees(rot_FIR), rot_FIR / (2 * math.pi)))
rot_half = coef_FR * lam_CO2 ** 2 * (0.5 * Int_neB)
print('2) ne0 减半 -> Int(ne B dl) = %.6g T m^-2 ; CO2 旋转角 = %.6g deg' % (0.5 * Int_neB, math.degrees(rot_half)))
print('3) Bp 加倍 -> CM(FIR) = %.6g deg（4 倍）' % math.degrees(4 * phi_CM_FIR))
vb2 = math.sqrt(2 * 150.0e3 * e / mD)
Em2 = vb2 * Bt
dE2 = 9.0 * e * a0 * Em2
dlam2 = lam_Dalpha ** 2 * dE2 / (h * c)
print('4) 150 keV -> v = %.6g m/s ; E_m = %.6g MV/m ; 分裂 = %.6g nm ; 比 100 keV 大 %.6g 倍' % (
    vb2, Em2 / 1e6, dlam2 / 1e-9, dlam2 / dlam_stark))
d_ne = ne0 * 4 * a / 3 - ne0 * 4 * a / 3 * (1 - (2 * a / 12 / a) ** 2) ** 1.5
print('5) 13 弦：弦距 = %.6g m ; 相邻弦 Int(ne dl) 差 = %.6g m^-2' % (2 * a / 12, d_ne))
print('6) 0.005 deg -> Int(ne B dl) 分辨率 = %.6g T m^-2（相对 %.6g%%）' % (
    math.radians(0.005) / (coef_FR * lam_CO2 ** 2), 100 * (math.radians(0.005) / (coef_FR * lam_CO2 ** 2)) / Int_neB))
