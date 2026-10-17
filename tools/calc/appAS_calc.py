# -*- coding: utf-8 -*-
"""附录 AS 数值基线：射线追踪与波束光学（wajue 工作点 D1）。

只 print，不写文件；输出重定向到 tmp/appAS_calc.txt 作为正文与门禁的唯一数值来源。
工作点（与其他附录同口径）：n_e0=1e20 m^-3, T_e=9.194 keV, B0=5.3 T, R0=6.2 m, a=2.0 m。
"""
import math

# ---- 声明输入 ----
n_e0 = 1.0e20          # m^-3 中心密度
T_e_keV = 9.194        # keV
B0 = 5.3               # T 轴上环向场
R0 = 6.2               # m 大半径
a = 2.0                # m 小半径（限制器半径）
f_EC = 140.0e9         # Hz 电子回旋加热频率
f_ref_lo = 50.0e9      # Hz 反射计扫频下限
f_ref_mid = 60.0e9
f_ref_hi = 75.0e9
f_ECE = 280.0e9        # Hz ECE 通道（二次谐波）
df_ECE = 1.0e9         # Hz 滤波器带宽
w0 = 0.02              # m 发射腰斑
f_mirror = 1.0                # m 聚焦镜焦距
D_mirror = 0.20        # m 镜面口径
P_EC = 1.0e6           # W 单束功率
theta_launch = 20.0    # deg 发射方向与 B 的夹角
theta_launch2 = 60.0   # deg 大角度发射（用于找折射转折点）
v_perp = 1.0e4         # m/s 湍流垂直速度（多普勒反射计）
theta_D = 20.0         # deg 多普勒反射计夹角
eta_abs = 0.95         # 声明：基频 O 模单程吸收份额

c = 2.99792458e8
e = 1.602176634e-19
m_e = 9.1093837015e-31
eps0 = 8.8541878128e-12
keV = 1.0e3 * e

T_e = T_e_keV * keV                       # J

print('==== 声明输入 ====')
print('n_e0 = %.6g m^-3, T_e = %.6g keV = %.6g J, B0 = %.6g T, R0 = %.6g m, a = %.6g m'
      % (n_e0, T_e_keV, T_e, B0, R0, a))
print('f_EC = %.6g Hz, 反射计 %.6g/%.6g/%.6g Hz, ECE %.6g Hz (df=%.6g), P_EC=%.6g W'
      % (f_EC, f_ref_lo, f_ref_mid, f_ref_hi, f_ECE, df_ECE, P_EC))
print('w0 = %.6g m, 镜 f=%.6g m / D=%.6g m, theta=%.6g deg, eta_abs=%.6g'
      % (w0, f_mirror, D_mirror, theta_launch, eta_abs))

# ================= AS.1 从波动到射线：WKB 判据 =================
print()
print('==== AS.1 从波动到射线：WKB 判据 ====')
omega_EC = 2 * math.pi * f_EC
n_c_EC = eps0 * m_e * omega_EC ** 2 / e ** 2
X0 = n_e0 / n_c_EC
k0 = omega_EC / c
lam0 = c / f_EC
n_p = math.sqrt(1 - X0)
k_p = k0 * n_p
lam_p = 2 * math.pi / k_p
print('omega_EC        = %.10g rad/s' % omega_EC)
print('n_c(140GHz)     = %.10g m^-3' % n_c_EC)
print('X0 = n_e0/n_c   = %.10g' % X0)
print('k0 = omega/c    = %.10g m^-1，lambda0 = %.10g m = %.6g mm' % (k0, lam0, lam0 * 1e3))
print('n_p = sqrt(1-X0)= %.10g，k_p = %.10g m^-1，lambda_p = %.10g mm'
      % (n_p, k_p, lam_p * 1e3))
print('k_p * a         = %.10g  （WKB 判据 kL >> 1）' % (k_p * a))
print('lambda_p/(2 pi a) = %.10g' % (lam_p / (2 * math.pi * a)))

n_e = lambda r: n_e0 * (1 - (r / a) ** 2)
Xr = lambda r: n_e(r) / n_c_EC
nr = lambda r: math.sqrt(max(1 - Xr(r), 1e-12))
kn = lambda r: k0 * nr(r)
# 局部密度梯度尺度 L_n = n_e / |dn_e/dr|
L_n = lambda r: (n_e0 * (1 - (r / a) ** 2)) / (2 * n_e0 * r / a ** 2) if r > 0 else float('inf')
print('L_n(0.2a)=%.10g m, L_n(0.5a)=%.10g m, L_n(0.9a)=%.10g m'
      % (L_n(0.2 * a), L_n(0.5 * a), L_n(0.9 * a)))
kLn_min = min(kn(r) * L_n(r) for r in [0.1 * a + 0.01 * a * i for i in range(90)])
r_min = min([0.1 * a + 0.01 * a * i for i in range(90)], key=lambda r: kn(r) * L_n(r))
print('kL_n 最小 = %.10g （在 r = %.6g m = %.6g a 处）' % (kLn_min, r_min, r_min / a))
print('1/(kL_n)_min = %.10g' % (1.0 / kLn_min))
print('1/(kL_n)_min 百分数 = %.10g %%' % (100.0 / kLn_min))

# ================= AS.2 射线方程与折射 =================
print()
print('==== AS.2 射线方程与折射 ====')
s20 = math.sin(math.radians(theta_launch))
n_turn_20 = n_c_EC * (1 - s20 ** 2)
print('sin(theta=%.6g deg) = %.10g' % (theta_launch, s20))
print('转折密度 n_e,turn(20 deg) = n_c*(1-sin^2) = %.10g m^-3 （> n_e0 -> 不转折）' % n_turn_20)
print('n_turn/n_e0 = %.10g' % (n_turn_20 / n_e0))

s60 = math.sin(math.radians(theta_launch2))
n_turn_60 = n_c_EC * (1 - s60 ** 2)
r_t60 = a * math.sqrt(1 - n_turn_60 / n_e0)
print('转折密度 n_e,turn(60 deg) = %.10g m^-3' % n_turn_60)
print('r_t(60 deg) = %.10g m = %.6g a （转折点）' % (r_t60, r_t60 / a))
print('转折点处 n(r_t) = %.10g，k(r_t) = %.10g m^-1' % (nr(r_t60), kn(r_t60)))
print('转折点处 L_n = %.10g m，kL_n = %.10g' % (L_n(r_t60), kn(r_t60) * L_n(r_t60)))

# 射线路径长度（数值积分 ds = dr / cos(theta(r))，theta 由 n sin(theta)=sin(theta0) 定）
def path_length(th_deg, r_start, r_end, npts=200000):
    s0 = math.sin(math.radians(th_deg))
    total = 0.0
    dr = (r_start - r_end) / npts
    for i in range(npts):
        r = r_start - (i + 0.5) * dr
        s = s0 / nr(r)
        s = min(s, 1.0)
        cs = math.sqrt(max(1 - s * s, 1e-12))
        total += dr / cs
    return total

L_path60 = path_length(theta_launch2, a, r_t60)
print('射线路径长度（a -> r_t, 60 deg） = %.10g m' % L_path60)
print('直线径向距离 = %.10g m，路径/直线 = %.10g' % (a - r_t60, L_path60 / (a - r_t60)))
# 单程相位
def phase_along(th_deg, r_start, r_end, npts=200000):
    s0 = math.sin(math.radians(th_deg))
    total = 0.0
    dr = (r_start - r_end) / npts
    for i in range(npts):
        r = r_start - (i + 0.5) * dr
        s = s0 / nr(r)
        s = min(s, 1.0)
        cs = math.sqrt(max(1 - s * s, 1e-12))
        total += kn(r) * dr / cs
    return total

phi60 = phase_along(theta_launch2, a, r_t60)
print('单程相位积分 = %.10g rad = %.10g 圈' % (phi60, phi60 / (2 * math.pi)))

# 小角度入射时的折射角（在 r=0 处）
for th in (10.0, 20.0, 30.0):
    s0 = math.sin(math.radians(th))
    s_in = s0 / nr(0.0)
    print('发射 %.6g deg -> 中心处 sin(theta)=%.10g，theta=%.6g deg'
          % (th, s_in, math.degrees(math.asin(min(s_in, 1.0)))))

# ================= AS.3 高斯波束与准光学 =================
print()
print('==== AS.3 高斯波束与准光学 ====')
z_R = math.pi * w0 ** 2 / lam0
theta_div = lam0 / (math.pi * w0)
print('z_R = pi w0^2/lambda = %.10g m' % z_R)
print('发散角 theta = lambda/(pi w0) = %.10g rad = %.6g deg' % (theta_div, math.degrees(theta_div)))
wz = lambda z: w0 * math.sqrt(1 + (z / z_R) ** 2)
for z in (1.0, 2.0, 4.0):
    print('w(%.6g m) = %.10g m' % (z, wz(z)))
Rz = lambda z: z * (1 + (z_R / z) ** 2) if z else float('inf')
print('波前曲率半径 R(2 m) = %.10g m' % Rz(2.0))
I_waist = 2 * P_EC / (math.pi * w0 ** 2)
print('腰斑处峰值强度 I0 = 2P/(pi w0^2) = %.10g W/m^2' % I_waist)
I_center = 2 * P_EC / (math.pi * wz(2.0) ** 2)
print('等离子体中心处 I = %.10g W/m^2（w=%.6g m）' % (I_center, wz(2.0)))
w_f = lam0 * f_mirror / (math.pi * (D_mirror / 2))
z_Rf = math.pi * w_f ** 2 / lam0
print('聚焦腰斑 w_f = lambda f/(pi D/2) = %.10g m = %.6g mm' % (w_f, w_f * 1e3))
print('聚焦后瑞利长度 z_Rf = %.10g m' % z_Rf)
I_f = 2 * P_EC / (math.pi * w_f ** 2)
E_f = math.sqrt(2 * I_f / (c * eps0))
print('焦斑峰值强度 I_f = %.10g W/m^2' % I_f)
print('焦斑电场 E_f = sqrt(2I/(c eps0)) = %.10g V/m = %.6g MV/m' % (E_f, E_f / 1e6))
print('焦斑面积 pi w_f^2 = %.10g m^2 = %.6g mm^2'
      % (math.pi * w_f ** 2, math.pi * w_f ** 2 * 1e6))
E_center = math.sqrt(2 * I_center / (c * eps0))
print('中心处电场 E = %.10g V/m = %.6g kV/m' % (E_center, E_center / 1e3))
print('焦斑/中心电场比 E_f/E_center = %.10g' % (E_f / E_center))

# ================= AS.4 反射计与密度剖面 =================
print()
print('==== AS.4 反射计与密度剖面 ====')
n_c_of = lambda f: eps0 * m_e * (2 * math.pi * f) ** 2 / e ** 2
n_c_prefactor = eps0 * m_e * (2 * math.pi * 10e9) ** 2 / e ** 2  # f 以 GHz 计
print('n_c 系数（f 单位 GHz）= %.10g m^-3/GHz^2' % n_c_prefactor)
r_t_of = lambda f: a * math.sqrt(max(1 - n_c_of(f) / n_e0, 0.0))
for f in (f_ref_lo, f_ref_mid, f_ref_hi):
    print('f=%.6g GHz: n_c=%.10g m^-3, r_t=%.10g m (%.6g a)'
          % (f / 1e9, n_c_of(f), r_t_of(f), r_t_of(f) / a))

def group_delay(f, npts=200000):
    """往返群延迟 tau = 2 * int_{r_t}^{a} dr / (c n(r))，O 模 n=sqrt(1-X)。"""
    rt = r_t_of(f)
    nc = n_c_of(f)
    def nn(r):
        return math.sqrt(max(1 - n_e(r) / nc, 1e-12))
    dr = (a - rt) / npts
    total = 0.0
    for i in range(npts):
        r = rt + (i + 0.5) * dr
        total += dr / nn(r)
    return 2 * total / c

for f in (f_ref_lo, f_ref_mid, f_ref_hi):
    tau = group_delay(f)
    print('f=%.6g GHz: 往返群延迟 tau = %.10g s = %.6g ns' % (f / 1e9, tau, tau * 1e9))
tau_lo, tau_mid, tau_hi = group_delay(f_ref_lo), group_delay(f_ref_mid), group_delay(f_ref_hi)
dtau = tau_hi - tau_lo
dr = r_t_of(f_ref_lo) - r_t_of(f_ref_hi)
print('50->75 GHz: Delta r = %.10g m, Delta tau = %.10g s = %.6g ns'
      % (dr, dtau, dtau * 1e9))
print('平均灵敏度 dr/df = %.10g m/GHz' % (dr / 25.0e9 * 1e9))
# 相位 -> 密度分辨率（dphi/df = 2 pi tau(f)，用最高频的群延迟）
dphi = 0.01  # rad（声明相位分辨）
df_from_phi = dphi / (2 * math.pi * tau_hi)
dr_from_phi = abs(dr / 25.0e9) * df_from_phi
dn_from_phi = (2 * n_e0 * r_t_of(f_ref_hi) / a ** 2) * dr_from_phi
print('相位分辨 %.6g rad 折合频率分辨 df = dphi/(2 pi tau) = %.10g Hz = %.6g kHz'
      % (dphi, df_from_phi, df_from_phi / 1e3))
print('折合半径分辨 dr = |dr/df| df = %.10g m = %.6g um' % (dr_from_phi, dr_from_phi * 1e6))
print('折合密度分辨 dn_e = |dn/dr| dr = %.10g m^-3（在 r_t(75GHz) 处，相对 %.6g）'
      % (dn_from_phi, dn_from_phi / n_e(r_t_of(f_ref_hi))))
f_D = 2 * v_perp * math.sin(math.radians(theta_D)) / lam0
print('多普勒反射计频移 f_D = 2 v_perp sin(theta)/lambda = %.10g Hz = %.6g MHz'
      % (f_D, f_D / 1e6))
df_D = 1.0e5  # Hz 声明多普勒谱频率分辨
dv_D = lam0 * df_D / (2 * math.sin(math.radians(theta_D)))
print('多普勒谱分辨 %.6g kHz 折合速度分辨 dv = lambda df/(2 sin) = %.10g m/s'
      % (df_D / 1e3, dv_D))

# ================= AS.5 电子回旋发射成像 =================
print()
print('==== AS.5 电子回旋发射成像 ====')
f_ce_of_B = lambda B: e * B / (2 * math.pi * m_e)
f_ce_coef = e / (2 * math.pi * m_e)
print('f_ce 系数 = %.10g Hz/T（= %.6g GHz/T）' % (f_ce_coef, f_ce_coef / 1e9))
B_of_r = lambda r: B0 * R0 / (R0 + r)
f_ce_r = lambda r: f_ce_of_B(B_of_r(r))
for r in (0.0, 1.0, 2.0):
    print('r=%.6g m: B=%.10g T, f_ce=%.10g GHz, 2f_ce=%.10g GHz'
          % (r, B_of_r(r), f_ce_r(r) / 1e9, 2 * f_ce_r(r) / 1e9))
B_ECE = f_ECE / 2 / (e / (2 * math.pi * m_e))
r_ECE = B0 * R0 / B_ECE - R0
print('ECE %.6g GHz -> B=%.10g T -> r=%.10g m' % (f_ECE / 1e9, B_ECE, r_ECE))
dr_df = 2 * (e / (2 * math.pi * m_e)) * B0 * R0 / f_ECE ** 2
print('|dr/df| 在 %.6g GHz = %.10g m/Hz -> 带宽 %.6g GHz 对应 dr = %.10g m = %.6g cm'
      % (f_ECE / 1e9, dr_df, df_ECE / 1e9, dr_df * df_ECE, dr_df * df_ECE * 100))
B_nu = 2 * f_ECE ** 2 * T_e / c ** 2
print('瑞利-金斯亮度 B_nu = 2 nu^2 k T/c^2 = %.10g W m^-2 Hz^-1 sr^-1' % B_nu)
lam_ECE = c / f_ECE
print('ECE 波长 lambda = %.10g m = %.6g mm' % (lam_ECE, lam_ECE * 1e3))
P_ECE = B_nu * df_ECE * lam_ECE ** 2
print('单模接收功率 P = B_nu df lambda^2 = %.10g W = %.6g uW' % (P_ECE, P_ECE * 1e6))
print('A Omega 截止：A=1e-2 m^2 时 Omega = %.10g sr' % (lam_ECE ** 2 / 1.0e-2))

# ================= AS.6 功率沉积定位 =================
print()
print('==== AS.6 功率沉积定位 ====')
f_ce_of_r_res = lambda r: f_EC
r_res = B0 * R0 / (f_EC / (e / (2 * math.pi * m_e))) - R0
print('基频共振 r_res(f_ce=f_EC) = %.10g m' % r_res)
n_par = n_p * math.cos(math.radians(theta_launch))
v_the = math.sqrt(2 * T_e / m_e)
print('热速度 v_the = sqrt(2T_e/m_e) = %.10g m/s' % v_the)
print('热速度（本书其他附录口径）sqrt(T_e/m_e) = %.10g m/s' % math.sqrt(T_e / m_e))
dv_ratio = n_par * v_the / c
print('n_par = n_p cos(theta) = %.10g' % n_par)
print('多普勒展宽 Delta omega/omega = n_par v_the/c = %.10g' % dv_ratio)
for sign, tag in ((+1, '+'), (-1, '-')):
    w_eff = 2 * math.pi * f_EC * (1 - sign * dv_ratio)
    B_shift = w_eff / (e / m_e)
    r_shift = B0 * R0 / B_shift - R0
    print('多普勒 %s 支：omega_eff=%.10g rad/s -> B=%.10g T -> r=%.10g m（偏移 %.6g m）'
          % (tag, w_eff, B_shift, r_shift, r_shift - r_res))
r_plus = B0 * R0 / ((2 * math.pi * f_EC * (1 - dv_ratio)) / (e / m_e)) - R0
print('1/e 半宽（+ 支）Delta r = %.10g m' % abs(r_plus - r_res))
print('共振层厚度（把两支都算上，- 支已在等离子体外）= %.10g m'
      % (2 * abs(r_plus - r_res)))
P_abs = eta_abs * P_EC
print('吸收功率 P_abs = eta_abs P = %.10g W = %.6g MW' % (P_abs, P_abs / 1e6))
w_res = wz(a - r_res)
print('共振处波束半径 w = %.10g m' % w_res)
A_beam = math.pi * w_res ** 2
V_dep = A_beam * abs(2 * (r_shift - r_res))
print('波束截面 = %.10g m^2，沉积体积 = %.10g m^3' % (A_beam, V_dep))
print('平均功率密度 = %.10g W/m^3' % (P_abs / V_dep))
print('功率密度上限（把功率压在一个焦斑体积） = %.10g W/m^3' % (P_abs / (A_beam * 2 * w_f)))
I_dep = P_abs / A_beam
print('共振处功率面密度 = %.10g W/m^2' % I_dep)
print('单程吸收份额 %.6g -> 双程 %.10g' % (eta_abs, 1 - (1 - eta_abs) ** 2))

# ================= AS.7 估算练习 =================
print()
print('==== AS.7 估算练习解答 ====')
# 1) 100 GHz O 模垂直入射，n_e=5e19，求折射率与相速度、k
f1 = 100.0e9
nc1 = n_c_of(f1)
X1 = 5.0e19 / nc1
n1 = math.sqrt(1 - X1)
k1 = 2 * math.pi * f1 / c * n1
vp1 = c / n1
print('[练1] 100 GHz: n_c=%.10g, X=%.10g, n=%.10g, k=%.10g m^-1, v_p=%.10g m/s (%.6g c)'
      % (nc1, X1, n1, k1, vp1, vp1 / c))
# 2) 腰斑减半 -> 瑞利长度与发散角
w_half = w0 / 2
z_R_h = math.pi * w_half ** 2 / lam0
print('[练2] w0/2=%.6g m -> z_R=%.10g m（原 %.6g m），theta=%.10g deg（原 %.6g deg）'
      % (w_half, z_R_h, z_R, math.degrees(lam0 / (math.pi * w_half)),
         math.degrees(theta_div)))
# 3) 70 GHz 反射计的转折半径与群延迟
f3 = 70.0e9
print('[练3] 70 GHz: n_c=%.10g, r_t=%.10g m, tau=%.10g s = %.6g ns'
      % (n_c_of(f3), r_t_of(f3), group_delay(f3), group_delay(f3) * 1e9))
# 4) 多普勒频移 -> 垂直速度
f_D4 = 5.0e6
v4 = f_D4 * lam0 / (2 * math.sin(math.radians(theta_D)))
print('[练4] f_D=%.6g MHz -> v_perp=%.10g m/s = %.6g km/s' % (f_D4 / 1e6, v4, v4 / 1e3))
# 5) ECE 260 GHz -> r 与分辨率
f5 = 260.0e9
B5 = f5 / 2 / (e / (2 * math.pi * m_e))
r5 = B0 * R0 / B5 - R0
dr5 = 2 * (e / (2 * math.pi * m_e)) * B0 * R0 / f5 ** 2 * df_ECE
print('[练5] 260 GHz -> B=%.10g T, r=%.10g m, df=%.6g GHz -> dr=%.10g m = %.6g cm'
      % (B5, r5, df_ECE / 1e9, dr5, dr5 * 100))
# 6) 2 MW、w_f=5 mm 的焦斑强度与电场
P6, w6 = 2.0e6, 5.0e-3
I6 = 2 * P6 / (math.pi * w6 ** 2)
E6 = math.sqrt(2 * I6 / (c * eps0))
print('[练6] 2 MW, w_f=5 mm: I=%.10g W/m^2, E=%.10g V/m = %.6g MV/m' % (I6, E6, E6 / 1e6))
print('[练6] 面积比 (5mm/w_f)^2 = %.10g，强度比 I6/I_f = %.10g'
      % ((w6 / w_f) ** 2, I6 / I_f))
print('[练6] 声明输入 P=%.6g W, w_f=%.6g m' % (P6, w6))

print()
print('==== 声明输入补充 ====')
print('练1 工作密度 5.0e19 m^-3；练4 多普勒频移 5.0e6 Hz；练5 ECE 260 GHz')
print('相位分辨 0.01 rad；多普勒谱分辨 1.0e5 Hz；接收面积口径 1.0e-2 m^2')
