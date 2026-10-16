import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# -------------------- 物理参数 --------------------
n0    = 1e18            # 背景密度 m^-3
Te_eV = 10.0            # 电子温度 (eV)
q_e   = 1.602e-19
eps0  = 8.854e-12
Theta = Te_eV * q_e     # 电子热能 (J) —— 必须转成焦耳！
lam_D = np.sqrt(eps0 * Theta / (n0 * q_e**2))
print("lambda_D [m] =", lam_D)

# ------------- 无量纲化: x = r/lam_D, y = e*phi/Theta -------------
# 球坐标泊松方程化为:  y'' + (2/x) y' = exp(y) - 1
x_min, x_max = 1e-3, 15.0      # 计算区间（以 lam_D 为单位）
def poisson_debye(x, y):
    """y = [y, y']，无量纲球坐标泊松方程"""
    y0, y1 = y
    expo = np.clip(y0, -60.0, 60.0)              # 防指数溢出
    return np.vstack([y1, np.exp(expo) - 1.0 - 2.0 * y1 / x])

def bc(ya, yb, amp):
    """边界条件恰好 2 个: y(x_min) = amp, y(x_max) = 0"""
    return np.array([ya[0] - amp, yb[0]])

def solve(amp, tol=1e-8, max_nodes=50000):
    def bc_amp(ya, yb):
        return bc(ya, yb, amp)
    x = np.linspace(x_min, x_max, 400)
    y_guess = np.zeros((2, x.size))
    y_guess[0] = amp * (x_min / x)               # 初始猜测: 1/x 型
    return solve_bvp(poisson_debye, bc_amp, x, y_guess, tol=tol, max_nodes=max_nodes)

# -------------------- 1) 小振幅（线性极限）求解 --------------------
sol = solve(amp=1e-2)
print("success:", sol.success)   # 唯一可信的收敛判据（布尔量）
print("status :", sol.status)    # 0=收敛, 1=网格节点超限, 2=雅可比奇异
print("message:", sol.message)   # 终止原因的文字说明
print("nodes  :", sol.x.size, "| niter:", sol.niter)
print("rms residual max:", np.max(sol.rms_residuals), "| tol = 1e-8")
print("BC residual max :", np.max(np.abs(bc(sol.y[:, 0], sol.y[:, -1], 1e-2))))

# -------------------- 2) 与线性德拜解逐点对比 --------------------
xx = np.linspace(x_min, x_max, 2000)
y_num = sol.sol(xx)[0]
y_lin = 1e-2 * (x_min / xx) * np.exp(-(xx - x_min))   # 同一 x_min 定标
for xp in (0.1, 1.0, 5.0):
    i = np.argmin(np.abs(xx - xp))
    print("x=%.0e: y_num=%.6e  y_lin=%.6e  相对偏差=%.3e"
          % (xp, y_num[i], y_lin[i], y_num[i] / y_lin[i] - 1.0))

# ------------- 3) 线性极限的证据: 相对偏差与 amp 成正比 -------------
for a in (1e-2, 1e-1, 1.0):
    s = solve(a)
    yv = float(s.sol(np.array([1.0]))[0][0])
    yl = a * x_min * np.exp(-(1.0 - x_min))
    dev = yv / yl - 1.0
    print("amp=%.0e: x=1 相对偏差=%+.3e, 偏差/amp=%+.3e (success=%s)"
          % (a, dev, dev / a, s.success))

# ------------- 4) 负例: 故意把节点上限压到 200 -------------
bad = solve(1e-2, tol=1e-12, max_nodes=200)
print("bad: success=%s status=%d message=%s" % (bad.success, bad.status, bad.message))

# -------------------- 5) 绘图 --------------------
plt.figure(figsize=(8, 5))
plt.loglog(xx, np.abs(y_num), lw=2, label='nonlinear Poisson (amp=1e-2)')
plt.loglog(xx, np.abs(y_lin), '--', lw=2, label=r'linear Debye: $e^{-x}/x$')
plt.xlabel('x = r / lambda_D')
plt.ylabel('|y| = |e phi / Theta|')
plt.legend(); plt.grid(True, which='both', ls='--', alpha=0.5)
plt.tight_layout(); plt.savefig('debye_shielding.png', dpi=150)
plt.show()
