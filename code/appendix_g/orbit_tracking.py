import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# -------------------- 物理参数 --------------------
q = 1.602e-19       # 电荷 (C)
m = 1.673e-27       # 质量 (kg, 质子)
E0 = np.array([0.0, 1e3, 0.0])   # 均匀电场 (V/m)，沿 y
B0 = 1.0                          # 均匀磁场幅值 (T)，沿 z
L  = 0.05                         # 磁镜梯度长度 (m)，用于非均匀场

def B_uniform(r):
    """均匀场 B = B0 * z_hat"""
    return np.array([0.0, 0.0, B0])

def B_mirror(r):
    """近轴磁镜场，满足 div B = 0:
       B_z = B0(1 + (z/L)^2),  B_r = -(r/2) dB_z/dz = -r*B0*z/L^2
       (r 为到轴的垂直距离)"""
    x, y, z = r
    Bz = B0 * (1.0 + (z/L)**2)
    Br_scale = -0.5 * B0 * (2.0*z/L**2)   # (1/2) dBz/dz 的负号
    return np.array([Br_scale*x, Br_scale*y, Bz])

# -------------------- 洛伦兹力方程 --------------------
def make_rhs(B_func):
    def rhs(t, y):
        r, v = y[0:3], y[3:6]
        a = (q/m) * np.cross(v, B_func(r))     # 仅磁力项 (E=0)
        return np.concatenate([v, a])
    return rhs

# -------------------- 情形 A: 均匀场 + E x B 漂移 --------------------
def rhs_EB(t, y):
    r, v = y[0:3], y[3:6]
    a = (q/m) * (E0 + np.cross(v, np.array([0, 0, B0])))
    return np.concatenate([v, a])

y0 = np.concatenate([[0, 0, 0], [1e5, 0, 1e4]])
t = np.linspace(0, 2e-6, 2000)
solA = solve_ivp(rhs_EB, [t[0], t[-1]], y0, t_eval=t, rtol=1e-9, atol=1e-12)
vdrift = np.cross(E0, np.array([0, 0, B0])) / B0**2
print("E x B drift speed [m/s]:", np.linalg.norm(vdrift))

# -------------------- 情形 B: 磁镜反射 --------------------
y0m = np.concatenate([[0.005, 0, 0], [0, 3e5, 1.0e5]])
t2 = np.linspace(0, 3e-6, 4000)
solB = solve_ivp(make_rhs(B_mirror), [t2[0], t2[-1]], y0m,
                 t_eval=t2, rtol=1e-9, atol=1e-12)
zB = solB.y[2]
print("mirror turning points z [m]:", zB.min(), zB.max())

# -------------------- 能量守恒检查 --------------------
def kinetic(sol):
    v2 = sol.y[3]**2 + sol.y[4]**2 + sol.y[5]**2
    return 0.5*m*v2
KE = kinetic(solB)
print("relative KE drift:", (KE.max()-KE.min())/KE.mean())

# -------------------- 绘图 --------------------
fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot(solA.y[0], solA.y[1], solA.y[2], lw=0.8)
ax1.set_title("Uniform B: orbit + ExB"); ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]'); ax1.set_zlabel('z [m]')
ax2 = fig.add_subplot(132)
ax2.plot(solA.y[0], solA.y[1], lw=0.8)
ax2.set_aspect('equal'); ax2.set_title("x-y projection (ExB drift)")
ax2.set_xlabel('x [m]'); ax2.set_ylabel('y [m]')
ax3 = fig.add_subplot(133)
ax3.plot(solB.t*1e6, zB, lw=0.8)
ax3.set_title("Mirror: z(t)"); ax3.set_xlabel('t [us]'); ax3.set_ylabel('z [m]')
plt.tight_layout()
plt.savefig('single_particle_orbit.png', dpi=150)
plt.show()
