import numpy as np
import matplotlib.pyplot as plt

# -------------------- 物理参数 --------------------
n    = 1e18            # m^-3
Te_eV = 10.0           # eV
q_e  = 1.602e-19
m_e  = 9.109e-31
eps0 = 8.854e-12
Theta = Te_eV * q_e    # J

# 等离子体频率与德拜长度
omega_pe = np.sqrt(n * q_e**2 / (eps0 * m_e))
v_te     = np.sqrt(Theta / m_e)          # 注意: 本书 v_te = sqrt(k_B T_e/m_e)
lam_D    = v_te / omega_pe
print("omega_pe [rad/s]:", omega_pe)
print("lambda_D [m]:", lam_D)

# -------------------- Bohm-Gross 色散关系 --------------------
# omega^2 = omega_pe^2 + 3 k^2 v_te^2 ;  近似有效范围 k*lambda_D <~ 0.3
k = np.linspace(0.0, 1.0 / lam_D, 500)
omega = np.sqrt(omega_pe**2 + 3.0 * k**2 * v_te**2)

# -------------------- 有效性回溯检查 --------------------
# 与冷等离子体对比: 热修正量 omega^2/omega_pe^2 - 1 = 3(k lambda_D)^2
kLD = k * lam_D
print("at k*lambda_D=0.1, thermal correction =", 3*0.1**2)
print("at k*lambda_D=0.3, thermal correction =", 3*0.3**2,
      " (Landau damping already noticeable)")
print("at k*lambda_D=1.0, thermal correction =", 3*1.0**2,
      " (Bohm-Gross no longer reliable)")

# -------------------- 绘图 --------------------
plt.figure(figsize=(8, 5))
plt.plot(kLD, omega / omega_pe, 'b-', lw=2, label='Bohm-Gross')
plt.axhline(y=1.0, color='r', linestyle='--', label='Cold: omega = omega_pe')
plt.axvspan(0, 0.3, color='green', alpha=0.08)
plt.text(0.32, 1.02, 'reliable region (k*lambda_D < 0.3)', fontsize=9)
plt.xlabel('k * lambda_D')
plt.ylabel('omega / omega_pe')
plt.title('Langmuir Wave Dispersion Relation')
plt.legend(); plt.grid(True, ls='--', alpha=0.5)
plt.xlim(0, 1.0); plt.ylim(0.9, 2.5)
plt.tight_layout()
plt.savefig('langmuir_dispersion.png', dpi=150)
plt.show()
