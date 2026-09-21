# -*- coding: utf-8 -*-
import sympy as sp
def h(s): print("\n== " + s + " ==")

# ch4 suppl L290 boxed formula: simplify the symbolic result
k,m_i,kB,T_e,T_i,eps0,e,n0,g_ = sp.symbols('k m_i k_B T_e T_i epsilon_0 e n_0 gamma_i', positive=True, real=True)
w2 = k**2*kB*(T_e*T_i*eps0*g_*k**2*kB + T_e*e**2*n0 + T_i*e**2*g_*n0)/(m_i*(T_e*eps0*k**2*kB + e**2*n0))
# divide numerator & denominator by e^2 n0 -> denominator = m_i T_e eps0 k^2 kB/(e^2 n0) + m_i
# note lamDe^2 = eps0 kB T_e/(n0 e^2)
lam = sp.sqrt(eps0*kB*T_e/(n0*e**2))
target = k**2*(g_*kB*T_i/m_i + (kB*T_e/m_i)/(1+k**2*lam**2))
d = sp.simplify(w2 - target)
print("w2 - target =", d)

# ch3 suppl tokamak drift combined formula
h("ch3 suppl L113-129: v = m/(qBR)(vpar^2 + vperp^2/2) from grad-B + curvature")
# v_gradB = (m vperp^2/(2 q B^3)) B x gradB,  |gradB| = B/R, Rhat x B  -> need direction/order
# v_curv  = (m vpar^2/(q B^2 R^2)) B x Rc
# book writes combined as m/(qB^2R) (vpar^2 + vperp^2/2) Rhat x B  and says |Rhat x B| = B
# so |v| = m/(qBR)(vpar^2 + vperp^2/2)
print("|v_gradB| = m vperp^2 |B x gradB|/(2 q B^3) = m vperp^2 (B * B/R)/(2 q B^3) = m vperp^2/(2 q B R) OK")
print("|v_curv|  = m vpar^2 |B x Rc|/(q B^2 R^2)  = m vpar^2 (B R)/(q B^2 R^2) = m vpar^2/(q B R) OK")
print("so |v| = m/(qBR)(vpar^2 + vperp^2/2) : CORRECT")

# but the dimensions of the book's vector formula
R,B,vpar2,vperp2,m,q = sp.symbols('R B v_par2 v_perp2 m q', positive=True)
print("book vector form: m/(q B^2 R)(vpar^2+vperp^2/2) Rhat x B ; |Rhat x B| = B")
print("=> magnitude m/(q B R)(...) : dimension check [kg/(C*T*m)]*[m^2/s^2]")
kg,C,T,M,S = sp.symbols('kg C T M S')
print("T = kg/(C s) ; so kg/(C*T*m)*m^2/s^2 = kg*C*s/(C*kg*m)*m^2/s^2 = m/s : OK")

# thermal average: <vpar^2> = v^2/3, <vperp^2> = 2v^2/3 for isotropic
h("thermal average")
print("<vpar^2 + vperp^2/2> = v^2/3 + (2v^2/3)/2 = v^2/3 + v^2/3 = 2v^2/3 OK")
print("v_drift = m/(qBR) * 2/3 * v^2 ; with mv^2 = 3kT -> = 2kT/(qBR)  OK")

# ch6 L2311 suppl L146-176 R cutoff check
h("ch6 R-wave cutoff: n_R^2 = 1 - wpe^2/(w(w-wce)) = 0 -> w(w-wce) = wpe^2")
print("quadratic w^2 - wce*w - wpe^2 = 0 -> w = (wce + sqrt(wce^2+4wpe^2))/2 : correct")

# ch6 L2317-2324 the L-wave counterexamples
h("ch6 L2317-2324: L-wave n^2 counterexamples in the book")
for rp,ro in [(0.1,0.5),(2.0,0.5)]:
    nL2 = 1 - rp**2/(ro*(ro+1))   # wpe^2/(w(w+wce)) with wce=1
    nR2 = 1 - rp**2/(ro*(ro-1))
    print(f"wpe/wce={rp}, w/wce={ro}: nL^2={nL2:.4f}  nR^2={nR2:.4f}")
print("book says nL^2 ~ 0.987 and -4.33 : check")

# ch6 suppl L173 whistler R cutoff numeric
h("R cutoff numeric")
import numpy as np
wce=8.794e6; wpe=5.6414e7
w=(wce+np.sqrt(wce**2+4*wpe**2))/2
print("w =", w, "f =", w/2/np.pi)

# ch2 suppl L200 7430 constant
h("Debye constant 7430")
import numpy as np
e=1.602176634e-19; eps0=8.8541878128e-12
print("sqrt(eps0/e) =", np.sqrt(eps0/e), " head 7430*sqrt(1.602e-19/eps0)")

# ch5 suppl L26 claimed 25/(2.513e-6)=9.95e6
h("25/(2.513e-6) check")
print(25/2.5129e-6)

# ch2 suppl L38 applies 7430 formula with floor of decimal
h("ch2 suppl L204: 7430*sqrt(2/1e16)")
print(7430*np.sqrt(2/1e16), " vs sqrt2*1e-8*7430 =", np.sqrt(2)*1e-8*7430)
