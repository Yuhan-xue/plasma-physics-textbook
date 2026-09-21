# -*- coding: utf-8 -*-
"""Symbolic verification of suspected errors."""
import sympy as sp
import numpy as np

def h(s): print("\n== " + s + " ==")

# ---------- ch2 L1367 area: Debye sphere particle number in book's own ch2 table ----------
h("ch2 table L1397-1400: N_D from the table's own n, lambda_D")
for name, n, T_eV, lam in [("ISM",1e6,1e-2,10.0),("ionosphere",1e12,1e-1,1e-3),
                           ("lab",1e18,10.0,1e-5),("ICF",1e32,1e4,1e-10)]:
    ND = 4/3*np.pi*n*lam**3
    print(f"{name:11s} N_D = {ND:.3g}")

# ---------- ch3: E x B with B along z, E along y -> E x B direction ----------
h("ch3 suppl L48-67 ExB direction: E=E yhat, B=B zhat -> E x B = ?")
# yhat x zhat = xhat  (right-handed). Correct.
print("yhat x zhat = xhat  -> +x, book says +x : OK")

# ---------- ch3 grad-B drift sign check ----------
h("ch3 L1676 grad-B drift sign from force F = -mu grad B")
# F = -mu grad B ; v_F = F x B /(q B^2) = -mu (grad B x B)/(qB^2) = mu (B x grad B)/(q B^2). OK matches book.

# ---------- ch4 suppl L250-313: u_i1 sign check via symbolic solve ----------
h("ch4 suppl ion-acoustic: solve the actual linearized system")
w,k,mi,ephi,n0,ee,Ti_,g_ = sp.symbols('w k m_i ephi n_0 e_ T_i gamma_i', positive=True)
vt2 = sp.Symbol('v_ti2', positive=True)
u = sp.symbols('u')
sol_u = sp.solve(sp.Eq(mi*w*u, k*ephi + g_*k**2*vt2/w*u), u)[0]
print("u =", sp.simplify(sol_u), "   book:", "k*e*phi1*w/(m_i*(w^2-gamma_i k^2 v_ti^2))")
print("identical?", sp.simplify(sol_u - ephi*k*w/(mi*(w**2-g_*k**2*vt2))) == 0)

# ---------- ch4 suppl: ambipolar Da ----------
h("ch4 suppl ambipolar: is Da = (mu_i D_e + mu_e D_i)/(mu_e+mu_i)?")
De,Di,mu_e,mu_i,dn,n = sp.symbols('D_e D_i mu_e mu_i dn n', positive=True)
E = -(De-Di)/n*dn/(mu_e+mu_i)
Gi = -Di*dn + mu_i*n*E
Da = sp.simplify(-Gi/dn)
print("Da = -Gamma_i/(dn) =", Da)
print("book  =", (mu_i*De+mu_e*Di)/(mu_e+mu_i))
print("equal?", sp.simplify(Da - (mu_i*De+mu_e*Di)/(mu_e+mu_i)) == 0)

# ---------- ch5 suppl Bennett ----------
h("ch5 suppl Bennett: check mu0 I^2 = 8 pi N k(T_e+T_i)")
mu0,I,N,kB,T = sp.symbols('mu_0 I N k_B T', positive=True)
print("standard Bennett: mu0 I^2 = 8 pi N k_B T (T = T_e+T_i)  ok")

# ---------- ch6 suppl whistler exact vs approx ----------
h("ch6 suppl whistler: exact R-wave dispersion vs book approx")
e=1.602176634e-19; me=9.1093837015e-31; eps0=8.8541878128e-12; c=2.99792458e8
def wpe_(n): return np.sqrt(n*e**2/(eps0*me))
def wce_(B): return e*B/me
n=1e12; B=5e-5
wp=wpe_(n); wc=wce_(B)
f=1e4; w=2*np.pi*f
nR2 = 1 - wp**2/(w*(w-wc))
print("nR^2 =", nR2, "-> n^2 >> 1 ok; n =", np.sqrt(nR2))
k_ex = w*np.sqrt(nR2)/c
k_ap = wp/c*np.sqrt(w/wc)
print("k exact =", k_ex, " k approx =", k_ap, " rel diff =", (k_ex-k_ap)/k_ex)
vp_ex = w/k_ex; vg_ex = c**2*k_ex/w
print("vp exact =", vp_ex, " book 3.95e6")
print("vg exact =", vg_ex, " 2*vp =", 2*vp_ex, " book 7.9e6")
# proper whistler: vg = 2 v_p holds exactly for w = c^2k^2 wc/wp^2 (k^2 law)
print("for k^2 law: vp = w/k, vg = 2w/k = 2 vp : correct relation")

# ---------- ch6 suppl wavepacket quote inconsistency ----------
h("ch6 suppl wavepacket: book quotes 0.311 then uses 0.315")
print("2.84/9 =", 2.84/9, " 2.8/9 =", 2.8/9, " true =", 2.8393/9)
x = 2.84/9
print("using 0.311: vp=", c/np.sqrt(1-0.311**2), " vg=", c*np.sqrt(1-0.311**2))
print("using 0.3155: vp=", c/np.sqrt(1-0.315478**2), " vg=", c*np.sqrt(1-0.315478**2))

# ---------- ch2 suppl L73 Gamma consistency claim ----------
h("ch2 suppl L73: claims Gamma ~ N_D^{-2/3} ~ 2e-3 AND Gamma=N_D^{-1}~1e-4")
nn=1e16; Te=2*e; lam=np.sqrt(eps0*Te/(nn*e**2))
a=(3/(4*np.pi*nn))**(1/3)
G = e**2/(4*np.pi*eps0*a*Te)
ND_sph = nn*4/3*np.pi*lam**3
print("true Gamma (a = Wigner-Seitz) =", G)
print("(1/3) N_D_sph^{-2/3} =", (1/3)*ND_sph**(-2/3))
print("N_D_sph^{-1} =", ND_sph**-1)
print("book claim 2e-3 vs true", G, " ratio =", G/2e-3)

# ---------- ch2 n^-1/2 scaling (book correct) ----------
h("ch2 L1467 ND ~ n^{-1/2} check")
print("ND = (4pi/3) n lam^3 ~ n * n^{-3/2} = n^{-1/2} : correct")

# ---------- ch5 suppl L67-71 vA/cs inconsistency: gamma ----------
h("ch5 suppl cs definition inconsistency")
print("book (1) uses cs = sqrt(2 kT/m) = 9.8e5 ; standard cs = sqrt(gamma kT_e/m_i)")

# ---------- ch1 suppl L132: v_t = sqrt(2kT/m) called '1D scale param' ----------
h("ch1 suppl L132: v_t = sqrt(2kT/m), claims v_t = sqrt(2) v_rms(1D)")
import numpy as np
print("1D rms = sqrt(kT/m); sqrt(2)*sqrt(kT/m) = sqrt(2kT/m) : consistent")
print("but NOTE v_t = sqrt(2 kT/m) equals 3D most-probable speed value coincidentally - book flags it")

# ---------- ch4 suppl L284-291 Poisson closure boxed formula ----------
h("ch4 suppl L290: boxed dispersion omega^2 = k^2(gamma_i kT_i/m_i + (kT_e/m_i)/(1+k^2 lamDe^2))")
# check by symbolic derivation
k,w,n0,ee,eps0_,kB,T_e,T_i,g_ = sp.symbols('k omega n_0 e epsilon_0 k_B T_e T_i gamma_i', positive=True)
lamDe = sp.sqrt(eps0_*kB*T_e/(n0*ee**2))
ni1 = ee*n0*k**2*sp.Symbol('phi1')/(sp.Symbol('m_i')*(w**2-g_*k**2*kB*T_i/sp.Symbol('m_i')))
ne1 = n0*ee*sp.Symbol('phi1')/(kB*T_e)
poisson = sp.Eq(-k**2*sp.Symbol('phi1'), ee/eps0_*(ne1-ni1))
sol = sp.solve(poisson, w**2)
print("w^2 =", [sp.simplify(s) for s in sol])
