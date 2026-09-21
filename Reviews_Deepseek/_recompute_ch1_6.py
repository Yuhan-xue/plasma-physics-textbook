# -*- coding: utf-8 -*-
"""Independent recomputation of every numeric example in chapters 1-6."""
import numpy as np
from scipy import constants as C

e = C.e; me = C.m_e; mp = C.m_p; eps0 = C.epsilon_0; mu0 = C.mu_0; kB = C.k
h = lambda s: print("\n== " + s + " ==")

h("ch2 L1334 Debye length n=1e18 Te=2eV")
n, Te = 1e18, 2*e
lam = np.sqrt(eps0*Te/(n*e**2)); print("lambda_De =", lam, "m  book says 1.05e-5")
VD = 4/3*np.pi*lam**3; print("V_D =", VD, "book 4.85e-15"); print("N_D =", n*VD, "book 4.85e3")
print("N_D via 7430 sqrt(T/n):", 7430*np.sqrt(2/1e18))

h("ch2 L1397-1400 table: lambda_D and omega_pe")
for name, nn, TT in [("ISM",1e6,1e-2),("ionosphere",1e12,1e-1),("lab",1e18,10.0),("ICF",1e32,1e4)]:
    L = np.sqrt(eps0*TT*e/(nn*e**2)); w = np.sqrt(nn*e**2/(eps0*me))
    print(f"{name:11s} n={nn:.0e} T={TT:g}eV -> lam={L:.3e} m  wpe={w:.3e} rad/s")

h("ch2 L1467-1475 n=1e3, T=1e6 K")
n1, T1 = 1e3, 1e6
lam1 = np.sqrt(eps0*kB*T1/(n1*e**2)); print("lambda_D =", lam1, "book 2.18e3")
ND1 = 4/3*np.pi*n1*lam1**3; print("N_D =", ND1, "book 4.4e13")
print("check ND for T in eV:", 1e6*kB/e)

h("ch2 L1437-1443 Gamma vs N_D relation")
# a = (3/(4 pi n))^{1/3};  Gamma = e^2/(4 pi eps0 a kB T);  lam_D from same n,T
for G in [1e-2, 1e-1, 1.0, 10.0]:
    # invert: Gamma = (1/3) N_D^{-2/3}
    ND = (1/(3*G))**1.5
    print(f"Gamma={G:g} -> N_D = {ND:.4g}  (book: Gamma=1e-2 -> N_D~1.9e2)")
# direct check
nn, TT = 1e20, 1e4*e
lam = np.sqrt(eps0*TT/(nn*e**2)); ND = 4/3*np.pi*nn*lam**3
a = (3/(4*np.pi*nn))**(1/3); G = e**2/(4*np.pi*eps0*a*TT)
print(f"direct: n=1e20 Te=10keV -> N_D={ND:.4g}, Gamma={G:.4g}, (1/3)ND^-2/3={(1/3)*ND**(-2/3):.4g}")
print(f"1/(3 sqrt3) * Gamma^-3/2 = {1/(3*np.sqrt(3))*G**-1.5:.4g} vs N_D={ND:.4g}")

h("ch3 L1591-1611 gyro params B=5T, Te=Ti=10keV")
B, T = 5.0, 1e4*e
wce = e*B/me; print("wce =", wce, "book 8.8e11 ; fce =", wce/2/np.pi)
vte = np.sqrt(T/me); print("vte =", vte, "book 4.2e7")
print("rho_Le =", vte/wce, "book 4.8e-5")
mi = 2*mp; wci = e*B/mi; print("wci =", wci, "book 2.4e8; fci =", wci/2/np.pi, "book 38MHz")
vti = np.sqrt(T/mi); print("vti =", vti, "book 6.9e5")
print("rho_Li =", vti/wci, "book 2.9e-3 ; ratio =", (vti/wci)/(vte/wce))

h("ch3 L1727-1730 mirror loss cone Rm=3 -> theta")
Rm = 3.0; print("theta_m =", np.degrees(np.arcsin(np.sqrt(1/Rm))), "deg  (book suppl 35.3)")

h("ch5 L2132-2147 coronal Alfven speed n=1e14, B=1e-3")
n, B = 1e14, 1e-3
rho = n*mp; vA = B/np.sqrt(mu0*rho)
print("rho_m =", rho, "book 1.67e-13")
print("vA =", vA, "book 2.2e6")
print("sqrt(mu0 rho) =", np.sqrt(mu0*rho), "book 4.58e-10")

h("ch5 suppl L9-48 beta")
B, n, T = 5.0, 1e20, 1e4*e
p = 2*n*T; print("p =", p, "book 3.204e5")
B2_2mu0 = B**2/(2*mu0); print("B^2/2mu0 =", B2_2mu0, "book 9.95e6")
print("beta =", 2*mu0*p/B**2, "book 3.2e-2")
B, n, T = 0.1, 1e18, 4*e
p = 2*n*T; print("low-T p =", p, "book 1.28"); print("low-T beta =", 2*mu0*p/B**2, "book 3.2e-4")

h("ch5 suppl L63-71 Alfven tokamak n=1e20 m_D=2mp B=5T")
mD = 2*mp; n, B = 1e20, 5.0
rho = n*mD; vA = B/np.sqrt(mu0*rho)
print("rho_m =", rho, "book 3.344e-7")
print("vA =", vA, "book 7.7e6")
cs = np.sqrt(2*1e4*e/mD); print("cs(sqrt(2kT/m)) =", cs, "book 9.8e5; ratio", vA/cs)
cs_iso = np.sqrt(1e4*e/mD); print("cs(iso, sqrt(kTe/mi)) =", cs_iso)

h("ch5 suppl L73-92 Alfven lab & corona")
n, B = 1e18, 0.1; rho = n*mp; print("lab rho=",rho,"book 1.673e-9; vA=",B/np.sqrt(mu0*rho),"book 2.2e6")
print("lab cs =", np.sqrt(4*e/mp), "book 2.0e4")
n, B = 1e14, 1e-3; rho = n*mp; print("corona vA=",B/np.sqrt(mu0*rho))
T = 1e6*kB; print("corona T(K)=1e6 ->", T/e, "eV; cs(gamma=5/3) =", np.sqrt(5/3*T/mp), "book 1.5e5")

h("ch5 suppl L107-123 magnetic Reynolds")
eta = 5e-5*1*10/(1e4)**1.5; print("eta(10keV) =", eta, "book 5e-10")
Rm = 4e-7*np.pi*1e4*1/eta; print("Rm tokamak =", Rm, "book 2.5e7 (they used 4pi*1e-7 exact)")
Rm2 = mu0*1e4*1/eta; print("Rm exact =", Rm2)
print("tau_eta = mu0*1^2/eta =", mu0/eta, "book 2500 s ; Rm*L/v =", Rm2*1e-4)
eta2 = 5e-5*10/(4)**1.5; print("eta(4eV) =", eta2, "book 6.25e-5")
Rm3 = mu0*1e3*0.1/eta2; print("Rm lab =", Rm3, "book 2.0")
print("tau_eta =", mu0*0.01/eta2, "s  book 0.2ms")

h("ch5 suppl L306-330 Bennett")
I, N = 1e5, 1e19
kT = mu0*I**2/(8*np.pi*N); print("kT =", kT, "J  book 5.0e-17")
print("T =", kT/kB, "K  book 3.6e6 ;", kT/e, "eV book 310")
print("I for T=100eV:", np.sqrt(8*np.pi*N*100*e/mu0), "book 57kA")

h("ch6 suppl L9-82 Langmuir wave n=1e18 Te=4eV")
n, Te = 1e18, 4*e
wpe = np.sqrt(n*e**2/(eps0*me)); print("wpe =", wpe, "book 5.64e10; fpe =", wpe/2/np.pi, "book 8.98GHz")
print("lambda_pe = 2pi c/wpe =", 2*np.pi*C.c/wpe, "book 3.34e-2")
vte = np.sqrt(Te/me); print("vte =", vte, "book 8.39e5")
k = 2*np.pi/1e-2; print("k =", k, "book 628.3")
eps = 3*k**2*vte**2/wpe**2; print("eps =", eps, "book 2.62e-4")
print("sqrt(1+eps)-1 =", np.sqrt(1+eps)-1, "book 1.31e-4")
print("w =", wpe*np.sqrt(1+eps), "book 5.6407e10")
lamDe = vte/wpe; print("lamDe =", lamDe, "book 1.49e-5 ; k*lamDe =", k*lamDe, "book 9.35e-3")

h("ch6 suppl L85-141 cutoff n=1e12")
n = 1e12
wpe = np.sqrt(n*e**2/(eps0*me)); fpe = wpe/2/np.pi
print("wpe =", wpe, "book 5.64e7 ; fpe =", fpe, "book 8.97e6")
for f in [5e6, 20e6, 30e6]:
    if f > fpe:
        k = 2*np.pi*np.sqrt(f**2-fpe**2)/C.c; print(f"f={f/1e6:g}MHz k={k:.4f} rad/m  lambda={2*np.pi/k:.3f} m")
    else:
        print(f"f={f/1e6:g}MHz -> cutoff")

h("ch6 suppl L144-215 whistler B=5e-5, n=1e12")
B, n = 5e-5, 1e12
wce = e*B/me; wpe = np.sqrt(n*e**2/(eps0*me))
print("wce =", wce, "book 8.79e6 ; fce =", wce/2/np.pi, "book 1.40e6")
print("fpe =", wpe/2/np.pi)
wc_R = (wce+np.sqrt(wce**2+4*wpe**2))/2; print("R cutoff w =", wc_R, "book 6.09e7 ; f =", wc_R/2/np.pi, "book 9.69MHz")
w = 2*np.pi*1e4
k = wpe/C.c*np.sqrt(w/wce); print("k =", k, "book 1.59e-2 (approx formula)")
vp = w/k; print("vp =", vp, "book 3.95e6 ; vg=2vp =", 2*vp, "book 7.9e6")
# exact R-wave
nR2 = 1 - wpe**2/(w*(w-wce)); print("exact nR^2 =", nR2)
if nR2>0:
    kex = w*np.sqrt(nR2)/C.c; print("exact k =", kex, "vp =", w/kex, "vg =", C.c**2*kex/w)

h("ch6 suppl L217-271 phase/group")
n = 1e12; wpe = np.sqrt(n*e**2/(eps0*me)); w = 2*np.pi*20e6
k = np.sqrt(w**2-wpe**2)/C.c; print("k =", k, "book 0.374")
vp = w/k; print("vp =", vp, "book 3.36e8")
vg = C.c**2/vp; print("vg =", vg, "book 2.68e8 ; vp*vg =", vp*vg, "c^2 =", C.c**2)
print("v_p/c-1 =", vp/C.c-1)

h("ch6 suppl L348-379 wavepacket n=1e17 f=9GHz")
n = 1e17; wpe = np.sqrt(n*e**2/(eps0*me)); fpe = wpe/2/np.pi
print("fpe =", fpe, "book 2.8e9 ; 8.98*sqrt(n) =", 8.98*np.sqrt(n))
w = 2*np.pi*9e9; x = wpe**2/w**2
print("wpe/w =", wpe/w, "book 0.311 (they wrote 2.84/9=0.315)")
vp = C.c/np.sqrt(1-x); vg = C.c*np.sqrt(1-x)
print("vp =", vp, "book 3.16e8 ; vg =", vg, "book 2.85e8")
print("t =", 10/vg, "book 3.5e-8 ; t0 =", 10/C.c, "book 33.3e-9 ; diff =", 10/vg-10/C.c)

h("ch2 suppl L9-48 low-T argon n=1e18 Te=4eV")
n, Te = 1e18, 4*e
lam = np.sqrt(eps0*Te/(n*e**2)); print("lam =", lam, "book 1.49e-5")
print("num 8.854e-12*6.408e-19 =", 8.854e-12*6.408e-19, "book 5.673e-30")
wpe = np.sqrt(n*e**2/(eps0*me)); print("wpe =", wpe, "book 5.64e10 ; fpe =", wpe/2/np.pi, "book 8.97GHz")
ND = n*4/3*np.pi*lam**3; print("ND =", ND, "book 1.4e4")

h("ch2 suppl L50-76 n=1e16 Te=2eV L=0.1")
n, Te = 1e16, 2*e
lam = np.sqrt(eps0*Te/(n*e**2)); print("lam =", lam, "book 1.05e-4")
print("L/lam =", 0.1/lam, "book 952")
print("N_D(n*lam^3) =", n*lam**3, "book 1.16e4")
print("N_D(4pi/3) =", n*4/3*np.pi*lam**3)
# Gamma consistency claim
a = (3/(4*np.pi*n))**(1/3); G = e**2/(4*np.pi*eps0*a*Te)
print("Gamma =", G, "book claims ~2e-3 from ND^-2/3")
print("(1/3)ND^-2/3 with ND=1.16e4:", (1/3)*(n*4/3*np.pi*lam**3)**(-2/3))
print("N_D^-1 =", (n*4/3*np.pi*lam**3)**-1)

h("ch2 suppl L78-119 env comparison")
for name, n, T in [("ISM",1e6,0.86*e),("lab",1e18,4*e),("fusion",1e20,1e4*e)]:
    lam = np.sqrt(eps0*T/(n*e**2)); w = np.sqrt(n*e**2/(eps0*me))
    ND = n*4/3*np.pi*lam**3
    print(f"{name:7s} lam={lam:.4g} (book ISM 6.89, lab 1.49e-5, fus 7.4e-5) ; wpe={w:.4g} ; ND={ND:.4g}")

h("ch2 suppl L122-147 Debye potential table")
lam = 1.49e-5; Q4 = e/(4*np.pi*eps0)
for r in [10e-6, 50e-6, 100e-6, 1e-3]:
    phi0 = Q4/r; print(f"r={r*1e6:7.0f}um phi0={phi0:.3e} (book 1.44e-4/2.88e-5/1.44e-5/1.44e-6) r/lam={r/lam:.3f} exp={np.exp(-r/lam):.3e} phi={phi0*np.exp(-r/lam):.3e}")
print("Q/(4pi eps0) =", Q4)

h("ch2 suppl L192-221 debye calc")
for label,n,T in [("A glow",1e16,2.0),("B tok",1e20,1e4)]:
    lam = 7430*np.sqrt(T/n); print(f"{label}: lam={lam:.4g} ; ND={4*np.pi/3*n*lam**3:.4g}")
lam = 7430*np.sqrt(1e4/1e20); print("B: r/lam at 1cm =", 0.01/lam, "; exp =", np.exp(-0.01/lam))

h("ch3 suppl L9-46 gyro params")
B, T = 5.0, 1e4*e
for name, m in [("e",me),("p",mp),("D",2*mp)]:
    vth = np.sqrt(2*T/m); wc = e*B/m
    print(f"{name}: vth={vth:.4g} (book e5.93e7 p1.38e6 D9.79e5) wc={wc:.4g} f={wc/2/np.pi:.4g} r={vth/wc:.4g}")
print("book r_e=67um r_p=2.88mm r_D=4.08mm; f_p=76.2MHz f_D=38.2MHz")

h("ch3 suppl L48-67 ExB")
print("vE =", 1e3/5.0, "m/s  book 200")

h("ch3 suppl L70-105 loss cone")
Rm = 3; print("theta_m =", np.degrees(np.arcsin(1/np.sqrt(Rm))), "book 35.3")
v, vpar = 2e6, 1.8e6; vperp = np.sqrt(v**2-vpar**2)
print("vperp =", vperp, "book 8.72e5 ; sin^2 =", (vperp/v)**2, "book 0.19")

h("ch3 suppl L107-143 tokamak drifts")
T = 1e4*e; print("2kT/(eBR) =", 2*T/(e*5*2), "book 2.0e3 for both e and D")
# check the derivation: v = m/(qBR)(vpar^2 + vperp^2/2); iso: <vpar^2>=v^2/3, <vperp^2>=2v^2/3
print("<vpar^2+vperp^2/2> = v^2/3+v^2/3 =", 2/3.0, "v^2")
print("=> m/(qBR)*(2/3)v^2 = 2kT/(qBR) if mv^2=3kT:", 2/3*3, "/3")
# numerical direction check
print("v for D at 10keV:", np.sqrt(3*T/(2*mp)), "vs book 1e6")

h("ch3 suppl L249-282 mirror Rm=5")
Rm=5; th=np.degrees(np.arcsin(1/np.sqrt(Rm))); print("theta_lc =", th, "book 26.6")
c=np.cos(np.arcsin(1/np.sqrt(Rm)))
print("f_one = (1-cos)/2 =", (1-c)/2, "book 5.3% ; f_both =", 1-c, "book 10.6%")

h("ch4 suppl L57-123 Ohm's law estimates")
B,n,Te,L,u,j,om = 5.0,1e20,1e4*e,1.0,1e4,1e5,1e3
eta = 5e-5*10/(1e4)**1.5; print("eta =", eta)
print("eta j =", eta*j, "book 5e-5 ; eps_eta =", eta*j/(u*B), "book 1e-9")
jB_en = j*B/(e*n); print("jB/en =", jB_en, "book 3.1e4 ; eps_H =", jB_en/(u*B), "book 0.62")
pe = n*Te; gradpe = pe/L
print("pe =", pe, "book 1.602e5 ; gradpe/en =", gradpe/(e*n), "book 1e4 ; eps_p =", gradpe/(e*n)/(u*B), "book 0.2")
coef = me/(n*e**2); print("m_e/(n e^2) =", coef, "book 3.55e-13")
print("inertia =", coef*om*j, "book 3.6e-5 ; eps_I =", coef*om*j/(u*B), "book 7e-10")

h("ch4 suppl L177-211 ion sound")
T = 4*e; mp_ = mp
cs = np.sqrt(T/mp_); print("cs =", cs, "book 1.96e4")
vti = np.sqrt(2*0.1*e/mp_); print("vti =", vti, "book 4.37e3 ; ratio =", cs/vti, "book 4.5")

h("ch4 suppl L250-313 ion acoustic dispersion")
# check (A) vs (B) and the substitution u_i1 = k e phi1/(m_i(w^2-gamma k^2 vti^2)) * w
# verify algebra symbolically-ish
import sympy as sp
w,k,ephi,mi,phi1,n0,ee,Te_,Ti_,e0,vt2 = sp.symbols('w k ephi m phi1 n0 e T_e T_i eps0 vt2', positive=True)
# from book: m w u = k e phi + gamma k^2 vti^2 /w * u
# => u (m w - g k^2 vt2/w) = k e phi => u = k e phi w/(m(w^2 - g k^2 vt2))
expr = sp.simplify(k*ephi*2/w)  # placeholder
u = sp.symbols('u')
sol = sp.solve(sp.Eq(mi*w*u, k*ephi + 3*k**2*vt2/w*u), u)
print("u solved =", sp.simplify(sol[0]))
gam = sp.symbols('gamma_', positive=True)
sol2 = sp.solve(sp.Eq(mi*w*u, k*ephi + gam*k**2*vt2/w*u), u)
print("u(gamma) =", sp.simplify(sol2[0]))

h("ch4 suppl L125-175 ambipolar")
mu_e, mu_i, De, Di, Tev, Tiv = sp.symbols('mu_e mu_i D_e D_i T_e T_i', positive=True)
E = -(De-Di)/(mu_e+mu_i)*sp.Symbol('dn')/sp.Symbol('n')
Gam_i = -Di*sp.Symbol('dn') + mu_i*sp.Symbol('n')*E
print("Gamma_i =", sp.simplify(sp.expand(Gam_i)))
Da = sp.simplify(-(Gam_i).subs(sp.Symbol('dn'),-1)/1)
print("Da =", sp.simplify((mu_i*De+mu_e*Di)/(mu_e+mu_i) - Da))
print("book Da =", (mu_i*De+mu_e*Di)/(mu_e+mu_i))
print("with Einstein: mu_i*De = mu_i*mu_e kTe/e ; mu_e*Di = mu_e*mu_i kTi/e")
print("book says Da ~ mu_i k(Te+Ti)/e = Di(1+Te/Ti)")
print("numeric: Di=1e-3, Te=4,Ti=1 -> Da =", 1e-3*(1+4))
