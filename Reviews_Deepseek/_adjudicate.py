# -*- coding: utf-8 -*-
"""Settle the (b.grad)b sign convention and re-verify F-02."""
import numpy as np

def h(s): print("\n== " + s + " ==")

h("A. Circular field line, R0=2, evaluate (b.grad)b = db/ds numerically")
# Field line: circle of radius R0 in xy-plane, b = phi-hat
# parametrize by arc length s: phi = s/R0
R0 = 2.0
def b_of_s(s):
    phi = s/R0
    return np.array([-np.sin(phi), np.cos(phi), 0.0])   # phi-hat
ds = 1e-6
for s0 in [0.0, 0.7*R0, 2.3*R0]:
    dbds = (b_of_s(s0+ds) - b_of_s(s0-ds))/(2*ds)
    phi = s0/R0
    pos = R0*np.array([np.cos(phi), np.sin(phi), 0.0])
    to_center = -pos/np.linalg.norm(pos)     # unit vector pointing toward origin
    print(f"s={s0:.3f}  (b.grad)b = {np.round(dbds,6)}")
    print(f"          to_center = {np.round(to_center,6)}")
    print(f"          dot = {np.dot(dbds, to_center):+.6f}   |dbds| = {np.linalg.norm(dbds):.6f} (1/R0={1/R0})")

h("B. Same via the operator (b.grad)b in cylindrical coords")
# b = phi-hat.  (b.grad) = (1/R) d/dphi.  d(phi-hat)/dphi = -R-hat
# => (b.grad)b = (1/R)(-R-hat) = -R-hat/R
# At phi=0: -R-hat = -(1,0,0) = (-1,0,0) ; to_center at (R,0,0) is also (-1,0,0).
print("(b.grad)b = -Rhat/R  ->  at phi=0 equals (-1,0,0)")
print("to_center at (2,0,0) = (-1,0,0)")
print("=> THESE POINT THE SAME WAY.  (b.grad)b points TOWARD the center.")
print("NOTE: -Rhat is NOT 'toward -x' generally; -Rhat IS the inward radial direction.")
print("My earlier script printed '-Rhat/R' and I mis-read '-Rhat' as 'away from center'.")

h("C. Cross-check with a concrete matrix/vector calculation of kappa")
phi = 0.0
Rhat = np.array([np.cos(phi), np.sin(phi), 0.0])
phat = np.array([-np.sin(phi), np.cos(phi), 0.0])
print("Rhat =", Rhat, " phat =", phat)
print("(b.grad)b = -Rhat/R0 =", -Rhat/R0, " -> parallel to -Rhat = toward center:",
      np.allclose(-Rhat/R0, -Rhat/R0))
print("R_c (pointing toward center, magnitude R_c=R0) = -R0*Rhat =", -R0*Rhat)
print("kappa = (b.grad)b =", -Rhat/R0)
print("Is kappa parallel to R_c?", np.allclose(np.cross(-Rhat/R0, -R0*Rhat), 0))

h("D. Now test the book's two forms on its own tokamak example")
B = 5.0; Rc = 2.0
def formA(phi):
    Rhat = np.array([np.cos(phi), np.sin(phi), 0.0])
    phat = np.array([-np.sin(phi), np.cos(phi), 0.0])
    Bv = B*phat
    Rcv = -Rc*Rhat                      # toward center
    return np.cross(Bv, Rcv)/(B**2*Rc**2)
def formB(phi):
    Rhat = np.array([np.cos(phi), np.sin(phi), 0.0])
    phat = np.array([-np.sin(phi), np.cos(phi), 0.0])
    b = phat
    kappa = -Rhat/Rc                    # = (b.grad)b
    return np.cross(b, kappa)/B
for phi in [0.0, 0.7, 2.3]:
    a, b_ = formA(phi), formB(phi)
    print(f"phi={phi}: formA={np.round(a,6)}  formB={np.round(b_,6)}  equal={np.allclose(a,b_)}")

h("E. Re-verify F-02 independently (ch2_supplement.tex:71-73)")
e=1.602176634e-19; eps0=8.8541878128e-12; kB=1.380649e-23
n=1e16; Te=2*e
lam=np.sqrt(eps0*Te/(n*e**2))
a=(3/(4*np.pi*n))**(1/3)
G_a   = e**2/(4*np.pi*eps0*a*Te)          # book's own ch2 definition (Wigner-Seitz)
G_lam = e**2/(4*np.pi*eps0*lam*Te)        # the supplement's stated formula
ND_sph = n*4/3*np.pi*lam**3               # the book's N_D convention (v2.tex:1429)
ND_bare= n*lam**3                         # what the supplement actually computed
print("lambda_D      =", lam)
print("N_D (sphere)  =", ND_sph, "  <- book's own convention")
print("N_D (bare)    =", ND_bare, "  <- what supplement line 71 computed")
print("ratio         =", ND_sph/ND_bare, " (should be 4pi/3 =", 4*np.pi/3, ")")
print()
print("Gamma via a       =", G_a)
print("Gamma via lambda  =", G_lam)
print("(1/3) N_D_sph^-2/3 =", (1/3)*ND_sph**(-2/3))
print("book claims Gamma ~ 2e-3")
print("=> ratio book/true =", 2e-3/G_a)
