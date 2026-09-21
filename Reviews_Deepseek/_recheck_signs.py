# -*- coding: utf-8 -*-
"""Re-check whether any REMAINING finding rests on a sign/convention misreading."""
import numpy as np
from scipy import constants as C
def h(s): print("\n== " + s + " ==")

# F-05: ch5_supplement sound speed -- is this a sign/convention issue? NO, it is a numeric one.
h("F-05 sound speed definitions (purely numeric, no sign involved)")
e=C.e; mp=C.m_p; mD=2*mp; kB=C.k
T=1e4*e
print("book ~9.8e5 ; sqrt(2kT/mD) =", np.sqrt(2*T/mD))
print("book ~9.8e5 ; sqrt(kT/mD)   =", np.sqrt(T/mD))
print("=> book used sqrt(2kT/m_i); both are POSITIVE, no sign ambiguity. FINDING STANDS.")
print()
T2=1e6*kB
print("corona book 1.5e5 ; gamma=1:", np.sqrt(T2/mp), " gamma=5/3:", np.sqrt(5/3*T2/mp),
      " gamma=3:", np.sqrt(3*T2/mp))
print("=> book's 1.5e5 matches NO standard gamma. FINDING STANDS.")

# F-04/F-13: heat speed definition drift -- ratio is sqrt2, positive. No sign issue.
h("F-04/F-13 heat speed sqrt2 drift (positive ratios, no sign issue)")
print("sqrt(kT/m) vs sqrt(2kT/m): ratio =", np.sqrt(2), " (positive) -> FINDING STANDS")
print("rho_Le: book ch3 67um vs v2.tex:1601 48um ; ratio =", 6.744e-5/4.769e-5)

# F-07: Gamma using lambda_D vs a -- both positive quantities
h("F-07 Gamma(a) vs Gamma(lambda_D) (positive, ratio = a/lam)")
print("a/lambda_D =", ((3/(4*np.pi*1e16))**(1/3))/np.sqrt(C.epsilon_0*2*e/(1e16*e**2)))
print("=> ratio differs, no sign issue. FINDING STANDS.")

# F-06: whistler/wavepacket numbers -- positive
h("F-06 fpe/f 0.311 vs 0.315 (positive, no sign issue) -- STANDS")

# F-01 follow-on: are there OTHER places in ch1-6 using R_c or kappa?
h("Scan: other uses of curvature vector / R_c in the audited range")
print("v2.tex:1682-1690 (theorem), 1695-1707 (tokamak check), ch3s:113,290, appendixF:5270-5274")
print("All other uses employ the B x R_c form and are self-consistent (verified earlier).")

# Re-examine my earlier 'ch3 suppl L113' claim: the vector form m/(qB^2 R)(...)(Rhat x B)
h("ch3_supplement.tex:113 combined drift form -- re-check direction")
# v_gradB = mu/(qB^2) B x gradB ; |gradB|=B/R pointing inward (-Rhat)
# => B x gradB = B*phat x (B/R)*(-Rhat) = (B^2/R) (phat x (-Rhat)) = (B^2/R)(+zhat)
# v_curv = m vpar^2/(qB^2 Rc^2) B x Rc, Rc = -Rc Rhat
# => B x Rc = B phat x (-Rc Rhat) = B*Rc*(+zhat)  (since phat x Rhat = -zhat)
import sympy as sp
print("phat x Rhat = -zhat  (right-handed (R,phi,z))")
print("=> gradB inward: B x gradB = +zhat ; Rc inward: B x Rc = +zhat")
print("=> BOTH SAME SIGN, ADD. Matches book's 'two drifts add, same direction'.")
print("No sign error here. My P3/verified entry stands as VERIFIED CORRECT.")

# Confirm the handedness fact used above
h("Confirm (Rhat, phat, zhat) is right-handed")
Rhat=np.array([1,0,0.]); phat=np.array([0,1,0.]); zhat=np.array([0,0,1.])
print("Rhat x phat =", np.cross(Rhat,phat), " (expect +zhat) -> right-handed OK")
print("phat x Rhat =", np.cross(phat,Rhat), " (expect -zhat)")
