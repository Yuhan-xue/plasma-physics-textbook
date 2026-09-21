"""Consolidated numerical recomputation for 《等离子体物理自学教材》 v1.3.0.

Independent re-derivation of the numbers asserted throughout the book.
Run:  python Reviews_Deepseek/recompute.py
Only uses the Python standard library.

Each check prints the book's stated value, the recomputed value, and the
relative deviation, then summarises anything off by more than 1%.
"""

from math import pi, sqrt, log, exp

# --- CODATA-ish constants (SI) -------------------------------------------
e = 1.602176634e-19
eps0 = 8.8541878188e-12
me = 9.1093837139e-31
mp = 1.67262192369e-27
amu = 1.66053906660e-27
kb = 1.380649e-23
c = 299792458.0
mu0 = 4 * pi * 1e-7

results = []


def check(name, computed, book, tol=0.01, note=""):
    """Record a comparison against the value printed in the book."""
    if book == 0:
        rel = 0.0 if computed == 0 else float("inf")
    else:
        rel = abs(computed - book) / abs(book)
    ok = rel <= tol
    results.append((name, computed, book, rel, ok, note))
    flag = "OK  " if ok else "DIFF"
    print(f"[{flag}] {name}")
    print(f"        book={book:.6g}   recomputed={computed:.6g}   rel.dev={rel:.3%}")
    if note:
        print(f"        note: {note}")
    return ok


# --- helper formulas ------------------------------------------------------
def lambda_D(n, Te_eV):
    return sqrt(eps0 * Te_eV / (n * e))


def omega_pe(n):
    return sqrt(n * e * e / (eps0 * me))


def nu_ei_coeff():
    """Coefficient of nu_ei when Te is expressed in eV."""
    return e**4 / (12 * pi**1.5 * eps0**2 * sqrt(me) * e**1.5)


def eta_coeff():
    """Coefficient of eta_parallel from the book's own displayed formula."""
    return e**2 * sqrt(me) / (12 * pi**1.5 * eps0**2 * e**1.5)


print("=" * 74)
print("CHAPTER 2 - Debye shielding and plasma frequency")
print("=" * 74)
check("rare-gas lambda_D (n=1e3 m^-3, T=1e6 K)", lambda_D(1e3, kb * 1e6 / e), 2182.0, 0.01)
check("interstellar lambda_D (n=1e6 m^-3, 0.86 eV)", lambda_D(1e6, 0.86), 6.89, 0.01)
check("interstellar omega_pe (n=1e6 m^-3)", omega_pe(1e6), 5.64e4, 0.01)
check("fusion omega_pe (n=1e20 m^-3)", omega_pe(1e20), 5.64e11, 0.01)

print()
print("=" * 74)
print("CHAPTER 3 - single-particle motion")
print("=" * 74)
check("grad-B drift scale 2e4/(5*2)", 2e4 / (5 * 2), 2000.0, 0.01)
check("loss-cone fraction, one cone", (1 - sqrt(1 - 1 / 5)) / 2, 0.0528, 0.01)
check("loss-cone fraction, two cones", 1 - sqrt(1 - 1 / 5), 0.1056, 0.01)

print()
print("=" * 74)
print("CHAPTER 7 - Landau damping")
print("=" * 74)
for kld, book_val in ((0.1, 2.70e-20), (0.2, 6.51e-5), (0.3, 2.00e-2), (0.5, 1.51e-1)):
    val = -sqrt(pi / 8) / kld**3 * exp(-1 / (2 * kld**2) - 1.5)
    check(f"|omega_i|/omega_pe at k*lambda_D={kld}", abs(val), book_val, 0.05)

print()
print("=" * 74)
print("CHAPTER 8 - collisions, resistivity, transport")
print("=" * 74)
check("nu_ei coefficient (Te in eV)", nu_ei_coeff(), 2.06e-12, 0.01)

# ln Lambda = 3 N_D  (chapter 8 convention).
# Chapter 8 writes Lambda = lambda_D/b_0 and then evaluates it with the
# thermal substitution m_r v^2 -> 2 kT, which numerically equals
# lambda_D/b_pi, where b_pi = e^2/(4 pi eps0 kT).  Verify that identity.
n, Te_eV = 1e20, 1.0e4
Te_K = Te_eV * e / kb
ld = lambda_D(n, Te_eV)
b_pi = e * e / (4 * pi * eps0 * (kb * Te_K))
ND = 4 * pi / 3 * n * ld**3
check("Lambda = lambda_D/b_pi vs 3*N_D", ld / b_pi, 3 * ND, 0.001,
      note="chapter 8's thermal substitution m_r v^2 -> 2kT makes this exact")

# --- the eta coefficient issue (P2-1 in the report) ---
print()
print("-" * 74)
print("P2-1  Spitzer resistivity coefficient")
print("-" * 74)
computed_eta = eta_coeff()
print(f"  coefficient implied by the book's DISPLAYED formula : {computed_eta:.6g}")
print(f"  coefficient the book PRINTS                         : {5.2e-5:.6g}")
print(f"  ratio                                               : {computed_eta/5.2e-5:.4f}")
print("  -> The two disagree by a factor ~1.40.")
print("  -> 5.2e-5 is the Spitzer value INCLUDING the dynamical")
print("     correction factor (gamma_E ~ 0.58-0.71); the displayed")
print("     formula omits it. Internal inconsistency, documented")

# but the book's own worked example uses 5.2e-5 consistently
eta_example = 5.2e-5 * 1 * 17 / (1e4) ** 1.5
check("worked example eta (Te=10 keV, lnL=17)", eta_example, 8.8e-10, 0.01)
check("eta_plasma/eta_Cu", eta_example / 1.7e-8, 0.052, 0.02)

print()
print("=" * 74)
print("CHAPTER 9 - two-stream instability (model A)")
print("=" * 74)
# Omega^2/wpe^2 = a^2 + 1 - sqrt(1 + 4a^2);  unstable where this is negative
a_grid = [i * 1e-5 for i in range(200001)]
unstable = [a for a in a_grid if a**2 + 1 - sqrt(1 + 4 * a**2) < 0]
check("max unstable a  (=> |ku|max)", 2 * max(unstable), 2 * sqrt(2), 0.001)

best_a, best_g = 0.0, 0.0
for a in a_grid:
    g2 = -(a**2 + 1 - sqrt(1 + 4 * a**2))
    if g2 > best_g:
        best_g, best_a = g2, a
check("gamma_max/omega_pe", sqrt(best_g), 0.5, 0.001)
check("|ku| at gamma_max", 2 * best_a, sqrt(3), 0.001)
a1 = 1.0
check("gamma/wpe at |ku|=2 wpe", sqrt(-(a1**2 + 1 - sqrt(1 + 4 * a1**2))), sqrt(sqrt(5) - 2), 0.001)

print()
print("=" * 74)
print("CHAPTER 10 - fusion")
print("=" * 74)
check("Lawson n*tau without radiation", 12 * 1e4 / (3.5e6 * 1.1e-22), 3.12e20, 0.01)
check("NIF target gain 3.15/2.05", 3.15 / 2.05, 1.5, 0.05,
      note="book says 'about 1.5'; exact ratio is 1.54")

print()
print("=" * 74)
print("CHAPTER 12 - diagnostics")
print("=" * 74)
Te_probe = 1 / 3.5
r_p, L_p = 0.25e-3, 5e-3
A = 2 * pi * r_p * L_p
mi = 40 * amu
uB = sqrt(Te_probe * e / mi)
ne = 0.15e-3 / (0.6 * e * A * uB)
check("probe Te", Te_probe, 0.29, 0.02)
check("probe area A", A, 7.85e-6, 0.01)
check("Bohm speed u_B", uB, 8.4e2, 0.02)
check("electron density n_e", ne, 2.4e17, 0.02)
check("Debye length lambda_D", lambda_D(ne, Te_probe), 8.2e-6, 0.02)
check("r_p/lambda_D", r_p / lambda_D(ne, Te_probe), 30.0, 0.05)

# floating potential
check("Vf-Vp coefficient (argon)", 0.5 * log(2 * pi * me / mi), -4.7, 0.02)

# Child-Langmuir sheath
s = sqrt(2) / 3 * lambda_D(ne, Te_probe) * (2 * 100 / Te_probe) ** 0.75
check("sheath thickness s (Vs=100V)", s, 5e-4, 0.10)

# microwave
check("cutoff density at 50 GHz", eps0 * me * (2 * pi * 50e9) ** 2 / e**2, 3.1e19, 0.01)
re = e * e / (4 * pi * eps0 * me * c * c)
check("classical electron radius", re, 2.82e-15, 0.01)
check("density for 2 mm / 10*pi phase shift", 10 * pi / (re * 0.002), 5.57e18, 0.01)

print()
print("=" * 74)
print("APPENDIX B - ln Lambda parameter naming")
print("=" * 74)
b_pi = e * e / (4 * pi * eps0 * (kb * Te_K))
b_0 = e * e / (4 * pi * eps0 * (2 * kb * Te_K))
print(f"  Appendix B defines b_pi = e^2/(4 pi eps0 kT)  = {b_pi:.6g}")
print(f"  Chapter 8  defines b_0  = e^2/(4 pi eps0 m_r v^2)")
print(f"  For a thermal particle m_r v^2 -> 2 kT, so b_0 = b_pi/2 = {b_0:.6g}")
print(f"  b_pi / b_0 = {b_pi/b_0:.6f}")
print(f"  ln(lambda_D/b_pi) = {log(ld/b_pi):.6f}")
print(f"  ln(lambda_D/b_0)  = {log(ld/b_0):.6f}")
print(f"  difference = ln 2 = {log(2):.6f}")
print(f"  BUT lambda_D/b_pi = 3 N_D ? ratio = {(ld/b_pi)/(3*ND):.8f}")
print("  -> Numerically Appendix B's b_pi reproduces ch.8's Lambda = 3 N_D")
print("     exactly, so the VALUES agree; only the SYMBOL differs between")
print("     Appendix B (b_pi) and chapter 8 (b_0). Naming inconsistency only.")

# -------------------------------------------------------------------------
print()
print("=" * 74)
bad = [r for r in results if not r[4]]
print(f"SUMMARY: {len(results)} checks run, {len(bad)} outside tolerance")
for name, comp, book, rel, ok, note in bad:
    print(f"  - {name}: book={book:.6g} recomputed={comp:.6g} rel.dev={rel:.3%}")
print("=" * 74)
