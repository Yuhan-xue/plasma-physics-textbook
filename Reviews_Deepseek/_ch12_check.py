from math import pi, sqrt, log

e = 1.602176634e-19
eps0 = 8.8541878188e-12
me = 9.1093837139e-31

print("=== ch12 argon probe example ===")
Te_eV = 1 / 3.5
print("Te =", round(Te_eV, 4), "eV   (book says 0.29)")

r = 0.25e-3
L = 5e-3
A = 2 * pi * r * L
print("A =", "%.3e" % A, "m^2   (book says 7.85e-6)")

mi = 40 * 1.66053906660e-27
uB = sqrt(Te_eV * e / mi)
print("uB =", "%.4e" % uB, "m/s   (book says 8.4e2)")

Iisat = 0.15e-3
ne = Iisat / (0.6 * e * A * uB)
print("ne =", "%.3e" % ne, "m^-3   (book says 2.4e17)")

ld = sqrt(eps0 * Te_eV / (ne * e))
print("lambda_D =", "%.3e" % ld, "m   (book says 8.2e-6)")
print("rp/lambda_D =", round(r / ld, 1), "   (book says ~30)")

print()
print("=== floating potential, argon ===")
ratio = 0.5 * log(2 * pi * me / mi)
print("Vf - Vp = (Te/2)*ln(2*pi*me/mi) =", round(ratio, 3), "* Te   (book says -4.7 Te)")

print()
print("=== child-langmuir sheath estimate in the example ===")
Vs = 100.0
s = sqrt(2) / 3 * ld * (2 * Vs / Te_eV) ** 0.75
print("s =", "%.3e" % s, "m =", round(s * 1e6, 1), "um   (book says ~5e2 um)")
print("s/rp =", round(s / r, 2), "   (book says >= 1)")
