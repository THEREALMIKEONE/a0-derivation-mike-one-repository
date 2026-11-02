# Open question 1: is our layer drifting? What the p-family means as an equation of state, against DESI DR2.
# rho_vac = rho_vac,0 (H/H0)^{2p}  =>  w_eff(a) = -1 - (1/3) dln rho/dln a = -1 - (2p/3) dlnH/dln a.
# dlnH/dln a < 0 always (H falls as the universe grows)  =>  w_eff > -1 for p > 0: never phantom.
import math
from scipy.optimize import brentq
Om0, Or0 = 0.315, 9.15e-5
def H2(a, p, OL0):  # solve the family self-consistently: E^2 = Om a^-3 + Or a^-4 + OL0 E^{2p}
    f = lambda E2: Om0 * a**-3 + Or0 * a**-4 + OL0 * E2**p - E2
    return brentq(f, 1e-6, 1e12)
def w_eff(z, p):
    OL0 = 1 - Om0 - Or0
    a = 1 / (1 + z); h = 1e-4
    lnE = lambda la: 0.5 * math.log(H2(math.exp(la), p, OL0))
    dlnH = (lnE(math.log(a) + h) - lnE(math.log(a) - h)) / (2 * h)
    return -1 - (2 * p / 3) * dlnH
print(" z    w_eff(p=0.011)  w_eff(p=0.016, +1 sigma)   DESI DR2+CMB+DESY5 CPL (w0=-0.752, wa=-0.86)")
for z in [0, 0.3, 0.5, 1.0, 2.0]:
    cpl = -0.752 + (-0.86) * z / (1 + z)
    print(f"{z:4.1f}   {w_eff(z, 0.011):+.4f}          {w_eff(z, 0.016):+.4f}                    {cpl:+.3f}")
# where the CPL best fit crosses -1
zc = brentq(lambda z: -0.752 - 0.86 * z / (1 + z) + 1, 0.01, 3)
print(f"\nCPL best fit crosses w = -1 at z = {zc:.2f}; the p-family (and v3's ghost-free S(t)) never goes below -1.")
