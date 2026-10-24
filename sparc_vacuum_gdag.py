# ============================================================
# SPARC test of g† = c H_Λ / 2π (horizon temperature of the vacuum, not the expansion rate)
# Rule FROZEN 2026-10-23 before running (same as the paper's GF3):
#   H_Λ = H0 √Ω_Λ with H0 = 67.4, Ω_Λ = 0.685 (Planck 2018)  ->  g† = 8.63e-11 m/s²
#   TEST: ΔBIC(fixed − free) < 10 on the paper's own SPARC likelihood (cell p04, unchanged)
#   DIAGNOSTIC: same ΔBIC for the paper's value cH0/2π (GF3 as published)
# ============================================================
import os, sys, math
suite = sys.argv[1]; os.environ["GM_ROOT"] = sys.argv[2]; os.chdir(sys.argv[2])
G = {}
for name in ["p01_setup", "p02_geometry", "p03_map", "p04_galactic_gradient"]:
    f = os.path.join(suite, "cells", name + ".py")
    exec(compile(open(f, encoding="utf-8").read(), f, "exec"), G)
np = G["np"]; minimize = G["minimize"]; C = G["C"]; marg_nll = G["marg_nll"]
law = lambda gb, gd: gb / (1 - np.exp(-np.sqrt(gb / gd)))
nll = lambda lt, lg: marg_nll(np.log10(law(C["gbar_grid"], 10**lg)), lt, C)
n = len(C["ly"])
free = minimize(lambda x: nll(x[1], x[0]) if -11 < x[0] < -9 else 1e30, [-10.0, math.log(0.1)], method="Nelder-Mead")
def dbic(gd):
    b = minimize(lambda x: nll(x[0], math.log10(gd)), [math.log(0.1)], method="Nelder-Mead")
    return (2 * b.fun + math.log(n)) - (2 * free.fun + 2 * math.log(n))
H0 = 67.4e3 / 3.0856775814913673e22; c = 299792458.0
gL = c * H0 * math.sqrt(0.685) / (2 * math.pi); g0 = c * H0 / (2 * math.pi)
print(f"\nFREE g† = {10**free.x[0]:.3e}   N points = {n}, galaxies = {C['n_gal']}")
d1 = dbic(gL); d0 = dbic(g0)
print(f"TEST  g† = cH_Λ/2π = {gL:.3e}:  ΔBIC = {d1:+.1f}  -> {'PASS' if d1 < 10 else 'FAIL'}")
print(f"DIAG  g† = cH0/2π  = {g0:.3e}:  ΔBIC = {d0:+.1f}")
