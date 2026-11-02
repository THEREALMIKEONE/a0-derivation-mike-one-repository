# Finishing the picture inside the framework's own cosmology (no outside input).
# Family: layer scale mu ~ H^p, rho_vac = rho_v0 (H/H0)^(2p), vacuum <-> matter exchange, radiation separate.
# Data: Pantheon+ (full cov), DESI DR2 BAO (full cov, r_d computed), Planck 2018 distance priors (R, l_A, wb).
# RULES fixed before running:
#  W1 (why 1/3): the clock postulate omega_J = mu = H today  <=>  Omega_vac = 2/3  <=>  the medium removes exactly 1/3
#      of the bending (Lambda_5 = -2 A'^2). It is consistent if, fitting everything else, dchi2(Omega_vac = 2/3) < 4.
#  W2 (why now): report the redshift at which omega_J = H in the best fit and how long the |omega_J/H - 1| < 1.3% window lasts.
#  W3 (speed): joint p from SN+BAO+CMB and from the 145 discs (family distances) must agree within 2 sigma.
import os, sys, io, contextlib, math
suite, root = sys.argv[1], sys.argv[2]; os.environ["GM_ROOT"] = root; os.chdir(root)
G = {"__name__": "x"}; sys.argv[0] = os.path.join(suite, "run.py")
with contextlib.redirect_stdout(io.StringIO()):
    for n in ["p01_setup", "p02_geometry", "p03_map", "p04_galactic_gradient"]:
        f = os.path.join(suite, "cells", n + ".py"); exec(compile(open(f).read(), f, "exec"), G)
    src = open(os.path.join(suite, "cells", "p05_vacuum_along_S.py")).read()
    src = src.split("# §6–§7 — Dark energy as the cosmic slice")[0].rsplit("# =====", 1)[0]
    exec(compile(src, "p05defs", "exec"), G)
    f = os.path.join(suite, "cells", "p08_redshift.py"); exec(compile(open(f).read(), f, "exec"), G)
np = G["np"]; from scipy.integrate import solve_ivp; from scipy.optimize import minimize
E_orig = G["E"]; OG, NEFF = G["OMEGA_G_H2"], G["NEFF"]
LX = np.linspace(0, math.log(1 + 1.2e7), 6000); _cache = {}
def table(Om, p, h):
    key = (round(Om, 7), round(p, 7), round(h, 6))
    if key not in _cache:
        Or = OG * (1 + 0.2271 * NEFF) / h**2; Ov = 1 - Om - Or
        rhs = lambda lx, y: [3 + Or * math.exp(4 * lx - y[0]) - 3 * Ov * math.exp((p - 1) * y[0])]   # d ln u / d ln(1+z)
        sol = solve_ivp(rhs, (0, LX[-1]), [0.0], t_eval=LX, rtol=1e-8, atol=1e-10, method="LSODA")
        if len(_cache) > 400: _cache.clear()
        _cache[key] = np.exp(0.5 * sol.y[0]) if sol.success and len(sol.y[0]) == len(LX) else None
    return _cache[key]
def E_new(z, p):
    if p.get("model") != "fam": return E_orig(z, p)
    T = table(p["Om"], p["pp"], p["h"])
    return np.interp(np.log1p(np.asarray(z, float)), LX, T)
G["E"] = E_new
chi_grid, chi2_sn, chi2_bao, chi2_cmb, SN = G["chi_grid"], G["chi2_sn"], G["chi2_bao"], G["chi2_cmb"], G["SN"]
def chi2(Om, pp, h, wb):
    if not (0.15 < Om < 0.5 and -0.5 < pp < 1.2 and 0.55 < h < 0.85 and 0.018 < wb < 0.027): return 1e10
    P = {"model": "fam", "Om": Om, "pp": pp, "h": h, "wb": wb}
    if table(Om, pp, h) is None: return 1e10
    chi = chi_grid(P)
    return chi2_sn(P, chi, SN) + chi2_bao(P, chi) + chi2_cmb(P, chi)
def fit(fixed):
    names = [n for n in ("Om", "pp", "h", "wb") if n not in fixed]
    x0s = [{"Om": 0.31, "pp": 0.1, "h": 0.68, "wb": 0.02237}, {"Om": 0.33, "pp": 0.3, "h": 0.66, "wb": 0.0224}, {"Om": 0.30, "pp": 0.0, "h": 0.69, "wb": 0.0223}]
    best = None
    for x0 in x0s:
        f = lambda x: chi2(**{**fixed, **dict(zip(names, x))})
        r = minimize(f, [x0[n] for n in names], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-3, "maxiter": 5000, "maxfev": 5000})
        if best is None or r.fun < best.fun: best = r
    return best.fun, {**fixed, **dict(zip(names, best.x))}

import numpy as np
c_free, P_free = fit({})
Or = lambda h: OG * (1 + 0.2271 * NEFF) / h**2
print(f"[free] chi2 {c_free:.2f} p {P_free['pp']:+.4f} Om {P_free['Om']:.4f}")
rows=[]
for Om in np.linspace(0.290,0.316,9):
    c,P=fit({"Om":Om}); rows.append((1-Om-Or(P['h']),c)); print(f"  Omega_vac={rows[-1][0]:.4f} dchi2={c-c_free:+.2f}")
rows=np.array(sorted(rows)); f=np.linspace(rows[0,0],rows[-1,0],2000); cc=np.interp(f,rows[:,0],rows[:,1]); ins=f[cc<cc.min()+1]
print(f"[Omega_vac profile] {0.5*(ins.min()+ins.max()):.4f} +- {0.5*(ins.max()-ins.min()):.4f};  medium removes {1-0.5*(ins.min()+ins.max()):.3f} of the bending")
