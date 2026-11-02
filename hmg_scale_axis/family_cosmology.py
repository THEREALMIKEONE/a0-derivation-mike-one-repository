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
c_free, P_free = fit({})
Or = lambda h: OG * (1 + 0.2271 * NEFF) / h**2
print(f"[free]   chi2 = {c_free:.2f}  Om = {P_free['Om']:.4f}  p = {P_free['pp']:+.3f}  h = {P_free['h']:.4f}  wb = {P_free['wb']:.5f}  -> Omega_vac today = {1-P_free['Om']-Or(P_free['h']):.4f}")
c_l, P_l = fit({"pp": 0.0})
print(f"[p = 0, LCDM] chi2 = {c_l:.2f} (dchi2 {c_l-c_free:+.2f})  Om = {P_l['Om']:.4f}  Omega_vac = {1-P_l['Om']-Or(P_l['h']):.4f}")
# W1: Omega_vac = 2/3 exactly, p free  (Om = 1/3 - Or(h): profile on h)
def chi2_23(pp, h, wb): return chi2(1/3 - Or(h), pp, h, wb)
best = min((minimize(lambda x: chi2_23(*x), x0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-3, "maxiter": 4000})
            for x0 in ([0.1, 0.67, 0.0224], [0.3, 0.65, 0.0224], [0.0, 0.68, 0.0223])), key=lambda r: r.fun)
print(f"[W1] Omega_vac = 2/3: chi2 = {best.fun:.2f}, dchi2 = {best.fun-c_free:+.2f} ({'PASS' if best.fun-c_free < 4 else 'FAIL'}); p = {best.x[0]:+.3f}, h = {best.x[1]:.4f}")
best_l = min((minimize(lambda x: chi2(1/3 - Or(x[0]), 0.0, x[0], x[1]), x0, method="Nelder-Mead") for x0 in ([0.67, 0.0224],)), key=lambda r: r.fun)
print(f"     same with p = 0 forced: dchi2 = {best_l.fun-c_free:+.2f}")
# profile in p (Om, h, wb free) for the expansion-side p
prof = []
for pp in np.linspace(-0.1, 0.5, 25):
    c, P = fit({"pp": pp}); prof.append((pp, c, P["Om"], P["h"]))
prof = np.array(prof); fine = np.linspace(-0.1, 0.5, 2401); cc = np.interp(fine, prof[:, 0], prof[:, 1]); m = cc.min()
ins = fine[cc < m + 1]; p_exp, s_exp = 0.5 * (ins.min() + ins.max()), 0.5 * (ins.max() - ins.min())
print(f"[expansion SN+BAO+CMB] p = {p_exp:+.3f} +- {s_exp:.3f}")
# Omega_vac along the p-profile (to show what 'the vacuum fraction' is inside the family)
print("     Omega_vac today along the profile:", ", ".join(f"p={a:+.2f}: {1-b-Or(c):.3f}" for a, _, b, c in prof[::4]))
# W2: why now
T = table(P_free["Om"], P_free["pp"], P_free["h"]); z = np.expm1(LX); Ov0 = 1 - P_free["Om"] - Or(P_free["h"])
Ov = Ov0 * T**(2 * P_free["pp"]) / T**2; ratio = np.sqrt(1.5 * Ov)       # omega_J / H
zeq = float(np.interp(1.0, ratio[::-1], z[::-1])) if ratio[0] >= 1 else None
print(f"[W2] omega_J/H today = {ratio[0]:.4f}", f"; omega_J = H at z = {zeq:.3f}" if zeq is not None else "; omega_J = H lies in the future")
# cosmic time window where |ratio-1| < 1.3%  (use the family's own a(t), integrate forward a bit with LCDM-like extrapolation)
tz = lambda zz: float(np.trapezoid(1 / ((1 + z[z <= zz]) * T[z <= zz]), z[z <= zz])) if zz > 0 else 0.0
H0s = 100 * P_free["h"] / 3.0857e19 * 3.156e16  # per Gyr
t0 = float(np.trapezoid(1 / ((1 + z) * T), z)) / H0s
print(f"     age t0 = {t0:.2f} Gyr; omega_J t0 = {ratio[0]*T[0]*t0*H0s:.3f}")
# W3: discs with family distances
ZG = G["ZG"]; Eg = E_new(ZG, {"model": "fam", "Om": P_free["Om"], "pp": P_free["pp"], "h": P_free["h"]})
LNZ = np.log1p(ZG)
def chi_of(E_):
    f_ = (1 + ZG) / E_; return np.concatenate([[0], np.cumsum(0.5 * (f_[1:] + f_[:-1]) * np.diff(LNZ))])
E_ref = np.sqrt(0.3 * (1 + ZG)**3 + 0.7); chi_ref = chi_of(E_ref)
locA, locB, _fit, _disc_gbar, _rows, _load = G["locA"], G["locB"], G["_fit"], G["_disc_gbar"], G["_rows"], G["_load"]
rawA = [(float(q[1]), 10**float(q[4]), float(q[5]), float(q[6])) for q in (l.split() for l in _rows("price21.txt"))]
rawB = {nm: _load(fn) for nm, fn in (("KROSS", "kross_rout.csv"), ("KMOS3D", "kmos3d_rout.csv"), ("KGES", "kges_rout.csv"))}
rng = np.random.default_rng(5)
def disc_R(Eg):
    chi = chi_of(Eg); fr = lambda zz: np.interp(zz, ZG, chi) / np.interp(zz, ZG, chi_ref); lnE = lambda zz: math.log(np.interp(zz, ZG, Eg))
    out = []; h = []
    for zz, Mb, Re, fdm in rawA:
        fd = 1 - (1 - fdm) * fr(zz)
        if 0.02 < fd < 0.98: h.append((lnE(zz), math.log(_disc_gbar(Mb, Re) / math.log(fd)**2)))
    h = np.array(h); RA = lambda hh, ll: float((hh[:, 0] * (hh[:, 1] - ll.mean())).sum() / (hh[:, 0]**2).sum())
    out.append(("Price+21", RA(h, locA), float(np.std([RA(h[rng.integers(0, len(h), len(h))], locA[rng.integers(0, len(locA), len(locA))]) for _ in range(2000)]))))
    for nm, d in rawB.items():
        d = d.copy(); d[:, 2] = d[:, 2] / np.array([fr(zz) for zz in d[:, 0]])
        lo, hi = np.percentile(d[:, 3], [10, 90]); L = locB[(locB[:, 2] >= lo) & (locB[:, 2] <= hi)]; zm = float(np.median(d[:, 0]))
        RB = lambda hh, LL: (_fit(hh[:, 1], hh[:, 2]) - _fit(LL[:, 0], LL[:, 1])) / lnE(zm)
        out.append((nm, RB(d, L), float(np.std([RB(d[rng.integers(0, len(d), len(d))], L[rng.integers(0, len(L), len(L))]) for _ in range(150)]))))
    w = np.array([1 / s**2 for *_, s in out]); Rv = np.array([r for _, r, _ in out])
    return out, float((w * Rv).sum() / w.sum()), float(1 / math.sqrt(w.sum()))
samp, pd, sd = disc_R(Eg)
t = abs(p_exp - pd) / math.hypot(s_exp, sd)
pj = (p_exp / s_exp**2 + pd / sd**2) / (1 / s_exp**2 + 1 / sd**2); sj = (1 / s_exp**2 + 1 / sd**2)**-0.5
print(f"[W3] discs p = {pd:+.3f} +- {sd:.3f}; expansion p = {p_exp:+.3f} +- {s_exp:.3f}; {t:.1f} sigma apart -> {'PASS' if t < 2 else 'FAIL'}")
print(f"     joint p = {pj:+.3f} +- {sj:.3f}  (p=0: {abs(pj)/sj:.1f} sigma, p=1: {abs(pj-1)/sj:.1f} sigma); per sample:", ", ".join(f"{n} {r:+.2f}+-{s:.2f}" for n, r, s in samp))
# predicted galactic scale at z = 1, 2 relative to today (g ~ E^p)
for zz in (1.0, 2.0): print(f"     g(z={zz:.0f})/g(0) = {np.interp(zz, ZG, Eg)**pj:.3f}  (expansion reading cH/2pi: {np.interp(zz, ZG, Eg):.2f})")
