# KiDS-1000 isolated lenses on the vacuum clock (Section 6 of "Hyperconical gravity applied to the renormalization scale").
# Mike One, 30 October 2026. Pass rules written before running.
# Usage: python kids_lens_tests.py <path to a lensing_RAR_KiDS checkout>   (github.com/robertmonjo/lensing_RAR_KiDS)
#
# Part 1. One scale per population with the derived relation g = gN / (1 - exp(-sqrt(gN/g_x))):
#   isolated late-type and blue lenses PASS if ln(g_x/g_dag) is within 2 sigma of 0 (g_dag = 1.056e-10).
#   Inner vs outer halves of each profile (one scale fits all radii if within 2 sigma). MICE LCDM mock as control.
# Part 2. Mass corrections for early types from Mistele et al. 2024 (JCAP 04, 020): stellar mass x 1.4, hot gas
#   M_hot/M* = 10^-5.414 M*^0.47 (their Eq. 4.1). Reports how much of the early/red excess they remove.
# Part 3. Changed clock vs added neighbour mass, same environment profile (matter correlation function,
#   r0 = 5 Mpc, gamma = 1.8, one free amplitude b). A reading wins if chi2 is lower by more than 4 on both early
#   and red lenses with 0 < b < 5. Robustness: r0 = 4, 6 Mpc; gamma = 1.6, 2.0; clock reading with neighbours
#   only beyond r_core = 0.1, 0.3, 0.6 Mpc.
import sys, math
import numpy as np
sys.path.insert(0, sys.argv[1]); import hmg_model as h


GDAG = 1.056e-10

def curve(gN, gx):
    x = np.sqrt(gN / gx)
    return gN / (-np.expm1(-x))

def fit(gN, go, ge):
    lg = np.linspace(math.log(3e-12), math.log(3e-9), 3000)
    chi = np.array([np.sum(((go - curve(gN, math.exp(l))) / ge) ** 2) for l in lg])
    i = chi.argmin(); best = lg[i]; m = chi - chi[i]
    inside = lg[m <= 1.0]
    sig = (inside.max() - inside.min()) / 2      # 1-sigma on ln g_x
    return math.exp(best), sig, chi[i], len(gN)

def report(name, gN, go, ge):
    gx, sig, c2, n = fit(gN, go, ge)
    dS = -math.log(gx / GDAG)
    print(f"  {name:22s} g_x = {gx:.3e}  ({gx / GDAG:5.2f} x g_dag)  ln(g_x/g_dag) = {-dS:+.2f} +- {sig:.2f}   chi2 = {c2:5.1f} / {n - 1}")
    return dS, sig

def halves(name, gN, go, ge):
    o = np.argsort(gN); k = len(gN) // 2
    lo, hi = o[:k], o[k:]
    a = report(name + " outer", gN[lo], go[lo], ge[lo])
    b = report(name + " inner", gN[hi], go[hi], ge[hi])
    z = abs(a[0] - b[0]) / math.hypot(a[1], b[1])
    print(f"  {'':22s} inner vs outer: {z:.1f} sigma -> {'PASS' if z < 2 else 'FAIL'}")
    return a, b

print("PART 1\nA. One clock per population, KiDS-1000 isolated lenses")
pops = {}
for name, f in [("late (Sersic n<2.5)", "Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt"),
                ("early (Sersic n>=2.5)", "Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt"),
                ("blue (u-r)", "Fig-8_RAR-KiDS-isolated_Colorbin_1.txt"),
                ("red (u-r)", "Fig-8_RAR-KiDS-isolated_Colorbin_2.txt")]:
    pops[name] = h.load_morph(f)
    pops[name + "_fit"] = report(name, *pops[name])

print("\nB. One scale for all radii: inner vs outer half")
for name in ["early (Sersic n>=2.5)", "red (u-r)", "late (Sersic n<2.5)", "blue (u-r)"]:
    halves(name, *pops[name])

print("\nC. Isolated lenses on the vacuum clock: late and blue at ln(g_x/g_dag) = 0")
for name in ["late (Sersic n<2.5)", "blue (u-r)"]:
    dS, s = pops[name + "_fit"]
    print(f"  {name:22s} {abs(dS) / s:.1f} sigma from 0 -> {'PASS' if abs(dS) < 2 * s else 'FAIL'}")

gap_k = pops["early (Sersic n>=2.5)_fit"][0] - pops["late (Sersic n<2.5)_fit"][0]
gap_k_s = math.hypot(pops["early (Sersic n>=2.5)_fit"][1], pops["late (Sersic n<2.5)_fit"][1])
gap_c = pops["red (u-r)_fit"][0] - pops["blue (u-r)_fit"][0]
gap_c_s = math.hypot(pops["red (u-r)_fit"][1], pops["blue (u-r)_fit"][1])
print(f"\n  KiDS gap in ln(g_x), early - late: {-gap_k:+.2f} +- {gap_k_s:.2f};  red - blue: {-gap_c:+.2f} +- {gap_c_s:.2f}")

print("\nD. Control: MICE LCDM mock, same fits")
mice = {}
for name, f in [("late", "morph_sersic_late_allbins.txt"), ("early", "morph_sersic_early_allbins.txt"),
                ("blue", "morph_color_blue_allbins.txt"), ("red", "morph_color_red_allbins.txt")]:
    d = np.loadtxt(h.DATA_DIR / f)
    gN, go = 10 ** d[:, 0], 10 ** d[:, 1]
    ge = go * np.log(10) * (d[:, 3] - d[:, 2]) / 2
    mice[name] = report("MICE " + name, gN, go, ge)
for a, b, kgap, ks in [("early", "late", gap_k, gap_k_s), ("red", "blue", gap_c, gap_c_s)]:
    g = mice[a][0] - mice[b][0]; gs = math.hypot(mice[a][1], mice[b][1])
    z = abs(g - kgap) / math.hypot(gs, ks)
    print(f"  MICE gap in ln(g_x) {a} - {b}: {-g:+.2f} +- {gs:.2f}  vs KiDS {-kgap:+.2f}: {z:.1f} sigma apart -> "
          f"{'LCDM reproduces the gap' if z < 2 else 'LCDM does not reproduce the gap'}")

print("\nPART 2. Mass corrections (Mistele et al. 2024)")
from scipy.optimize import brentq
def mstar_from_mbar(mb): return brentq(lambda m: m*(1+h.f_cold(m)) - mb, 1e8, 1e13)
for name, f, mbar in [("early (Sersic)", "Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt", 2.666e10),
                      ("red (u-r)", "Fig-8_RAR-KiDS-isolated_Colorbin_2.txt", 1.768e10)]:
    ms = mstar_from_mbar(mbar); mhot = ms * 10**(-5.414) * ms**0.47
    k = (1.4*ms + mhot) / mbar
    gN, go, ge = h.load_morph(f)
    print(f"{name}: M* = {ms:.2e}, hot gas = {mhot/ms:.2f} M*, gN multiplied by k = {k:.2f}")
    a0_ = report(name + " raw", gN, go, ge)
    b0_ = report(name + " corrected", k*gN, go, ge)
    print(f"   corrections remove {100*(a0_[0]-b0_[0])/a0_[0]:.0f}% of the excess")

print("\nPART 3. Changed clock (MED) vs added neighbour mass (NB)")
from scipy.special import gamma as Gf
G = 6.674e-11; MSUN = 1.989e30; MPC = 3.0857e22
GDAG = 1.056e-10
H0 = 67.4e3 / MPC; RHO_CRIT = 3 * H0**2 / (8 * math.pi * G)
RHO_VAC = 0.685 * RHO_CRIT; RHO_M = 0.315 * RHO_CRIT
R0 = 5.0 * MPC; GAM = 1.8

def xibar(r): return 3 / (3 - GAM) * (r / R0) ** (-GAM)
def sigma_xi(R):  # projected xi for a power law
    return R0**GAM * R**(1 - GAM) * math.sqrt(math.pi) * Gf((GAM - 1) / 2) / Gf(GAM / 2)

def m_med(gN, b, R):
    ratio = RHO_M * (1 + b * xibar(R)) / RHO_VAC
    return curve(gN, GDAG * np.sqrt(np.maximum(1.0, ratio)))
def m_nb(gN, b, R):
    dsig = b * RHO_M * sigma_xi(R) * (GAM - 1) / (3 - GAM)
    return curve(gN, GDAG) + 4 * G * dsig

def chi2(model, gN, go, ge, b, R): return np.sum(((go - model(gN, b, R)) / ge) ** 2)
def best(model, gN, go, ge, R):
    bs = np.linspace(0, 8, 4001); c = np.array([chi2(model, gN, go, ge, b, R) for b in bs])
    i = c.argmin(); ok = bs[c - c[i] <= 1]
    return bs[i], (ok.max() - ok.min()) / 2, c[i]

MBAR = {"late": 1.479e10, "early": 2.666e10, "blue": 1.277e10, "red": 1.768e10}
FILES = {"late": "Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt", "early": "Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt",
         "blue": "Fig-8_RAR-KiDS-isolated_Colorbin_1.txt", "red": "Fig-8_RAR-KiDS-isolated_Colorbin_2.txt"}
res = {}
for p in ["early", "red", "late", "blue"]:
    gN, go, ge = h.load_morph(FILES[p])
    R = np.sqrt(G * MBAR[p] * MSUN / gN)
    c0 = np.sum(((go - curve(gN, GDAG)) / ge) ** 2)
    bm, sm, cm = best(m_med, gN, go, ge, R)
    bn, sn, cn = best(m_nb, gN, go, ge, R)
    res[p] = (cm, cn)
    print(f"{p:5s}  R = {R.min()/MPC*1e3:.0f} kpc - {R.max()/MPC:.1f} Mpc   N = {len(gN)}")
    print(f"       M0  no environment     chi2 = {c0:6.1f}")
    print(f"       MED medium clock       chi2 = {cm:6.1f}   b = {bm:.2f} +- {sm:.2f}")
    print(f"       NB  neighbour mass     chi2 = {cn:6.1f}   b = {bn:.2f} +- {sn:.2f}")
    # residual shape: where each model misses (outer third vs inner third)
    o = np.argsort(R); k = len(R) // 3
    for lab, mod, bb in [("MED", m_med, bm), ("NB", m_nb, bn)]:
        r = (go - mod(gN, bb, R)) / ge
        print(f"       {lab:3s} mean residual  inner {r[o[:k]].mean():+.2f}  middle {r[o[k:2*k]].mean():+.2f}  outer {r[o[2*k:]].mean():+.2f}  (sigma)")
print()
for p in ["early", "red"]:
    cm, cn = res[p]; d = cn - cm
    print(f"{p:5s}: chi2(NB) - chi2(MED) = {d:+.1f} -> " + ("MED wins" if d > 4 else "NB wins" if d < -4 else "no decision"))

print("\nPART 3 robustness")
MBAR2 = {"early": 2.666e10, "red": 1.768e10}
for r0, gm in [(4.0,1.8),(6.0,1.8),(5.0,1.6),(5.0,2.0)]:
    R0=r0*MPC; GAM=gm
    for p in ["early","red"]:
        gN,go,ge=h.load_morph(FILES[p]); R=np.sqrt(G*MBAR2[p]*MSUN/gN)
        _,_,cm=best(m_med,gN,go,ge,R); bn,_,cn=best(m_nb,gN,go,ge,R)
        print(f"r0={r0} gamma={gm} {p:5s}: MED {cm:6.1f}  NB {cn:6.1f} (b={bn:.2f})")
R0=5.0*MPC; GAM=1.8
print("cored medium: environment = neighbours beyond r_core only, enclosed mean density")
for rc in [0.1,0.3,0.6]:
    def m_med_c(gN,b,R,rc=rc*MPC):
        x=np.where(R>rc, xibar(R)*(1-(rc/R)**(3-GAM)), 0.0)   # mass between rc and R, over volume to R
        ratio=RHO_M*(1+b*x)/RHO_VAC
        return curve(gN,GDAG*np.sqrt(np.maximum(1.0,ratio)))
    for p in ["early","red"]:
        gN,go,ge=h.load_morph(FILES[p]); R=np.sqrt(G*MBAR2[p]*MSUN/gN)
        bm,_,cm=best(m_med_c,gN,go,ge,R); _,_,cn=best(m_nb,gN,go,ge,R)
        print(f"  r_core={rc} Mpc {p:5s}: MED {cm:6.1f} (b={bm:.2f})  NB {cn:6.1f}")
