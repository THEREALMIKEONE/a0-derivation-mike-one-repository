# Map 1: the edge, and which layer is active at each moment.
# Mike One, 30 October 2026.
# Edge = the universe's size, taken as the Hubble radius c/H (HMG: c/H = ct exactly).
# Active layer = the thermal layer, S_T = ln(hbar c / (k T) / R_floor): the size of a typical particle at that temperature.
# Two histories:
#   HMG: a ~ t at all times, H = 1/t, T = T0 t0 / t  (Monjo 2024, arXiv:2408.08339: linear except the Planck epoch)
#   Standard (LCDM with radiation, Planck 2018; the family used in paper 5, Section 7, with p ~ 0)
# Usage: python map1_edge_and_heat.py -> map1_edge_and_heat.txt, fig_map1.pdf/.png
import math
import numpy as np
from scipy.integrate import quad

HBARC_eVm = 1.97327e-7; kB = 8.617333e-5; c = 2.99792458e8
LP = 1.616255e-35; RF = math.sqrt(2) * LP
T0 = 2.7255; H0 = 67.4e3 / 3.0857e22; Om = 0.315; Orad = 9.15e-5; OL = 1 - Om - Orad
GYR = 3.156e16
S_of_R = lambda R: math.log(R / RF)
S_T = lambda T: S_of_R(HBARC_eVm / (kB * T))

def H_std(z): return H0 * math.sqrt(Orad * (1 + z) ** 4 + Om * (1 + z) ** 3 + OL)
def t_std(z):
    a = 1 / (1 + z)
    f = lambda la: math.exp(la) / (math.exp(la) * H0 * math.sqrt(Orad * math.exp(-4 * la) + Om * math.exp(-3 * la) + OL))
    # t(a) = integral_0^a da'/(a' H) = integral over ln a' of 1/H; radiation era below a = 1e-12 done analytically
    lo = math.log(1e-12); t_lo = (1e-12) ** 2 / (2 * H0 * math.sqrt(Orad))
    if a <= 1e-12: return a ** 2 / (2 * H0 * math.sqrt(Orad))
    return t_lo + quad(lambda la: 1 / (H0 * math.sqrt(Orad * math.exp(-4 * la) + Om * math.exp(-3 * la) + OL)), lo, math.log(a), limit=400)[0]
t0_std = t_std(0); t0_hmg = 1 / H0          # HMG age = 1/H0

events = [("today", 0.0), ("galaxy discs at z = 2", 2.0), ("CMB (recombination)", 1090.0),
          ("nucleosynthesis, T = 0.1 MeV", 0.1e6 / (kB * T0)), ("protons form, T = 150 MeV", 150e6 / (kB * T0)),
          ("electroweak, T = 100 GeV", 100e9 / (kB * T0))]
lines = [f"{'moment':32s} {'S_heat':>7s} | {'standard: t':>12s} {'S_edge':>7s} {'gap':>6s} | {'HMG: t':>12s} {'S_edge':>7s} {'gap':>6s}"]
def fmt_t(sec):
    for u, d in [("Gyr", GYR), ("yr", GYR / 1e9), ("day", 86400), ("s", 1)]:
        if sec >= d or u == "s": return f"{sec / d:9.3g} {u}"
for name, z in events:
    T = T0 * (1 + z); sT = S_T(T)
    ts = t_std(z); se = S_of_R(c / H_std(z))
    th = t0_hmg / (1 + z); sh = S_of_R(c * th)
    lines.append(f"{name:32s} {sT:7.1f} | {fmt_t(ts):>12s} {se:7.1f} {se - sT:6.1f} | {fmt_t(th):>12s} {sh:7.1f} {sh - sT:6.1f}")
lines.append("")
lines.append("gap = edge minus active layer, in e-folds. HMG: constant, the active layer always 67 e-folds below the edge.")
lines.append("Standard: the gap shrinks into the past (radiation era: H ~ T^2).")
lines.append(f"Today both edges agree: standard c/H0 at S = {S_of_R(c / H0):.2f}, HMG ct0 at S = {S_of_R(c * t0_hmg):.2f}.")
z2 = 2.0
lines.append(f"At z = 2 the edges differ by {S_of_R(c / H_std(z2)) - S_of_R(c * t0_hmg / (1 + z2)):+.2f} e-folds: the late sheet is the same in both.")
open("map1_edge_and_heat.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
zs = np.logspace(-3, 15, 400)
fig, ax = plt.subplots(figsize=(7.2, 3.9))
ts = np.array([t_std(z) for z in zs]); th = t0_hmg / (1 + zs)
ax.plot([S_of_R(c / H_std(z)) for z in zs], np.log10(ts / GYR), color="#1f5a99", lw=2.2, label="edge ($c/H$; $=ct$ at late times)")
ax.plot([S_T(T0 * (1 + z)) for z in zs], np.log10(ts / GYR), color="#b3432b", lw=2.2, label="active layer ($\\hbar c/kT$)")
for name, z in events[2:5]:
    T = T0 * (1 + z)
    ax.plot([S_T(T)], [math.log10(t_std(z) / GYR)], "o", color="#b3432b", ms=4)
    ax.annotate(name.split(",")[0], (S_T(T), math.log10(t_std(z) / GYR)), xytext=(-8, 4), textcoords="offset points",
                fontsize=6.8, ha="right")
ax.set_xlabel("$S$"); ax.set_ylabel("$\\log_{10}$ cosmic time (Gyr)")
ax.set_xlim(30, 145); ax.legend(loc="lower right", fontsize=7, frameon=False)
ax.set_title("map 1: the edge and the active layer through time", fontsize=9, loc="left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("fig_map1.pdf"); fig.savefig("fig_map1.png", dpi=160)
