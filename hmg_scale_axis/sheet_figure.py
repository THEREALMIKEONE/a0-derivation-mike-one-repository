# The (S, t) sheet of "Hyperconical gravity applied to the renormalization scale", Figure 1.
# Mike One, 30 October 2026.
#   x: S = ln(R / R_floor), R_floor = sqrt(2) l_P.   y: cosmic time t.
#   Edge of the universe, R = c t (HMG's slice).  Today's cut: the layers of Table 1 at t0.
#   Right panel: the last 13 Gyr near our layer. Dark-energy layer from the measured p = 0.011
#   (S_vac(t) = S_vac,0 - p ln(H/H0), with H = 1/t on the edge), against p = 1 (vacuum carried along the edge).
# Usage: python sheet_figure.py  -> fig_sheet.pdf, fig_sheet.png
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

c = 2.99792458e8; lP = 1.616255e-35; Rf = math.sqrt(2) * lP
tP = 5.391e-44; GYR = 3.156e16; t0 = 13.8 * GYR
BLUE, RED, GREY = "#1f5a99", "#b3432b", "0.35"

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.6, 3.9), gridspec_kw=dict(width_ratios=[1.6, 1]))

# left: the whole sheet
# standard history edge c/H(t): radiation + matter + vacuum (Planck 2018), times from the scale factor
H0 = 67.4e3 / 3.0857e22; Om, Orad = 0.315, 9.15e-5; OLm = 1 - Om - Orad
from scipy.integrate import quad
def H_a(a): return H0 * math.sqrt(Orad * a ** -4 + Om * a ** -3 + OLm)
def t_a(a):
    a1 = 1e-12
    if a <= a1: return a * a / (2 * H0 * math.sqrt(Orad))
    return a1 * a1 / (2 * H0 * math.sqrt(Orad)) + quad(lambda la: 1 / H_a(math.exp(la)), math.log(a1), math.log(a), limit=400)[0]
avals = np.logspace(-31, 0, 500)
t = np.array([t_a(a) for a in avals]); S_edge = np.log(c / np.array([H_a(a) for a in avals]) / Rf)
t_hmg = np.logspace(math.log10(tP), math.log10(1 / H0), 400); S_hmg = np.log(c * t_hmg / Rf)   # seconds; HMG age = 1/H0
ax.fill_betweenx(np.log10(t / GYR), -2, S_edge, color="#e3ebf5", lw=0)
ax.plot(S_edge, np.log10(t / GYR), color=BLUE, lw=2.2)
ax.axhline(math.log10(t0 / GYR), color=GREY, lw=1.0, ls="--")
layers = [("floor", 0.0, 0), ("W", 39.2, 0), ("proton", 45.1, 1), ("electron", 51.2, 0), ("atoms", 56.1, 1),
          ("Sun", 100.1, 0), ("Milky Way", 128.5, 1), ("ours", 139.9, 0)]
for name, S, up in layers:
    ax.plot([S], [math.log10(t0 / GYR)], "o", color="0.15", ms=3.5, zorder=5)
    ax.annotate(name, (S, math.log10(t0 / GYR)), xytext=(0, 5 + 9 * up), textcoords="offset points",
                ha="center", fontsize=7)
ax.text(68, -32, "everything below the edge\nexists at that time", fontsize=7.5, color=BLUE, ha="center")
ax.text(112, -44, "edge of the universe, $c/H$\n($=ct$, HMG's slice, at late times)", fontsize=7.5, color=BLUE, ha="center")
ax.text(2, -4.5, "today's cut across all layers", fontsize=7, color=GREY)
ax.set_xlim(-2, 146); ax.set_ylim(math.log10(tP / GYR), 6)
ax.set_xlabel("$S$ (e-folds of radius above the floor)")
ax.set_ylabel("$\\log_{10}$ cosmic time (Gyr)")
ax.set_title("the sheet", fontsize=9)

# right: the crossing, last 13 Gyr, linear time
tt = np.linspace(1.0, 13.8, 300) * GYR
S_e = np.log(c * tt / Rf)
p = 0.011; Svac0 = 139.92
S_v = Svac0 + p * np.log(tt / t0)          # H = 1/t on the edge: -p ln(H/H0) = p ln(t/t0)
bx.plot(S_e, tt / GYR, color=BLUE, lw=2.2, label="edge, $R=ct$")
bx.plot(S_v, tt / GYR, color=RED, lw=2.2, label="dark-energy layer, $p=0.011$")
bx.plot(S_e + (Svac0 - math.log(c * t0 / Rf)), tt / GYR, color=RED, lw=1.2, ls=(0, (2, 2)),
        label="$p=1$: carried by the edge")
bx.plot([139.90], [13.8], "*", color=RED, ms=12, zorder=6)
bx.annotate("why now", (139.90, 13.8), xytext=(4, -24), textcoords="offset points", fontsize=7.5)
bx.set_xlim(137.4, 140.6); bx.set_ylim(0, 15.5)
bx.set_xlabel("$S$"); bx.set_ylabel("cosmic time (Gyr)")
bx.set_title("near our layer", fontsize=9)
bx.legend(loc="center left", fontsize=7, frameon=False)
for a in (ax, bx):
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("fig_sheet.pdf"); fig.savefig("fig_sheet.png", dpi=160)
print(f"edge today S = {math.log(c * t0 / Rf):.2f}; edge at Planck time S = {math.log(c * tP / Rf):.2f}; "
      f"dark-energy layer moved {p * math.log(1 / 13.8):+.3f} in S since t = 1 Gyr")
