# Map 3 drawn: the film of the sheet. Every frame is map 2 at one time; only layers below the edge exist.
# Bands = what holds the structures at each S (map 2). Markers = where a layer has been checked through time.
# Mike One, 30 October 2026. Usage: python map3_film.py -> fig_map3_film.pdf/.png
import math
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

c = 2.99792458e8; LP = 1.616255e-35; RF = math.sqrt(2) * LP; HBARC = 1.97327e-16
GYR = 3.156e16; t0 = 13.8; tP = 5.391e-44 / GYR
S_mu = lambda mu: math.log(HBARC / mu / RF)
S_R = lambda R: math.log(R / RF)
S_OURS = math.log(c / (67.4e3 / 3.0857e22) / RF)

bands = [(0, S_mu(80.4), "#d9cfe9", "weak"), (S_mu(80.4), S_R(3.86e-13), "#cfdcef", "strong"),
         (S_R(3.86e-13), S_R(5.88e7), "#bfe0cc", "electromagnetic"),
         (S_R(5.88e7), S_OURS + 1.5, "#ddd0ee", "gravity")]
fig, ax = plt.subplots(figsize=(7.4, 4.3))
ly0, ly1 = math.log10(tP), math.log10(60)
for a, b, colr, lab in bands:
    ax.axvspan(a, b, color=colr, lw=0)
    xm = (a + b) / 2; ye = math.log10(RF * math.exp(xm) / c / GYR)
    narrow = (b - a) < 12
    ax.text(xm, (ye + math.log10(t0)) / 2, lab, ha="center", va="center", fontsize=6.8,
            rotation=90 if narrow else 0, zorder=1)
# above the edge: not yet inside the universe (S > S_u(t))
lt = np.linspace(ly0, math.log10(t0), 300)
Su = np.log(c * (10 ** lt) * GYR / RF)
poly = [(s, y) for s, y in zip(Su, lt)] + [(146, math.log10(t0)), (146, ly0)]
ax.add_patch(Polygon(poly, closed=True, color="white", zorder=2))
ax.plot(Su, lt, color="#1f5a99", lw=2.2, zorder=3, label="edge of the universe (map 1)")
ax.axhline(math.log10(t0), color="0.3", lw=1, ls="--", zorder=3)
ax.text(1, math.log10(t0) + 0.6, "today: one frame = map 2", fontsize=7, color="0.3")
ax.text(122, -20, "not yet inside\nthe universe", fontsize=8, color="0.45", ha="center", zorder=4)
checks = [(S_mu(0.286), math.log10(t0), "G (lunar ranging)"),
          (S_mu(0.286) + 3, math.log10(6.3), "$m_p/m_e$, 7.5 Gyr ago"),
          (S_R(5.29e-11), math.log10(5.6), "$\\alpha$, 8 Gyr ago"),
          (S_R(5.29e-11) + 1, math.log10(3.8e5 / 1e9), "$\\alpha$, $m_e$: CMB"),
          (S_OURS - 0.4, math.log10(6.0), "our vacuum layer, $p$")]
key = []
for i, (s, y, lab) in enumerate(checks, 1):
    ax.plot([s], [y], "o", color="#b3432b", ms=9, zorder=5)
    ax.text(s, y, str(i), color="white", fontsize=6, ha="center", va="center", zorder=6, fontweight="bold")
    key.append(f"{i}  {lab}")
ax.text(86, -28, "held still, checked:\n" + "\n".join(key), fontsize=7, color="#7a2a1a", va="top", zorder=4,
        linespacing=1.4)
ax.set_xlim(-1, 146); ax.set_ylim(ly0, ly1)
ax.set_xlabel("$S$ (e-folds of radius above the floor)"); ax.set_ylabel("$\\log_{10}$ cosmic time (Gyr)")
ax.legend(loc="lower right", fontsize=7.5, frameon=False)
ax.set_title("map 3: every layer holds still while the edge climbs through them", fontsize=9, loc="left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("fig_map3_film.pdf"); fig.savefig("fig_map3_film.png", dpi=160)
