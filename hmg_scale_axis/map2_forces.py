# Map 2: one frame of the sheet. All four interactions across S, from the floor to our layer, today.
# Mike One, 30 October 2026.
#
# S = ln(R / R_floor), R_floor = sqrt(2) l_P; a layer's energy scale is mu = hbar c / R.
# Gravity: alpha_G = (mu / M_P)^2  (one slice; equals the slice vacuum over the Planck density).
# Gauge interactions: the stack of media below each layer, d(1/alpha_i)/dS = +b_i(S) / 2 pi
# (S grows with size, so this is the renormalisation-group equation read along S).
#   Above the Z mass: alpha_1 (GUT-normalised), alpha_2, alpha_3 with Standard Model coefficients, top threshold.
#   Below the Z: alpha_s two loops with quark thresholds until confinement (alpha_s -> bound into hadrons);
#   alpha_em one loop with charged-lepton and quark thresholds, frozen below the electron.
#   Weak interaction below the W mass: short range, effective strength alpha_2 (mu / M_W)^2.
# Baselines at the Z mass (PDG): alpha_s = 0.1180, 1/alpha_em = 127.951, sin^2 theta_W = 0.23122.
# Usage: python map2_forces.py   -> map2_forces.txt, fig_map2_forces.pdf/.png  (no data needed)
import math, json
import numpy as np

HBARC = 1.97327e-16; LP = 1.616255e-35; RF = math.sqrt(2) * LP; MP = 1.220890e19   # GeV
MZ, MW, MT = 91.1876, 80.377, 172.69
LCONF = 0.286
S_of_mu = lambda mu: math.log(HBARC / mu / RF)
mu_of_S = lambda S: HBARC / (RF * math.exp(S))
c = 2.99792458e8; H0 = 67.4e3 / 3.0857e22
S_OURS = math.log(c / H0 / RF)

QUARKS = {'t': (MT, 2/3), 'b': (4.18, -1/3), 'c': (1.27, 2/3), 's': (0.095, -1/3), 'd': (0.0047, -1/3), 'u': (0.0022, 2/3)}
LEPTONS = {'tau': 1.77686, 'mu': 0.105658, 'e': 0.000511}

aem_Z = 1 / 127.951; s2w = 0.23122; as_Z = 0.1180
a1_Z = (5/3) * aem_Z / (1 - s2w)       # GUT-normalised hypercharge
a2_Z = aem_Z / s2w

def b_above(mu):
    """One-loop SM coefficients (b1 GUT-normalised, b2, b3) above M_Z, top included above its mass."""
    if mu > MT: return 41/10, -19/6, -7.0
    # below the top: remove the top quark's contributions (b1: 17/30, b2: 1/2... standard 5-flavour values)
    return 41/10 - 17/30, -19/6 - 1/2, -7.0 - 2/3

def nf(mu): return sum(1 for m, _ in QUARKS.values() if mu > m and mu > LCONF)

def run_up(a0, which, S1, n=20000):
    """From M_Z toward the floor (decreasing S)."""
    S0 = S_of_mu(MZ); inv = 1 / a0; h = (S1 - S0) / n
    for i in range(n):
        mu = mu_of_S(S0 + (i + 0.5) * h)
        b = b_above(mu)[which]
        inv += b / (2 * math.pi) * h
    return 1 / inv

def run_down_s(S1, n=20000):
    """alpha_s from M_Z toward larger S, two loops, until confinement."""
    S0 = S_of_mu(MZ); a = as_Z; h = (S1 - S0) / n
    for i in range(n):
        mu = mu_of_S(S0 + (i + 0.5) * h)
        if mu < LCONF: return None
        f = nf(mu); b0 = (33 - 2 * f) / (12 * math.pi); b1 = (153 - 19 * f) / (24 * math.pi ** 2)
        a += 2 * a * a * (b0 + b1 * a) * h          # d alpha/dS = -d alpha/d ln mu
        if a > 5: return None
    return a

def run_down_em(S1, n=20000):
    S0 = S_of_mu(MZ); inv = 1 / aem_Z; S1 = min(S1, S_of_mu(LEPTONS['e'])); h = (S1 - S0) / n
    for i in range(n):
        mu = mu_of_S(S0 + (i + 0.5) * h)
        q = [k for k, (m, Q) in QUARKS.items() if mu > m and mu > LCONF]
        l = [k for k, m in LEPTONS.items() if mu > m]
        b = 4 / 3 * (sum(3 * QUARKS[k][1] ** 2 for k in q) + len(l))
        inv += b / (2 * math.pi) * h
    return 1 / inv

SZ = S_of_mu(MZ)
def frame(S):
    mu = mu_of_S(S)
    out = dict(S=S, mu_GeV=mu, R_m=RF * math.exp(S), gravity=(mu / MP) ** 2)
    if S <= SZ:
        a1, a2, a3 = run_up(a1_Z, 0, S), run_up(a2_Z, 1, S), run_up(as_Z, 2, S)
        out.update(strong=a3, weak=a2, hypercharge=a1, em=1 / ((5 / 3) / a1 + 1 / a2))
    else:
        out.update(strong=run_down_s(S), em=run_down_em(S),
                   weak=a2_Z * (mu / MW) ** 2 if mu < MW else a2_Z, hypercharge=None)
    return out

if __name__ == "__main__":
    grid = np.linspace(0, S_OURS, 1401)
    rows = [frame(S) for S in grid]
    json.dump(rows, open("map2_forces_grid.json", "w"))
    named = [("floor", 0.0), ("alpha_2 = alpha_3", None), ("alpha_1 = alpha_2", None), ("top quark", S_of_mu(MT)),
             ("Z mass (baseline)", SZ), ("W mass", S_of_mu(MW)), ("bottom quark", S_of_mu(4.18)), ("tau", S_of_mu(1.77686)),
             ("confinement (proton at 45.1)", S_of_mu(LCONF)), ("muon", S_of_mu(0.105658)), ("electron", S_of_mu(0.000511)),
             ("Bohr radius (atoms)", math.log(5.29e-11 / RF)), ("largest cold body", math.log(5.88e7 / RF)),
             ("Sun", math.log(6.96e8 / RF)), ("Milky Way with halo", math.log(1.5e21 / RF)), ("our layer, c/H0", S_OURS)]
    # crossings above the Z
    def cross(i, j):
        f = lambda S: (lambda r: (1 / r[i]) - (1 / r[j]))(frame(S))
        lo, hi = 0.0, SZ
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(lo) * f(mid) <= 0: hi = mid
            else: lo = mid
        return (lo + hi) / 2
    S23 = cross("weak", "strong"); S12 = cross("hypercharge", "weak")
    def gcross(key):
        lo, hi = 0.0, SZ
        g = lambda S: math.log(frame(S)["gravity"]) - math.log(frame(S)[key])
        for _ in range(60):
            mid = (lo + hi) / 2
            if g(lo) * g(mid) <= 0: hi = mid
            else: lo = mid
        return (lo + hi) / 2
    Sg_em = gcross("em")
    named[1] = ("alpha_2 = alpha_3", S23); named[2] = ("alpha_1 = alpha_2", S12)
    named.insert(1, ("gravity = electromagnetism", Sg_em))
    f = lambda x: "    -    " if x is None else (f"{x:9.4f}" if x > 1e-3 else f"{x:9.2e}")
    lines = [f"{'layer':30s} {'S':>6s} {'mu (GeV)':>10s} {'gravity':>9s} {'strong':>9s} {'weak':>9s} {'1/a_em':>8s}"]
    for name, S in named:
        r = frame(S)
        lines.append(f"{name:30s} {S:6.1f} {r['mu_GeV']:10.3e} {f(r['gravity'])} {f(r['strong'])} {f(r['weak'])} "
                     f"{(1 / r['em']):8.2f}")
    r0 = frame(0.0)
    lines.append("")
    lines.append(f"at the floor: gravity {r0['gravity']:.3f}, strong {r0['strong']:.4f}, weak {r0['weak']:.4f}, "
                 f"hypercharge (GUT) {r0['hypercharge']:.4f}: the three gauge couplings do not meet in one point")
    lines.append(f"checks: alpha_s(m_b) = {frame(S_of_mu(4.18))['strong']:.3f} (0.225 +- 0.010), "
                 f"alpha_s(m_tau) = {frame(S_of_mu(1.77686))['strong']:.3f} (0.325 +- 0.015), "
                 f"1/alpha_em(electron) = {1 / frame(S_of_mu(0.000511))['em']:.2f} (137.036; hadronic part from data)")
    open("map2_forces.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, (bx, ax) = plt.subplots(2, 1, figsize=(7.4, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[0.55, 4]))
    # what holds the structures of each stretch of S (Table 1 of paper 5 / v3)
    # the four force layers (Mike One, 1 Nov 2026); sign flips at the floor and at our layer
    S_E = math.log(3.86e-13 / RF)
    bands = [(0, S_of_mu(MW), "#d9cfe9", "weak"),
             (S_of_mu(MW), S_E, "#cfdcef", "strong"),
             (S_E, math.log(5.88e7 / RF), "#bfe0cc", "electromagnetic: atoms, chemistry, life"),
             (math.log(5.88e7 / RF), S_OURS, "#ddd0ee", "gravity: planets, stars, galaxies")]
    for a, b, colr, lab in bands:
        bx.axvspan(a, b, color=colr, lw=0)
        if lab: bx.text((a + b) / 2, 0.5, lab, ha="center", va="center", fontsize=6.2)
    for xs in (0, S_OURS): bx.axvline(xs, color="#b3432b", lw=1.6)
    bx.text(1.2, 0.12, "sign flip", fontsize=5.6, color="#b3432b", ha="left", va="bottom")
    bx.text(S_OURS - 1.2, 0.12, "sign flip", fontsize=5.6, color="#b3432b", ha="right", va="bottom")
    bx.set_yticks([]); bx.set_ylim(0, 1)
    bx.set_title("the four force layers, each with its own vacuum", fontsize=8, loc="left")
    for s in ("top", "right", "left"): bx.spines[s].set_visible(False)
    S = np.array([r["S"] for r in rows])
    def col(k): return np.array([np.nan if (r[k] is None or r[k] <= 0) else r[k] for r in rows])
    ax.plot(S, np.log10(col("gravity")), color="#7a4fb5", lw=2, label="gravity, one slice")
    ax.plot(S, np.log10(col("strong")), color="#1f5a99", lw=2, label="strong")
    ax.plot(S, np.log10(col("weak")), color="#c27c0e", lw=2, label="weak")
    ax.plot(S, np.log10(col("em")), color="#2e8b57", lw=2, label="electromagnetic")
    for x, lab in [(S_of_mu(MW), "W"), (S_of_mu(LCONF), "protons"), (S_of_mu(0.000511), "electron"),
                   (math.log(5.29e-11 / RF), "atoms"), (math.log(5.88e7 / RF), "planets"),
                   (math.log(6.96e8 / RF), "stars"), (math.log(1.5e21 / RF), "galaxies"), (S_OURS, "ours")]:
        ax.axvline(x, color="0.8", lw=0.8, zorder=0)
        ax.text(x, 0.6, lab, rotation=90, fontsize=6.5, ha="right", va="bottom", color="0.35")
    ax.set_xlim(0, 141); ax.set_ylim(-12, 2.5)
    ax.set_xlabel("$S$ (e-folds of radius above the floor)"); ax.set_ylabel("$\\log_{10}$ coupling")
    ax.legend(loc="lower left", fontsize=7.5, frameon=False, ncol=2)
    ax.text(92, -6.5, "gravity keeps falling: $10^{-45}$ at the electron, $10^{-122}$ at our layer", fontsize=7, color="#7a4fb5", ha="center")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig("fig_map2_forces.pdf"); fig.savefig("fig_map2_forces.png", dpi=160)
