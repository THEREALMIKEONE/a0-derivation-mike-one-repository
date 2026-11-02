# The slice + stack model on the S axis.
# Mike One, 29 October 2026. Companion to "Hyperconical gravity applied to the renormalization scale", Section 2.
#
# Every layer of S is a slice with HMG's geometry: radius R, clock t = R/c,
# vacuum Omega_L * 3 / (8 pi G t^2).
#   Gravity: the projection of one slice. Coupling (l_P/R)^2 at every layer.
#   Gauge interactions: the projection of the stack. Each e-fold of S adds the medium of the
#   species active there; gauge-boson media anti-screen (HMG's neighbourhood branch), matter
#   media screen (HMG's Hubble branch). Read along S this is the renormalisation-group equation,
#   since S grows with size: d(1/alpha)/dS = +b/(2 pi) at one loop.
# Inputs: one baseline per coupling at the Z mass. Everything else follows from the species
# active in each stretch of S. Strong medium at two loops, electromagnetic medium at one loop.
# Pass rule (fixed before running): the stack's values lie within the measured errors.
#
# Usage: python stack_model.py      (no data needed; writes stack_model_rows.json)
import math, json

HBARC = 1.97327e-16          # GeV m
L_P = 1.616255e-35           # m
R_FLOOR = math.sqrt(2) * L_P  # where the quantum and black-hole lines meet
C = 2.99792458e8
G = 6.6743e-11
OMEGA_L = 0.6937             # HMG's projected vacuum fraction
MZ = 91.1876
L_CONF = 0.286               # GeV, where alpha_s reaches 1: below it the strong medium is bound into hadrons

def S_of_mu(mu): return math.log(HBARC / mu / R_FLOOR)
def mu_of_S(S): return HBARC / (R_FLOOR * math.exp(S))

QUARKS = {'t': (172.7, 2/3), 'b': (4.18, -1/3), 'c': (1.27, 2/3),
          's': (0.095, -1/3), 'd': (0.0047, -1/3), 'u': (0.0022, 2/3)}
LEPTONS = {'tau': 1.777, 'mu': 0.10566, 'e': 0.000511}

def medium(mu):
    """Species active in the medium at scale mu (GeV)."""
    q = [k for k, (m, _) in QUARKS.items() if mu > m and mu > L_CONF]
    l = [k for k, m in LEPTONS.items() if mu > m]
    return q, l

def dinv_alpha_s(mu, a):
    """d(1/alpha_s)/dS from the strong medium, two loops."""
    nf = len(medium(mu)[0])
    b0 = (33 - 2 * nf) / (12 * math.pi)
    b1 = (153 - 19 * nf) / (24 * math.pi ** 2)
    return -2 * (b0 + b1 * a)

def dinv_alpha_em(mu, a):
    """d(1/alpha_em)/dS from the charged medium below M_Z, one loop."""
    q, l = medium(mu)
    b = 4 / 3 * (sum(3 * QUARKS[k][1] ** 2 for k in q) + len(l))
    return b / (2 * math.pi)

def stack(a0, S0, S1, rate, n=20000):
    """Sum the medium layer by layer from S0 to S1 (midpoint rule)."""
    inv, h = 1 / a0, (S1 - S0) / n
    for i in range(n):
        S = S0 + (i + 0.5) * h
        inv += rate(mu_of_S(S), 1 / inv) * h
        if inv <= 0: return None
    return 1 / inv

SZ = S_of_mu(MZ)
AS_Z, AEM_Z = 0.1180, 1 / 127.951     # MSbar baselines at M_Z

layers = [('Z mass (baseline)', SZ), ('bottom quark', S_of_mu(4.18)), ('tau', S_of_mu(1.777)),
          ('charm', S_of_mu(1.27)), ('confinement', S_of_mu(L_CONF)), ('muon', S_of_mu(0.10566)),
          ('electron', S_of_mu(0.000511)), ('Bohr radius', math.log(5.29e-11 / R_FLOOR)),
          ('largest cold body', math.log(5.88e7 / R_FLOOR)), ('Milky Way', math.log(1.5e21 / R_FLOOR)),
          ('our layer', math.log(C / (67.4e3 / 3.0857e22) / R_FLOOR))]

def fmt(x): return '   -   ' if x is None else (f'{x:.4f}' if x > 1e-3 else f'{x:.2e}')

rows = []
print(f"{'layer':20s} {'S':>6s} {'radius m':>9s} {'slice vac':>10s} {'gravity':>9s} {'alpha_s':>8s} {'1/a_em':>8s}")
for name, S in layers:
    R = R_FLOOR * math.exp(S); t = R / C
    rho = OMEGA_L * 3 / (8 * math.pi * G * t * t)
    grav = (L_P / R) ** 2
    mu = mu_of_S(S)
    a_s = stack(AS_Z, SZ, S, dinv_alpha_s) if mu > L_CONF else None   # bound into hadrons below
    S_e = S_of_mu(LEPTONS['e'])
    a_em = stack(AEM_Z, SZ, min(S, S_e), dinv_alpha_em)                # frozen below the electron
    rows.append(dict(layer=name, S=round(S, 1), R_m=R, slice_vacuum=rho, gravity=grav,
                     alpha_s=a_s, alpha_em=a_em))
    print(f"{name:20s} {S:6.1f} {R:9.2e} {rho:10.2e} {grav:9.2e} {fmt(a_s):>8s} "
          f"{(1 / a_em):8.2f}")
json.dump(rows, open('stack_model_rows.json', 'w'), indent=1)

print("\nchecks (pass rule: within the measured errors):")
amb = stack(AS_Z, SZ, S_of_mu(4.18), dinv_alpha_s)
amt = stack(AS_Z, SZ, S_of_mu(1.777), dinv_alpha_s)
aem_e = stack(AEM_Z, SZ, S_of_mu(0.000511), dinv_alpha_em)
print(f" alpha_s at the bottom-quark scale: stack {amb:.3f}, measured 0.225 +- 0.010 -> {'pass' if abs(amb - 0.225) <= 0.010 else 'fail'}")
print(f" alpha_s at the tau scale:          stack {amt:.3f}, measured 0.325 +- 0.015 -> {'pass' if abs(amt - 0.325) <= 0.015 else 'fail'}")
print(f" 1/alpha_em at the electron:        stack {1 / aem_e:.2f}, measured 137.036"
      f" ({100 * (137.036 - 1 / aem_e) / 137.036:.2f}% short: the hadronic part, which needs data)")
