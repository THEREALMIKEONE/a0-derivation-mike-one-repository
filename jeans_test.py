# ============================================================
# Independent test of g† = c·sqrt(4πGρ_Λ)/2π  (Jeans rate of the vacuum medium)
# Rules FROZEN 2026-10-24 before running:
#   Data: Mistele+2024 galaxy-galaxy weak-lensing circular velocities (Table 1, via Monjo's repo), 20 points.
#   Model: g_obs = g_bar / (1 − exp(−sqrt(g_bar/g†))), g_bar = G M_bar / R², fitted in v with σ_v (std ⊕ sys).
#   Candidates (no parameter fitted to these data):
#     J  Jeans rate of the vacuum, constant:          c·sqrt(4πGρ_Λ)/2π = 1.056e-10
#     E  expansion rate at each bin's redshift:        c·H(z)/2π
#     L3 standard vacuum form, constant:               c²·sqrt(Λ/3)/2π = 0.863e-10
#   TEST: a candidate passes if its χ² is within Δχ² < 4 (2σ) of the best-fit single g†.
#   Known caveat stated in advance: stacks mix early and late types; early types show excess lensing (KiDS),
#   which biases the fitted g† high.
# ============================================================
import sys, math, numpy as np
d = np.loadtxt(sys.argv[1], comments="#")
z, lMb, lR, v, es, ey = d[:,1], d[:,2], d[:,3], d[:,4], d[:,5], d[:,6]
G=6.674e-11; MS=1.989e30; KPC=3.0857e19; c=299792458.0
H0=67.4e3/3.0856775814913673e22; Om=0.315; OL=1-Om
R=10**lR*KPC; gb=G*10**lMb*MS/R**2; sv=np.hypot(es,ey)*1e3
def vpred(gd): g=gb/(1-np.exp(-np.sqrt(gb/gd))); return np.sqrt(g*R)
chi=lambda gd: float(np.sum(((vpred(gd)-v*1e3)/sv)**2))
lg=np.linspace(-10.6,-9.4,2401); cs=np.array([chi(10**x) for x in lg]); b=cs.argmin()
best=10**lg[b]; m=cs[b]; one=10**lg[cs<m+1]; two=10**lg[cs<m+4]
print(f"best-fit g† = {best:.3e}  (1σ {one.min():.3e}–{one.max():.3e}; 2σ {two.min():.3e}–{two.max():.3e})  χ²min = {m:.1f} / 19 dof")
J=c*H0*math.sqrt(1.5*OL)/(2*math.pi); L3=c*H0*math.sqrt(OL)/(2*math.pi)
E=lambda zz: c*H0*math.sqrt(Om*(1+zz)**3+OL)/(2*math.pi)
gE=np.array([E(x) for x in z])
def chi_arr(gdarr):
    g=gb/(1-np.exp(-np.sqrt(gb/gdarr))); return float(np.sum(((np.sqrt(g*R)-v*1e3)/sv)**2))
for name,val in [("J  Jeans, constant", chi(J)), ("E  expansion H(z)", chi_arr(gE)), ("L3 Λ/3, constant", chi(L3))]:
    dc=val-m; print(f"  {name:22s} χ² = {val:6.1f}  Δχ² = {dc:+5.1f}  {'PASS' if dc<4 else 'FAIL'}")
print(f"  (values: J={J:.3e}, L3={L3:.3e}, E ranges {gE.min():.3e}–{gE.max():.3e} over z={z.min()}–{z.max()})")
