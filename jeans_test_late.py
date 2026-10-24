# ============================================================
# Test 2 of g† = c·sqrt(4πGρ_Λ)/2π on KiDS-1000 late-type and blue isolated lenses (Brouwer+2021 Fig. 8)
# Rules FROZEN 2026-10-24 before running:
#   Late/blue only (no hot halo, clean baryons). Masses: M* late 10^9.745, blue 10^9.610 + cold gas (Boselli 2014).
#   g_bar = given per point; model g = g_bar/(1−exp(−sqrt(g_bar/g†))); log residuals with ESD errors; points with g_obs>0.
#   Lens redshift for the expansion candidate: z = 0.25 (KiDS-bright lenses, 0.1<z<0.5).
#   Candidates: J = 1.056e-10 (constant); E = cH(0.25)/2π; L3 = 0.863e-10 (constant).
#   TEST: pass if Δχ² < 4 from the best-fit single g† (late and blue fitted jointly).
#   Caveat stated in advance: KiDS lenses overlap with Mistele+2024's, so this is not fully independent of test 1.
# ============================================================
import sys, math, numpy as np
D=sys.argv[1]; ESD2G=4.0*4.52e-30*3.086e16
c=299792458.0; H0=67.4e3/3.0856775814913673e22; Om=0.315; OL=1-Om
pts=[]
for f in ("Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt","Fig-8_RAR-KiDS-isolated_Colorbin_1.txt"):
    d=np.loadtxt(f"{D}/{f}",comments="#"); gb=d[:,0]; go=ESD2G*d[:,1]/d[:,4]; ge=ESD2G*d[:,3]/d[:,4]; ok=go>0
    pts.append((gb[ok],go[ok],ge[ok]/(go[ok]*math.log(10))))
def chi(gd):
    s=0.0
    for gb,go,sg in pts: s+=np.sum(((np.log10(gb/(1-np.exp(-np.sqrt(gb/gd))))-np.log10(go))/sg)**2)
    return float(s)
lg=np.linspace(-10.6,-9.4,2401); cs=np.array([chi(10**x) for x in lg]); b=cs.argmin(); m=cs[b]
one=10**lg[cs<m+1]; two=10**lg[cs<m+4]; n=sum(len(p[0]) for p in pts)
print(f"best-fit g† = {10**lg[b]:.3e} (1σ {one.min():.3e}–{one.max():.3e}; 2σ {two.min():.3e}–{two.max():.3e}); χ²min = {m:.1f} / {n-1} dof")
J=c*H0*math.sqrt(1.5*OL)/(2*math.pi); L3=c*H0*math.sqrt(OL)/(2*math.pi); E=c*H0*math.sqrt(Om*1.25**3+OL)/(2*math.pi)
for name,val in [("J  Jeans, constant",J),("E  expansion at z=0.25",E),("L3 Λ/3, constant",L3)]:
    dc=chi(val)-m; print(f"  {name:24s} g†={val:.3e}  Δχ² = {dc:+5.1f}  {'PASS' if dc<4 else 'FAIL'}")
