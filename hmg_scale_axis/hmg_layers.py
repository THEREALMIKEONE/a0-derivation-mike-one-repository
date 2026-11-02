# HMG's slice formulas on every layer of the scale axis, and the checks at our layer.
# S = ln(R/R_floor), R_floor = sqrt(2) l_P (v3 convention). A slice of radius R has clock t = R/c and
# vacuum rho = Omega_L * 3/(8 pi G t^2) with HMG's Omega_L = 0.6937 (Monjo 2018, PRD 98, 043508).
# Usage: python hmg_layers.py
import math
c=2.99792458e8; G=6.67430e-11; hbar=1.054571817e-34; lP=math.sqrt(hbar*G/c**3); Rf=math.sqrt(2)*lP
OL_HMG=0.6937; OL_FIT=0.6972; sOL=0.0037         # OL_FIT: vacuum_fraction_profile.py (Pantheon+, DESI DR2, Planck priors)
H0=67.4e3/3.0857e22
S=lambda R: math.log(R/Rf)
rows=[("floor",Rf),("W mass",hbar*c/(80.4e9*1.602e-19)),("confinement (alpha_s=1)",hbar*c/(0.286e9*1.602e-19)),
      ("electron",hbar*c/(0.511e6*1.602e-19)),("Bohr radius",5.29e-11),("largest cold body",5.88e7),
      ("Milky Way with halo",1.5e21),("our layer, c/H0",c/H0)]
print(f"{'layer':26s} {'S':>6s} {'radius m':>10s} {'vacuum kg/m3':>13s} {'accel m/s2':>11s} {'rho/rho_P':>10s} {'(l_P/R)^2':>10s}")
rhoP=OL_HMG*3/(8*math.pi*G*(lP/c)**2)
for n,R in rows:
    t=R/c; rho=OL_HMG*3/(8*math.pi*G*t**2); a=c*math.sqrt(4*math.pi*G*rho)/(2*math.pi)
    print(f"{n:26s} {S(R):6.1f} {R:10.2e} {rho:13.2e} {a:11.2e} {rho/rhoP:10.2e} {(lP/R)**2:10.2e}")
print(f"\nour layer: HMG vacuum = {OL_HMG:.4f} rho_crit; fitted {OL_FIT} +- {sOL} -> {(OL_FIT-OL_HMG)/sOL:.1f} sigma")
w=math.sqrt(1.5*OL_HMG); g=c*w*H0/(2*math.pi); g0=4*math.pi/w
from math import cos, pi
lo,hi=1.0,1.56
for _ in range(80):
    m=(lo+hi)/2; lo,hi=(m,hi) if m/cos(m)<g0 else (lo,m)
print(f"galactic scale g = c*sqrt(4 pi G rho)/2pi = {g:.3e} m/s2 ; gamma0 = 4pi/sqrt(3 Omega_L/2) = {g0:.2f} ; angle = {m/pi:.4f} pi (Monjo 2024 fits: 0.460+-0.002, 0.468+-0.003)")
