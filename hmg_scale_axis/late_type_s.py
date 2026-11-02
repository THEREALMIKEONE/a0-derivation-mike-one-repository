# Your neighbourhood-branch parameter s for late-type lenses, predicted by setting HMG's g_s to the vacuum clock.
# Uses HMG's own model code (hmg_model.py from github.com/robertmonjo/lensing_RAR_KiDS).
# Usage: python late_type_s.py <path to lensing_RAR_KiDS checkout>
import sys; sys.path.insert(0, sys.argv[1]); import hmg_model as h
from scipy.optimize import brentq
for lab, g in [("vacuum clock, Planck 2018 vacuum", 1.0556e-10), ("vacuum clock, HMG Omega_L 0.6937", 1.0625e-10)]:
    s = brentq(lambda s: h.g_s_dl(s) - g, 0.05, 0.999)
    print(f"{lab}: g = {g:.4e} -> s = {s:.3f}   (fitted: late-type 0.329 +- 0.023, blue 0.314 +- 0.023)")
