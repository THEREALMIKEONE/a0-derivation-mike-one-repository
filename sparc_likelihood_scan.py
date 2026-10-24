# ============================================================
# SPARC likelihood scan: ΔBIC(fixed g† − free g†) on a grid of g† values, using the unchanged
# SPARC likelihood of the companion paper's test suite (cell p04). Output: dbic_curve.json.
# Usage: python sparc_likelihood_scan.py <GM_physical_plane_paper dir> <run dir with data/sparc> <out.json>
# ============================================================
import os, sys, math, io, contextlib, json
suite = sys.argv[1]; os.environ["GM_ROOT"] = sys.argv[2]; os.chdir(sys.argv[2])
G = {}
with contextlib.redirect_stdout(io.StringIO()):
    for name in ["p01_setup", "p02_geometry", "p03_map", "p04_galactic_gradient"]:
        f = os.path.join(suite, "cells", name + ".py"); exec(compile(open(f).read(), f, "exec"), G)
np = G["np"]; minimize = G["minimize"]; C = G["C"]; marg_nll = G["marg_nll"]
law = lambda gb, gd: gb / (1 - np.exp(-np.sqrt(gb / gd)))
nll = lambda lt, lg: marg_nll(np.log10(law(C["gbar_grid"], 10**lg)), lt, C)
n = len(C["ly"])
free = minimize(lambda x: nll(x[1], x[0]) if -11 < x[0] < -9 else 1e30, [-10.0, math.log(0.1)], method="Nelder-Mead")
grid = np.linspace(math.log10(3e-11), math.log10(5e-10), 70)
db = []
for lg in grid:
    b = minimize(lambda x: nll(x[0], lg), [math.log(0.1)], method="Nelder-Mead")
    db.append((2*b.fun + math.log(n)) - (2*free.fun + 2*math.log(n)))
json.dump({"grid": grid.tolist(), "dbic": db, "free": 10**free.x[0]}, open(sys.argv[3], "w"))
print("free", 10**free.x[0])
