# The gauge couplings as the stack of medium layers below: alpha_s summed layer by layer (two loops, flavour thresholds)
# from alpha_s(M_Z) = 0.1180 to the bottom-quark and tau scales. Pass rule: within the measured errors.
# Usage: python stacked_running.py
import math
def run(a,mu0,mu1,nf,loops=2,n=4000):
    b0=(33-2*nf)/(12*math.pi); b1=(153-19*nf)/(24*math.pi**2); h=(math.log(mu1)-math.log(mu0))/n
    for _ in range(n): a+=-2*a*a*(b0+(b1*a if loops==2 else 0))*h
    return a
for loops in (1,2):
    amb=run(0.1180,91.19,4.18,5,loops); amt=run(amb,4.18,1.777,4,loops)
    print(f"{loops}-loop: alpha_s(m_b) = {amb:.3f} (measured 0.225 +- 0.010); alpha_s(m_tau) = {amt:.3f} (measured 0.325 +- 0.015)")
