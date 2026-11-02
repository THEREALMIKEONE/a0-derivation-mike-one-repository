# Usage: python btfr_sparc.py <data dir containing SPARC_Lelli2016c.mrt and sparc/*_rotmod.dat>
# BTFR test. Rule set before running: compare observed Vflat^4/(G Mb) with the value our curve
# predicts for the same galaxy (its own baryons, outer 3 points, Ydisk 0.5, Ybul 0.7, nothing fitted).
# Pass if mean log difference is within 2 sigma (bootstrap) AND within 0.1 dex (the M/L systematic).
import numpy as np, os, sys
G=6.674e-11; kpc=3.0857e19; Msun=1.989e30; g=1.0556e-10
law=lambda gb: gb/(1-np.exp(-np.sqrt(gb/g)))
D=sys.argv[1].rstrip('/')+'/'
obs=[];pre=[];gasd=[]
for L in open(D+'SPARC_Lelli2016c.mrt'):
    p=L.split()
    if len(p)<18: continue
    try: name=p[0]; Lum=float(p[7]); MHI=float(p[13]); V=float(p[15]); Q=int(p[17])
    except: continue
    if V<=0 or Q>2: continue
    f=D+'sparc/'+name+'_rotmod.dat'
    if not os.path.exists(f): continue
    a=np.loadtxt(f,comments='#'); a=a[-3:]
    R=a[:,0]*kpc; vg,vd,vb=a[:,3]*1e3,a[:,4]*1e3,a[:,5]*1e3
    gN=(vg*np.abs(vg)+0.5*vd**2+0.7*vb**2)/R
    if np.any(gN<=0): continue
    Vp=np.mean(np.sqrt(law(gN)*R))
    Mb=(0.5*Lum+1.33*MHI)*1e9*Msun
    obs.append(np.log10((V*1e3)**4/(G*Mb))); pre.append(np.log10(Vp**4/(G*Mb))); gasd.append(1.33*MHI>0.5*Lum)
obs,pre,gasd=map(np.array,(obs,pre,gasd)); d=obs-pre
rng=np.random.default_rng(2)
for lab,m in [("all",np.ones(len(d),bool)),("gas-dominated",gasd),("star-dominated",~gasd)]:
    bs=np.std([np.mean(rng.choice(d[m],m.sum())) for _ in range(4000)])
    print(f"{lab:15s} N={m.sum():3d} observed-predicted {np.mean(d[m]):+.3f} dex (+-{bs:.3f}) = {np.mean(d[m])/bs:+.1f} sigma, scatter {np.std(d[m]):.2f} dex")
