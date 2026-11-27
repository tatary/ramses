#!/usr/bin/env python3
"""photon_budget.py L ixHII output_dir...: number of ionized H atoms and of photons in flight on a uniform level-L grid
(ixHII = 1-based index of xHII among the hydro variables; for tests/ifront_thick.nml: photon_budget.py 6 6 output_00005)."""
import sys, glob, struct, os, array, math, re
def rec(f, t):
    (nb,) = struct.unpack('i', f.read(4)); a = array.array({'i4':'i','f8':'d'}[t]); a.frombytes(f.read(nb)); f.read(4); return a
def levsum(out, pre, L, fn):
    tot=None
    for path in sorted(glob.glob(os.path.join(out,pre+'_*.out*'))):
        icpu=int(path.split('out')[-1])
        with open(path,'rb') as f:
            ncpu=rec(f,'i4')[0]; nvar=rec(f,'i4')[0]; ndim=rec(f,'i4')[0]; nl=rec(f,'i4')[0]; nb=rec(f,'i4')[0]; rec(f,'f8')
            for il in range(nl):
                for ib in range(nb+ncpu):
                    rec(f,'i4'); nc=rec(f,'i4')[0]
                    if nc==0: continue
                    for ind in range(8):
                        v=[rec(f,'f8') for iv in range(nvar)]
                        if ib+1==icpu and il+1==L:
                            r=fn(v); tot=r if tot is None else [a+b for a,b in zip(tot,r)]
    return tot
L=int(sys.argv[1]); ix=int(sys.argv[2])-1; mH=1.6726e-24
for out in sys.argv[3:]:
    i=int(out[-5:]); info=open(f'{out}/info_{i:05d}.txt').read(); irt=open(f'{out}/info_rt_{i:05d}.txt').read()
    g=lambda s,k: float(re.search(k+r'\s*=\s*(\S+)',s).group(1))
    ul,ud,ut,t,box=g(info,'unit_l'),g(info,'unit_d'),g(info,'unit_t'),g(info,'time'),g(info,'boxlen'); unp=g(irt,'unit_np')
    dV=(ul*box/2**L)**3
    nion=levsum(out,'hydro',L,lambda v:[math.fsum(d*x for d,x in zip(v[0],v[ix]))])[0]*ud/mH*dV
    ph=levsum(out,'rt',L,lambda v:[math.fsum(v[4*k]) for k in range(len(v)//4)]); ph=[p*unp*dV for p in ph]
    print(f'{out:28s} t={t*ut/3.156e13:.4f} Myr  N_ion {nion:.5e}  N_flight {sum(ph):.5e}  sum {nion+sum(ph):.5e}')
