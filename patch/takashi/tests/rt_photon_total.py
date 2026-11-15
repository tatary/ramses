#!/usr/bin/env python3
"""rt_photon_total.py output_dir...: total number of photons in flight per group.
Sums the photon density over the level-5 cells of the rt files (the coarse level must cover the box; non-leaf cells
hold the volume average of their children, so one fully covered level gives the box total). Edit L for other levelmin."""
import sys, glob, struct, os, array, math, re
def rec(f, t):
    (nb,) = struct.unpack('i', f.read(4)); a = array.array({'i4':'i','f8':'d'}[t]); a.frombytes(f.read(nb)); f.read(4); return a
for out in sys.argv[1:]:
    i=int(out[-5:]); info=open(f'{out}/info_{i:05d}.txt').read(); irt=open(f'{out}/info_rt_{i:05d}.txt').read()
    g=lambda s,k: float(re.search(k+r'\s*=\s*(\S+)',s).group(1))
    ul=g(info,'unit_l'); unp=g(irt,'unit_np'); L=5
    tot=None
    for fn in sorted(glob.glob(os.path.join(out,'rt_*.out*'))):
        icpu=int(fn.split('out')[-1])
        with open(fn,'rb') as f:
            ncpu=rec(f,'i4')[0]; nvar=rec(f,'i4')[0]; ndim=rec(f,'i4')[0]; nl=rec(f,'i4')[0]; nb=rec(f,'i4')[0]; rec(f,'f8')
            if tot is None: tot=[0.0]*nvar
            for il in range(nl):
                for ib in range(nb+ncpu):
                    rec(f,'i4'); nc=rec(f,'i4')[0]
                    if nc==0: continue
                    for ind in range(8):
                        for iv in range(nvar):
                            a=rec(f,'f8')
                            if ib+1==icpu and il+1==L: tot[iv]+=math.fsum(a)
    N=[tot[4*k]*unp*(ul/2**L)**3 for k in range(len(tot)//4)]
    print(out, 'a=%.5f'%g(info,'aexp'), 'photons in flight g1..g5:', ' '.join('%.4e'%x for x in N), '| ionizing (g2-4): %.4e'%sum(N[1:4]))
