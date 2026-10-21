#!/usr/bin/env python3
"""Order-independent statistics of a RAMSES hydro output (all cpu files).
Usage: hydro_stats.py output_dir
Prints, per variable: sum, sum of squares, min, max over all leaf+non-leaf cells
of every level in every cpu file (ghost grids are skipped: only ibound==myid)."""
import sys, glob, struct, os, array, math

def rec(f, dtype, n=None):
    (nb,) = struct.unpack('i', f.read(4))
    data = array.array({'i4':'i','f8':'d'}[dtype]); data.frombytes(f.read(nb))
    (nb2,) = struct.unpack('i', f.read(4))
    assert nb == nb2
    return data

out = sys.argv[1]
files = sorted(glob.glob(os.path.join(out, 'hydro_*.out*')))
acc = None
for fn in files:
    icpu = int(fn.split('out')[-1])
    with open(fn, 'rb') as f:
        ncpu = rec(f, 'i4')[0]; nvar = rec(f, 'i4')[0]; ndim = rec(f, 'i4')[0]
        nlevelmax = rec(f, 'i4')[0]; nboundary = rec(f, 'i4')[0]; gamma = rec(f, 'f8')[0]
        if acc is None:
            acc = {'sum': [0.0]*nvar, 'sq': [0.0]*nvar, 'min': [math.inf]*nvar,
                   'max': [-math.inf]*nvar, 'ncell': 0}
        for ilevel in range(nlevelmax):
            for ibound in range(nboundary + ncpu):
                il = rec(f, 'i4')[0]; ncache = rec(f, 'i4')[0]
                if ncache == 0:
                    continue
                for ind in range(2 ** ndim):
                    for ivar in range(nvar):
                        a = rec(f, 'f8')
                        if ibound + 1 == icpu:
                            acc['sum'][ivar] += math.fsum(a); acc['sq'][ivar] += math.fsum(x*x for x in a)
                            acc['min'][ivar] = min(acc['min'][ivar], min(a))
                            acc['max'][ivar] = max(acc['max'][ivar], max(a))
                            if ivar == 0:
                                acc['ncell'] += len(a)
print('ncell', acc['ncell'])
for ivar in range(len(acc['sum'])):
    print(f"var{ivar+1} sum={acc['sum'][ivar]:.16e} sq={acc['sq'][ivar]:.16e} min={acc['min'][ivar]:.16e} max={acc['max'][ivar]:.16e}")
