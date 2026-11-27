"""ion_by_density.py <run dir> <output number>: ionised / neutral hydrogen mass and thermal energy in density bins,
photons in flight per group (needs yt)."""
import sys, yt, numpy as np
yt.set_log_level(40)
d,i=sys.argv[1],int(sys.argv[2])
ds=yt.load(f'{d}/output_{i:05d}'); ad=ds.all_data()
mH=1.6726e-24; X=0.76; Msun=1.989e33
sc=sorted(n for t,n in ds.field_list if t=='ramses' and 'scalar' in n)
rho=ad['ramses','Density'].to('g/cm**3').v; dx=ad['index','dx'].to('cm').v; V=dx**3
P=ad['ramses','Pressure'].to('erg/cm**3').v
xHI=ad['ramses',sc[0]].v; xHII=ad['ramses',sc[1]].v; xHeII=ad['ramses',sc[2]].v; xHeIII=ad['ramses',sc[3]].v
nH=X*rho/mH; mHy=X*rho*V/Msun
print('%s out%d t=%.6f ncell %d'%(d,i,float(ds.current_time.to('code_time').v) if hasattr(ds.current_time,'to') else 0,rho.size))
edges=[0,1e-3,1e-2,1e-1,1,10,100,1e30]
print('nH bin [cm^-3]      M_H [Msun]    M_HII [Msun]   M_HI [Msun]    M_H2 [Msun]   M_HeIII*      Eth [erg]')
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(nH>=lo)&(nH<hi)
    print('%8.0e - %-8s %.5e  %.5e  %.5e  %.5e  %.5e  %.5e'%(lo,('%.0e'%hi) if hi<1e29 else 'inf',mHy[m].sum(),(mHy*xHII)[m].sum(),(mHy*xHI)[m].sum(),(mHy*np.clip(1-xHI-xHII,0,1))[m].sum(),(mHy*xHeIII)[m].sum(),1.5*(P*V)[m].sum()))
print('total               %.5e  %.5e  %.5e  %.5e  %.5e  %.5e'%(mHy.sum(),(mHy*xHII).sum(),(mHy*xHI).sum(),(mHy*np.clip(1-xHI-xHII,0,1)).sum(),(mHy*xHeIII).sum(),1.5*(P*V).sum()))
ph=[]
for g in range(1,8):
    names=[n for t,n in ds.field_list if t=='ramses-rt' and n.lower().startswith('photon_density') and n.rstrip().endswith(str(g))]
    ph.append((ad['ramses-rt',names[0]].v*V).sum())
print('photons in flight g1..g7: '+' '.join('%.4e'%x for x in ph))
