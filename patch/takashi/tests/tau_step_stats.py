"""Optical depth per RT step in a production output, weighted by photon absorption (needs yt).
usage: tau_step_stats.py <run dir> <output number>
For the HI-ionising groups (5-7 of the 7-group DICE setup): the share of the photon absorption by tau_step =
sigma*c_red*nHI*dt_RT (= 0.27 x cell optical depth for rt_courant_factor=0.8), and the same split by the photon supply
a = sigma*c_red*Np*dt_RT (a >> 1: photon-rich cells where absorption balances recombination; a < 1: front side, where
skipping the convergence check over-ionises by about a factor 1+tau_step). TODO.md 3b / docs/merge_ks.md."""
import sys, yt, numpy as np
yt.set_log_level(40)
d,i=sys.argv[1],int(sys.argv[2])
ds=yt.load(f'{d}/output_{i:05d}'); ad=ds.all_data()
mH=1.6726e-24; X=0.76; c_red=0.01*2.99792458e10; cfl=0.8
sig={5:3.007e-18,6:5.687e-19,7:7.889e-20}       # HI cross sections of the ionizing groups (namelist)
sc=sorted(n for t,n in ds.field_list if t=='ramses' and 'scalar' in n)
rho=ad['ramses','Density'].to('g/cm**3').v if ('ramses','Density') in ds.field_list else ad['gas','density'].to('g/cm**3').v
dx=ad['index','dx'].to('cm').v; V=dx**3
xHI=ad['ramses',sc[0]].v; xHII=ad['ramses',sc[1]].v
nH=X*rho/mH; nHI=nH*xHI
print('ncell %d, levels dx[pc] min %.2f, M_gas %.3e Msun, M_HI %.3e, M_HII %.3e, M_H2 %.3e'%(rho.size,dx.min()/3.0857e18,(rho*V).sum()/1.989e33,(X*rho*V*xHI).sum()/1.989e33,(X*rho*V*xHII).sum()/1.989e33,(X*rho*V*np.clip(1-xHI-xHII,0,1)).sum()/1.989e33))
dt=cfl*dx/(3*c_red)                 # RT step of the cell's level
edges=[0,0.03,0.1,0.3,1,3,10,30,1e30]
aw={}; tot_abs=0; hist=np.zeros(len(edges)-1); a_num=0
for g,s in sig.items():
    names=[n for t,n in ds.field_list if t=='ramses-rt' and n.lower().startswith('photon_density') and n.rstrip().endswith(str(g))]
    Np=ad['ramses-rt',names[0]].v
    if g==5: print('group 5 field',names[0],'Np max %.3e median %.3e (yt units %s)'%(Np.max(),np.median(Np),ad['ramses-rt',names[0]].units))
    tau=s*c_red*nHI*dt              # optical depth of one RT step = 0.27 x cell optical depth
    a=s*c_red*Np*dt                 # ionisation strength per atom in one step (this group only)
    ab=Np*V*tau/(1+tau)             # photons absorbed in one step (takashi: old abundances)
    h,_=np.histogram(tau,bins=edges,weights=ab)
    # split by photon supply a (a>>1: photon-rich, absorption balances recombination, no ordering error)
    for lo,hi in ((0,0.1),(0.1,1),(1,10),(10,1e30)):
        m=(a>=lo)&(a<hi)
        aw.setdefault((lo,hi),[0,0,0]); aw[(lo,hi)][0]+=ab[m].sum(); aw[(lo,hi)][1]+=(ab[m]*tau[m]).sum(); aw[(lo,hi)][2]+=(ab[m]*(tau[m]>1)).sum()
    hist+=h; tot_abs+=ab.sum()
    print('group %d: absorbed/step %.3e, in-flight %.3e'%(g,ab.sum(),(Np*V).sum()))
print('fraction of the photon absorption by optical depth per RT step (all ionising groups):')
for lo,hi,h in zip(edges[:-1],edges[1:],hist): print('  tau_step %5.2f - %-7s: %6.2f %%'%(lo,('%.2f'%hi) if hi<1e29 else 'inf',100*h/tot_abs))
# cells: where are the thick, partially ionised cells
part=(xHII>0.05)&(xHII<0.95)
tau5=sig[5]*c_red*nHI*dt
print('partially ionised cells (0.05<xHII<0.95): %d (%.2f %% of cells), of which tau_step>1: %.1f %%, median nH %.3g'%(part.sum(),100*part.mean(),100*(tau5[part]>1).mean(),np.median(nH[part])))
print('absorption split by photon supply a = sigma*c*Np*dt (per group):')
for (lo,hi),v in sorted(aw.items()): print('  a %5.1f - %-6s: %6.2f %% of the absorption, mean tau_step %.3f, share with tau_step>1 %.2f %%'%(lo,('%.1f'%hi) if hi<1e29 else 'inf',100*v[0]/tot_abs,v[1]/max(v[0],1e-300),100*v[2]/tot_abs))
