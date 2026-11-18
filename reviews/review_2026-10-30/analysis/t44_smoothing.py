# Smoothing-window comparison: notebook (R=2 for gg/kg, fitted kp_kk) vs v2 CoLoRe-derived windows
# exp(-(R^2 + c dx^2) k^2), c = 1/12 (gg), 13/84 (kg), 19/84 (kk). Ten-box averaged Cls, ensemble_1 cosmology.
import t30
t30._REF='ensemble_1'
exec("""
cd=json.load(open(D+_REF+'/cosmo_params.json')); cc=CoLoReCosmo(cd,D+_REF+'/pk_linear.dat'); h=cd['h']; s0=cd['sigma_8']
z=np.linspace(1e-4,3.8,600); chi=cc.chi(z); H=cc.H(z); cs=cc.chi(cd['z_out'])
k=(ells[None,:]+.5)/chi[:,None]; P=cc.P(np.broadcast_to(z[:,None],k.shape).ravel(),k.ravel()).reshape(k.shape)
_,_,Wk,_=kernels(cc,z,cs)
""",t30.__dict__)
from t30 import *
import t30 as T
R=2.0; dx=2*cc.chi(3.8)*h/4096   # CoLoRe grid: box = 2 chi(z_max), n_grid=4096
kp=lambda c: h/np.sqrt(R**2+c*dx**2)
V={'A_notebook':dict(gg=h/R,kg=h/R,kk=0.29574),
   'B_v2':dict(gg=kp(1/12),kg=kp(13/84),kk=kp(19/84)),
   'C_v2_kg_only':dict(gg=h/R,kg=kp(13/84),kk=kp(19/84))}
print(f'dx={dx:.4f} Mpc/h; kp: '+'; '.join(f'{n}: gg {v["gg"]:.4f} kg {v["kg"]:.4f} kk {v["kk"]:.4f}' for n,v in V.items()))
boxes=[f'ensemble_{i}' for i in range(1,11)]; N=10
KKm=np.mean([data(b,'desi_qso_bin0')['kk'] for b in boxes],axis=0)
for n,v in V.items():
    T.kk_th=lim(Wk,Wk,v['kk'])
    for lmx in [500,1000]:
        a,ea,c2,nd=fit_kk(KKm,edges[1:]<=lmx,N=N); print(f'kk-only {n} lmax={lmx}: s8/fid={a:.4f}±{ea:.4f} chi2={c2:.1f}/{nd}')
samples=[('BGS','desi_bgs_bin0',0.4,[80,100,40]),('LRG','desi_lrg_bin0',1.0,[100,150,40]),('ELG','desi_elg_bin0',1.2,[60,80,40]),
         ('QSO','desi_qso_bin0',2.0,[60,100,40]),('LSST','lsst_gal_bin0',1.0,[40,50])]
res={}
for name,cat,z0,rmins in samples:
    B=load_box('ensemble_1',cat); _,_,_,Wg=kernels(cc,z,cs,B['nz'],B['bz'])
    Ds=[data(b,cat) for b in boxes]; gm=np.mean([d['gg'] for d in Ds],0); km=np.mean([d['kg'] for d in Ds],0); snm=np.mean([d['sn'] for d in Ds])
    for rmin in rmins:
        lmax=int(2*np.pi/rmin*cc.chi(z0)*h-0.5); mask=edges[1:]<=lmax; line=f'{name:4s} rmin={rmin:3d} lmax={lmax:4d}:'
        for n,v in V.items():
            T.kk_th=lim(Wk,Wk,v['kk']); tg=lim(Wg,Wg,v['gg'])*pw**2; tk=lim(Wk,Wg,v['kg'])*pw
            x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=N); res[(name,rmin,n)]=(x,err,c2,nd)
            line+=f' | {n[:1]}: s8={x[1]:.4f}±{err[1]:.4f} ({100*(x[1]/s0-1):+.2f}%, {(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} χ2={c2:.0f}/{nd}'
        a=res[(name,rmin,'A_notebook')]; b=res[(name,rmin,'B_v2')]
        line+=f' || ΔB-A = {100*(b[0][1]/a[0][1]-1):+.2f}% = {(b[0][1]-a[0][1])/a[1][1]:+.2f}σ10 ({(b[0][1]-a[0][1])/(a[1][1]*np.sqrt(N)):+.2f}σ1box)'
        print(line,flush=True)
