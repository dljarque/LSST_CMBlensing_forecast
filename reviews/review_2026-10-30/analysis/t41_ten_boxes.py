# Ten-box (ensemble_1-10) QSO + LSST fits; test_box_11 excluded.
# t30.py takes its reference cosmology from test_box_11; re-run its setup with ensemble_1 (identical in all ten boxes).
import t30
t30._REF='ensemble_1'
exec("""
cd=json.load(open(D+_REF+'/cosmo_params.json')); cc=CoLoReCosmo(cd,D+_REF+'/pk_linear.dat'); h=cd['h']; s0=cd['sigma_8']
z=np.linspace(1e-4,3.8,600); chi=cc.chi(z); H=cc.H(z); cs=cc.chi(cd['z_out'])
k=(ells[None,:]+.5)/chi[:,None]; P=cc.P(np.broadcast_to(z[:,None],k.shape).ravel(),k.ravel()).reshape(k.shape)
_,_,Wk,_=kernels(cc,z,cs)
kk_th=lim(Wk,Wk,KP_KK)
""",t30.__dict__)
from t30 import *
import sys
boxes=[f'ensemble_{i}' for i in range(1,11)]; N=len(boxes)
KK={b:data(b,'desi_qso_bin0')['kk'] for b in boxes}
lb=Bt@ells
def T_of(excl):
    m=np.mean([KK[b] for b in boxes if b!=excl],axis=0)
    r=(Bt@m)/(Bt@kk_th); x=np.log(lb); w=1/((Bt**2)@(2*kk_th**2/(2*ells+1)))*(Bt@kk_th)**2
    c=np.polyfit(x,r,3,w=np.sqrt(w)); return np.polyval(c,np.log(ells))
out={'kk':{},'fits':{}}
# kk ratios by ell range
R=np.array([(Bt@KK[b])/(Bt@kk_th) for b in boxes])
for lo,hi in [(30,100),(100,200),(200,300),(300,500),(500,1000)]:
    m=(lb>=lo)&(lb<hi); v=R[:,m].mean(1); print(f'kk data/theory {lo}-{hi}: {v.mean():.4f} ± {v.std(ddof=1)/np.sqrt(N):.4f}  per-box: {" ".join(f"{q:.3f}" for q in v)}')
    out['kk'][f'{lo}-{hi}']=[v.mean(),v.std(ddof=1)/np.sqrt(N),list(v)]
for lmx in [300,1000]:
    mask=edges[1:]<=lmx; a,ea,_,_=fit_kk(np.mean([KK[b] for b in boxes],axis=0),mask,N=N); print(f'kk-only s8/s8fid lmax={lmx}: {a:.4f} ± {ea:.4f}')
T_all=T_of(None)
print('T(ell) all at 30,50,100,200,300,500,1000:',' '.join(f'{np.interp(l,ells,T_all):.4f}' for l in [30,50,100,200,300,500,1000]))
sys.stdout.flush()
for name,cat,z0,rmins in [('QSO','desi_qso_bin0',2.0,[60,50,40]),('LSST','lsst_gal_bin0',1.0,[60,40])]:
    tg,tk=gal_theory('ensemble_1',cat); lmaxs={}
    for rmin in rmins:
        lmax=int(2*np.pi/rmin*cc.chi(z0)*h-0.5); mask=edges[1:]<=lmax; Ds=[]; s8s=[]; es=[]; s8T=[]
        for b in boxes:
            dd=data(b,cat); Ds.append(dd); tgb,tkb=gal_theory(b,cat)
            x,err,c2,nd=fit_joint(dd['gg'],dd['kg'],tgb,tkb,dd['sn'],mask); s8s.append(x[1]); es.append(err[1])
            x1,_,_,_=fit_joint(dd['gg'],dd['kg'],tgb,tkb*np.sqrt(T_of(b)),dd['sn'],mask); s8T.append(x1[1])
            print(f'  {name} {b:11s} rmin={rmin}: s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd} | withT {x1[1]:.4f}')
        gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
        x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=N)
        xT,eT,c2T,_=fit_joint(gm,km,tg,tk*np.sqrt(T_all),snm,mask,N=N)
        sc=np.std(s8s,ddof=1)
        print(f'{name} rmin={rmin} lmax={lmax}: scatter {sc:.4f} (±{sc/np.sqrt(2*(N-1)):.4f}) vs Knox {np.mean(es):.4f}; mean per-box {np.mean(s8s):.4f}; AVERAGED({N}): s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ, {100*(x[1]/s0-1):+.2f}%) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd} | mean LOO-T {np.mean(s8T):.4f} ({(np.mean(s8T)-s0)/err[1]:+.1f}σ, {100*(np.mean(s8T)/s0-1):+.2f}%) | avg with T_all {xT[1]:.4f} ({(xT[1]-s0)/eT[1]:+.1f}σ)')
        out['fits'][f'{name}{rmin}']=dict(lmax=lmax,s8=s8s,e=es,s8T=s8T,avg=[x[1],err[1]],A=x[0],Asn=x[2],chi2=[c2,nd],avgT=[xT[1],eT[1]],scatter=sc)
        sys.stdout.flush()
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'t41_out.json'),'w'),default=float)
