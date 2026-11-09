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
from scipy.optimize import minimize_scalar
boxes=[f'ensemble_{i}' for i in range(1,11)]; N=len(boxes)
KK={b:data(b,'desi_qso_bin0')['kk'] for b in boxes}
lb=Bt@ells
def T_of(excl):
    m=np.mean([KK[b] for b in boxes if b!=excl],axis=0)
    r=(Bt@m)/(Bt@kk_th); x=np.log(lb); w=1/((Bt**2)@(2*kk_th**2/(2*ells+1)))*(Bt@kk_th)**2
    c=np.polyfit(x,r,3,w=np.sqrt(w)); return np.polyval(c,np.log(ells))
T_all=T_of(None)
km=np.mean([KK[b] for b in boxes],axis=0)
R=np.array([(Bt@KK[b])/(Bt@kk_th) for b in boxes])
th=Bt@kk_th; thT=Bt@(kk_th*T_all)
knox=np.sqrt(np.diag(Bt@np.diag(2*kk_th**2/(2*ells+1))@Bt.T))/th/np.sqrt(N)
out={'lb':lb.tolist(),'r':((Bt@km)/th).tolist(),'rT':((Bt@km)/thT).tolist(),'sem':(R.std(0,ddof=1)/np.sqrt(N)).tolist(),'knox':knox.tolist(),'T':(thT/th).tolist()}
# kk chi2: model a*theory (a = (s8/s8fid)^2), cov Knox/N
def chi2kk(dkk,t,mask,a):
    C=(Bt@np.diag(2*kk_th**2/(2*ells+1))@Bt.T/N)[np.ix_(mask,mask)]; d=(Bt@dkk)[mask]-a*(Bt@t)[mask]; return d@np.linalg.solve(C,d)
for lmx in [300,1000]:
    mask=edges[1:]<=lmx
    for lab,t in [('noT',kk_th),('T_all',kk_th*T_all)]:
        r=minimize_scalar(lambda a: chi2kk(km,t,mask,a),bounds=(0.8,1.2),method='bounded')
        c0=chi2kk(km,t,mask,1.0)
        print(f'kk lmax={lmx} {lab}: best s8/fid={np.sqrt(r.x):.4f} chi2_min={r.fun:.1f}/{mask.sum()-1}  chi2(fid)={c0:.1f}/{mask.sum()}  dchi2={c0-r.fun:.1f}')
        out[f'kk{lmx}{lab}']=[float(np.sqrt(r.x)),float(r.fun),int(mask.sum()-1),float(c0),int(mask.sum())]
# joint fits: chi2 at best fit vs sigma8 fixed to fiducial (A, A_sn refit)
def fit_fixed(dgg,dkg,tg,tk,sn,mask,N):
    one=Bt@np.ones(len(ells)); d=np.concatenate([(Bt@dgg)[mask],(Bt@dkg)[mask]])
    m=lambda p: np.concatenate([((Bt@tg)*p[0]**2+p[1]*sn*one)[mask], ((Bt@tk)*p[0])[mask]])
    A0,N0=1.,0.
    for it in range(2):
        Cgg,Ckg,X=cov_full(tg*A0**2+N0,tk*A0,kk_th,sn,N); Cf=np.block([[Cgg[np.ix_(mask,mask)],X[np.ix_(mask,mask)]],[X[np.ix_(mask,mask)].T,Ckg[np.ix_(mask,mask)]]]); Ci=np.linalg.inv(Cf)
        f=lambda p:(d-m(p))@Ci@(d-m(p))
        r=minimize(f,[1,0],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000)); r=minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000))
        A0,N0=r.x[0],r.x[1]*sn
    return r.x,r.fun,len(d)-2
for name,cat,z0,rmins in [('QSO','desi_qso_bin0',2.0,[60,40]),('LSST','lsst_gal_bin0',1.0,[60,40])]:
    tg,tk=gal_theory('ensemble_1',cat)
    Ds=[data(b,cat) for b in boxes]
    gm=np.mean([d['gg'] for d in Ds],axis=0); kgm=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
    for rmin in rmins:
        lmax=int(2*np.pi/rmin*cc.chi(z0)*h-0.5); mask=edges[1:]<=lmax
        for lab,tkk in [('noT',tk),('T_all',tk*np.sqrt(T_all))]:
            x,err,c2,nd=fit_joint(gm,kgm,tg,tkk,snm,mask,N=N)
            xf,c2f,ndf=fit_fixed(gm,kgm,tg,tkk,snm,mask,N)
            print(f'{name} rmin={rmin} {lab}: best s8={x[1]:.4f} chi2_min={c2:.1f}/{nd} | s8=fid: chi2={c2f:.1f}/{ndf} (A={xf[0]:.3f}, A_sn={xf[1]:.3f}) dchi2={c2f-c2:.1f} sqrt={np.sqrt(max(c2f-c2,0)):.1f}')
            out[f'{name}{rmin}{lab}']=[float(x[1]),float(c2),int(nd),float(c2f),int(ndf)]
        sys.stdout.flush()
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'t42_out.json'),'w'))
