from pipeline import *
from scipy.optimize import minimize
import sys
boxes=['test_box_11','ensemble_1','ensemble_2','ensemble_3']
LMIN,LMAX,NB=30,1000,40
ells=np.arange(LMIN,LMAX+1)
Bt=make_hybrid_bin_matrix(ells,nbins=NB,ell_min=LMIN,ell_max=LMAX)
edges=np.logspace(np.log10(LMIN),np.log10(LMAX),NB+1)   # same edges as make_hybrid_bin_matrix (lmin>=ell_switch)
pw=pw1024[LMIN:LMAX+1]
cd=json.load(open(D+'test_box_11/cosmo_params.json')); cc=CoLoReCosmo(cd,D+'test_box_11/pk_linear.dat'); h=cd['h']; s0=cd['sigma_8']
z=np.linspace(1e-4,3.8,600); chi=cc.chi(z); H=cc.H(z); cs=cc.chi(cd['z_out'])
k=(ells[None,:]+.5)/chi[:,None]; P=cc.P(np.broadcast_to(z[:,None],k.shape).ravel(),k.ravel()).reshape(k.shape)
sw=simpson(np.eye(len(z)),x=z,axis=1) if False else np.array([simpson(np.eye(len(z))[i],x=z) for i in range(len(z))])
_,_,Wk,_=kernels(cc,z,cs)
def lim(WA,WB,kp): return sw@((C/H)[:,None]*WA[:,None]*WB[:,None]/chi[:,None]**2*P*np.exp(-(k/kp)**2))
KP_KK=0.29574; kk_th=lim(Wk,Wk,KP_KK)
def data(box,cat):
    out={}
    for xx,f in [('gg',f'cl_gg_{cat}'),('kg',f'cl_kg_{cat}'),('kk','cl_kk')]:
        dd=np.load(D+f'{box}/cls/{f}.npz'); out[xx]=dd['cls'][LMIN:LMAX+1]; out['sn']=float(dd['shot_noise']) if 'shot_noise' in dd.files else out.get('sn',0)
    return out
def gal_theory(box,cat):
    B=load_box(box,cat)
    _,_,_,Wg=kernels(cc,z,cs,B['nz'],B['bz'])
    return lim(Wg,Wg,h/2)*pw**2, lim(Wk,Wg,h/2)*pw
def cov_full(gg,kg,kk,sn,N=1):
    L=ells; G=gg+sn
    vgg=2*G**2/(2*L+1); vkg=(kg**2+kk*G)/(2*L+1); xc=2*kg*G/(2*L+1)
    return [Bt@np.diag(v)@Bt.T/N for v in (vgg,vkg,xc)]
def fit_joint(dgg,dkg,tg,tk,sn,mask,N=1):
    one=Bt@np.ones(len(ells)); d=np.concatenate([(Bt@dgg)[mask],(Bt@dkg)[mask]])
    m=lambda p: np.concatenate([((Bt@tg)*p[0]**2*(p[1]/s0)**2+p[2]*sn*one)[mask], ((Bt@tk)*p[0]*(p[1]/s0)**2)[mask]])
    A0,N0=1.,0.
    for it in range(2):
        Cgg,Ckg,X=cov_full(tg*A0**2+N0,tk*A0,kk_th,sn,N); Cf=np.block([[Cgg[np.ix_(mask,mask)],X[np.ix_(mask,mask)]],[X[np.ix_(mask,mask)].T,Ckg[np.ix_(mask,mask)]]]); Ci=np.linalg.inv(Cf)
        f=lambda p:(d-m(p))@Ci@(d-m(p))
        r=minimize(f,[1,s0,0],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000)); r=minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000))
        A0,N0=r.x[0],r.x[2]*sn
    J=np.array([(m(r.x+e)-m(r.x-e))/(2*e.max()) for e in np.diag([1e-5,1e-5,1e-3])]).T; err=np.sqrt(np.diag(np.linalg.inv(J.T@Ci@J)))
    return r.x,err,r.fun,len(d)-3
def fit_kk(dkk,mask,N=1):
    t=(Bt@kk_th)[mask]; d=(Bt@dkk)[mask]; C=(Bt@np.diag(2*kk_th**2/(2*ells+1))@Bt.T/N)[np.ix_(mask,mask)]; Ci=np.linalg.inv(C)
    a=(t@Ci@d)/(t@Ci@t); ea=1/np.sqrt(t@Ci@t); return np.sqrt(a),0.5*ea/np.sqrt(a),(d-a*t)@Ci@(d-a*t),mask.sum()-1
if __name__=='__main__':
    print('== kappa-kappa per box (kp_kk=0.2957, no pixwin) ==')
    allkk=[]
    for b in boxes:
        dk=data(b,'desi_qso_bin0')['kk']; allkk.append(dk)
        row=[]
        for a_,b_ in [(30,100),(100,200),(200,300),(300,500),(500,700),(700,1000)]:
            m=(ells>=a_)&(ells<b_); w=2*ells[m]+1; row.append(f'{np.sum(w*dk[m])/np.sum(w*kk_th[m]):.4f}')
        s,e,c2,nd=fit_kk(dk,edges[1:]<=1000.01)
        s5,e5,_,_=fit_kk(dk,edges[1:]<=500)
        print(f'{b:12s} data/theory by ell range [30-100,100-200,200-300,300-500,500-700,700-1000]: {" ".join(row)} | s8/fid lmax500={s5:.4f}±{e5:.4f} lmax1000={s:.4f}±{e:.4f} χ2={c2:.0f}/{nd}')
    dm=np.mean(allkk,axis=0); s,e,c2,nd=fit_kk(dm,edges[1:]<=1000.01,N=4); s5,e5,_,_=fit_kk(dm,edges[1:]<=500,N=4)
    row=[]
    for a_,b_ in [(30,100),(100,200),(200,300),(300,500),(500,700),(700,1000)]:
        m=(ells>=a_)&(ells<b_); w=2*ells[m]+1; row.append(f'{np.sum(w*dm[m])/np.sum(w*kk_th[m]):.4f}')
    print(f'{"mean of 4":12s} data/theory: {" ".join(row)} | s8/fid lmax500={s5:.4f}±{e5:.4f} lmax1000={s:.4f}±{e:.4f} χ2={c2:.0f}/{nd}')
    sys.stdout.flush()
    for cat,rmins in [('desi_qso_bin0',[60,40]),('lsst_gal_bin0',[60,40])]:
        tg,tk=gal_theory('test_box_11',cat)
        chi_h=cc.chi(2.0 if 'qso' in cat else 1.0)*h
        for rmin in rmins:
            lmax=int(2*np.pi/rmin*chi_h-0.5); mask=edges[1:]<=lmax
            res=[]; Ds=[]
            for b in boxes:
                dd=data(b,cat); Ds.append(dd)
                tgb,tkb=gal_theory(b,cat)
                x,err,c2,nd=fit_joint(dd['gg'],dd['kg'],tgb,tkb,dd['sn'],mask); res.append((x,err))
                print(f'{cat} rmin={rmin} lmax={lmax} nbins={mask.sum()} {b:12s}: s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.2f} χ2={c2:.0f}/{nd}')
            s8s=np.array([r[0][1] for r in res]); es=np.array([r[1][1] for r in res])
            print(f'   -> mean of per-box s8 = {s8s.mean():.4f} ± {es.mean()/2:.4f} (Knox/√4) ; scatter across boxes = {s8s.std(ddof=1):.4f} vs mean Knox err {es.mean():.4f}; mean shift {(s8s.mean()-s0)/(es.mean()/2):+.1f}σ')
            gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
            x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=4)
            print(f'   -> fit to AVERAGED Cls (cov/4): s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.2f} χ2={c2:.0f}/{nd}')
            sys.stdout.flush()
