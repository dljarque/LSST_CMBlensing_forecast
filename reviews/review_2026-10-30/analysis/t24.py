from pipeline import *
from t23 import build
from scipy.optimize import minimize
import sys
box,cat='test_box_11','desi_qso_bin0'
for rmin in [40,20]:
    S=Setup(box,cat,rmin=rmin); h=S.cd['h']; s0=S.s8; T=build(S,False)
    cl=lambda kp,X: T['sw']@(T[X]*np.exp(-(T['k']/kp)**2))
    g0=cl(h/2,'gg'); k0=cl(h/2,'kg'); kk0=cl(h/2,'kk')
    Cgg,Ckg,X=cov_theory(S,g0*S.pw**2,k0*S.pw**2,kk0); Ci=np.linalg.inv(np.block([[Cgg,X],[X.T,Ckg]])); d=np.concatenate([S.d['gg'],S.d['kg']])
    one=S.Bt@S.pw**2
    def model(p,mode):
        A,s,N=p[0],p[1],p[2]*1e-7; f=(s/s0)**2
        kgg,kkg=(p[3],p[4]) if mode=='free' else (h/2,h/2)
        return np.concatenate([S.Bt@(cl(kgg,'gg')*A**2*f*S.pw**2)+N*one, S.Bt@(cl(kkg,'kg')*A*f*S.pw)])
    for mode in ['fixed','free']:
        x0=[1,s0,0]+([h/2,h/2] if mode=='free' else [])
        fn=lambda p:(d-model(p,mode))@Ci@(d-model(p,mode)) if (mode=='fixed' or (0.05<p[3]<5 and 0.05<p[4]<5)) else 1e12
        r=minimize(fn,x0,method='Nelder-Mead',options=dict(xatol=1e-7,fatol=1e-9,maxiter=60000,maxfev=60000))
        r=minimize(fn,r.x,method='Nelder-Mead',options=dict(xatol=1e-7,fatol=1e-9,maxiter=60000,maxfev=60000))
        st=np.array([1e-4,1e-4,1e-2,1e-4,1e-4])[:len(r.x)]
        J=np.array([(model(r.x+e,mode)-model(r.x-e,mode))/(2*e.max()) for e in np.diag(st)]).T
        Cv=np.linalg.inv(J.T@Ci@J); err=np.sqrt(np.diag(Cv)); cr=Cv/np.outer(err,err)
        ex=f' kp_gg={r.x[3]:.3f}±{err[3]:.3f} kp_kg={r.x[4]:.3f}±{err[4]:.3f} corr(N,kp_gg)={cr[2,3]:+.2f} corr(s8,kp_kg)={cr[1,4]:+.2f}' if mode=='free' else ''
        print(f'rmin={rmin} linear+N, kp {mode}: s8={r.x[1]:.4f}±{err[1]:.4f} ({(r.x[1]-s0)/err[1]:+.1f}σ) A={r.x[0]:.3f} N={r.x[2]*1e-7:.2e}±{err[2]*1e-7:.1e} χ2={r.fun:.1f}{ex}'); sys.stdout.flush()
# kappa-kappa: sigma8 & kp free vs fixed
B=load_box(box,cat); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat'); h=cd['h']
zz=np.linspace(1e-4,3.8,400); _,_,W,_=kernels(cc,zz,cc.chi(cd['z_out']))
L,dx=colore_grid(cc); kp_code=h/np.sqrt(4+0.2261*dx**2)
for lmax in [500,1000,1500]:
    e=np.arange(20,lmax+1); Bt=make_hybrid_bin_matrix(e,nbins=30,ell_min=20,ell_max=lmax); Bd=make_hybrid_bin_matrix(B['kk']['ell'],nbins=30,ell_min=20,ell_max=lmax); d=Bd@B['kk']['cl']
    chi=cc.chi(zz); k=(e[None,:]+.5)/chi[:,None]; P=cc.P(np.broadcast_to(zz[:,None],k.shape).ravel(),k.ravel()).reshape(k.shape)
    base=(C/cc.H(zz))[:,None]*W[:,None]**2/chi[:,None]**2*P; sw=np.array([simpson(np.eye(len(zz))[i],x=zz) for i in range(len(zz))])
    th=lambda kp: Bt@(sw@(base*np.exp(-(k/kp)**2)))
    full=B['kk']['cl'].copy(); full[20:lmax+1]=sw@(base*np.exp(-(k/kp_code)**2)); Ci=np.linalg.inv(Bd@np.diag(2*full**2/(2*B['kk']['ell']+1))@Bd.T)
    m=lambda p:(p[0]/cd['sigma_8'])**2*th(p[1])
    f=lambda p:(d-m(p))@Ci@(d-m(p))
    r=minimize(f,[cd['sigma_8'],kp_code],method='Nelder-Mead',options=dict(xatol=1e-8,fatol=1e-10))
    J=np.array([(m(r.x+ee)-m(r.x-ee))/(2*ee.max()) for ee in np.diag([1e-5,1e-5])]).T; Cv=np.linalg.inv(J.T@Ci@J); er=np.sqrt(np.diag(Cv))
    a=(th(kp_code)@Ci@d)/(th(kp_code)@Ci@th(kp_code)); ea=1/np.sqrt(th(kp_code)@Ci@th(kp_code))
    print(f'kk lmax={lmax}: kp fixed (code {kp_code:.4f}): s8/fid={np.sqrt(a):.4f}±{0.5*ea/np.sqrt(a):.4f} | kp free: s8/fid={r.x[0]/cd["sigma_8"]:.4f}±{er[0]/cd["sigma_8"]:.4f}, kp={r.x[1]:.4f}±{er[1]:.4f}, corr={Cv[0,1]/er[0]/er[1]:+.2f}')
