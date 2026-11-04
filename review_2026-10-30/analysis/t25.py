from pipeline import *
from t23 import build
from scipy.optimize import minimize
import sys
box,cat='test_box_7','lsst_gal_bin0'
for rmin in [80,60,40]:
  for lmin in [30]:
    S=Setup(box,cat,rmin=rmin,lmin=lmin,nbins=25); h=S.cd['h']; s0=S.s8; T=build(S,False)
    cl=lambda kp,X: T['sw']@(T[X]*np.exp(-(T['k']/kp)**2))
    g0=cl(h/2,'gg')*S.pw**2; k0=cl(h/2,'kg')*S.pw; kk0=cl(h/2,'kk')
    sn=S.B['gg']['sn']
    for it in range(2):
        A0,N0=(1.,0.) if it==0 else (r.x[0],r.x[2]*sn)
        Cgg,Ckg,X=cov_theory(S,g0*A0**2+N0,k0*A0,kk0); Ci=np.linalg.inv(np.block([[Cgg,X],[X.T,Ckg]])); d=np.concatenate([S.d['gg'],S.d['kg']])
        m=lambda p: np.concatenate([S.Bt@(g0*p[0]**2*(p[1]/s0)**2)+p[2]*sn*(S.Bt@np.ones_like(g0)), S.Bt@(k0*p[0]*(p[1]/s0)**2)])
        f=lambda p:(d-m(p))@Ci@(d-m(p))
        r=minimize(f,[1,s0,0],method='Nelder-Mead',options=dict(xatol=1e-8,fatol=1e-10,maxiter=40000))
    J=np.array([(m(r.x+e)-m(r.x-e))/(2*e.max()) for e in np.diag([1e-5,1e-5,1e-3])]).T; err=np.sqrt(np.diag(np.linalg.inv(J.T@Ci@J)))
    # sigma8-fixed single-probe amplitudes
    Cgi=np.linalg.inv(Cgg); Cki=np.linalg.inv(Ckg); tg=S.Bt@g0; tk=S.Bt@k0
    Akg=(tk@Cki@S.d['kg'])/(tk@Cki@tk)
    res=(S.d['gg']-(S.Bt@g0)*r.x[0]**2*(r.x[1]/s0)**2-r.x[2]*sn)/np.sqrt(np.diag(Cgg))
    print(f'{cat} rmin={rmin} lmax={S.lmax} lmin={lmin}: s8={r.x[1]:.4f}±{err[1]:.4f} ({(r.x[1]-s0)/err[1]:+.1f}σ) A={r.x[0]:.3f} A_sn={r.x[2]:.2f}±{err[2]:.2f} χ2={r.fun:.1f}/{len(d)-3}  A_kg(s8 fixed)={Akg:.3f}')
    print('   gg pulls:',' '.join('%+.1f'%x for x in res))
    sys.stdout.flush()
