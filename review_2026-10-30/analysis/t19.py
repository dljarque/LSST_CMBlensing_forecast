from pipeline import *
from scipy.optimize import minimize
for rmin in [100,60,40,30,20]:
    S=Setup('test_box_11','desi_qso_bin0',rmin=rmin); S.cm=CambCosmo(S.cd,kmax=20.); h=S.cd['h']
    tg,tk=S.theory('notebook',kp=h/2); tkk=kk_theory(S,'colore',kp=h/2)*S.pw**2
    Cgg,Ckg,X=cov_theory(S,tg,tk,tkk); C=np.block([[Cgg,X],[X.T,Ckg]]); Ci=np.linalg.inv(C)
    d=np.concatenate([S.d['gg'],S.d['kg']]); bg=S.Bt@tg; bk=S.Bt@tk; one=S.Bt@(S.pw**2)
    def m(p):
        A,s,N=p; g=(s/S.s8)**2; return np.concatenate([bg*A**2*g+N*1e-7*one,bk*A*g])
    f=lambda p:(d-m(p))@Ci@(d-m(p))
    r=minimize(f,[1,S.s8,0],method='Nelder-Mead',options=dict(xatol=1e-7,fatol=1e-9,maxiter=20000))
    eps=np.array([1e-4,1e-4,1e-3]); J=np.array([(m(r.x+e)-m(r.x-e))/(2*e.max()) for e in np.diag(eps)]).T
    err=np.sqrt(np.diag(np.linalg.inv(J.T@Ci@J)))
    print(f'rmin={rmin:3d} lmax={S.lmax:4d}  linear + free constant N: s8={r.x[1]:.4f}±{err[1]:.4f} ({(r.x[1]-S.s8)/err[1]:+.1f}σ) A={r.x[0]:.3f} N={r.x[2]*1e-7:.2e} (shot noise {S.B["gg"]["sn"]:.2e}) chi2={r.fun:.0f}/{len(d)-3}')
