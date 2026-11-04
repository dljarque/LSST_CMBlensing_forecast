from pipeline import *
from t23 import build
from scipy.optimize import minimize
import sys
box='test_box_11'
def run(cat, rmin=40, lmin=30, ref='real'):
    S=Setup(box,cat,rmin=rmin,lmin=lmin,nbins=25); h=S.cd['h']; s0=S.s8
    if ref=='cross':
        gc=cat.rsplit('_bin',1)[0]; I=D+'colore_setup/inputs/'
        S.gf=GalaxyField(S.cc)
        zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
        zz=np.linspace(max(0.02,S.B['zn'].min()),S.B['zn'].max(),40)
        bx=np.array([S.gf.b_cross(x,np.interp(x,zi,bi),np.interp(x,zt,tt)) for x in zz]); S.B['bz']=interp1d(zz,bx,bounds_error=False,fill_value=np.nan)
    T=build(S,False)
    cl=lambda kp,X: T['sw']@(T[X]*np.exp(-(T['k']/kp)**2))
    g0=cl(h/2,'gg')*S.pw**2; k0=cl(h/2,'kg')*S.pw; kk0=cl(0.2957,'kk'); sn=S.B['gg']['sn']
    d=np.concatenate([S.d['gg'],S.d['kg']]); one=S.Bt@np.ones_like(g0)
    m=lambda p: np.concatenate([S.Bt@(g0*p[0]**2*(p[1]/s0)**2)+p[2]*sn*one, S.Bt@(k0*p[0]*(p[1]/s0)**2)])
    A0,N0=1.,0.
    for it in range(2):
        Cgg,Ckg,X=cov_theory(S,g0*A0**2+N0,k0*A0,kk0); Ci=np.linalg.inv(np.block([[Cgg,X],[X.T,Ckg]]))
        f=lambda p:(d-m(p))@Ci@(d-m(p))
        r=minimize(f,[1,s0,0],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000)); r=minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000))
        A0,N0=r.x[0],r.x[2]*sn
    J=np.array([(m(r.x+e)-m(r.x-e))/(2*e.max()) for e in np.diag([1e-5,1e-5,1e-3])]).T; err=np.sqrt(np.diag(np.linalg.inv(J.T@Ci@J)))
    Cgi=np.linalg.inv(Cgg); Cki=np.linalg.inv(Ckg); tk=S.Bt@k0
    Akg=(tk@Cki@S.d['kg'])/(tk@Cki@tk); eAkg=1/np.sqrt(tk@Cki@tk)
    # gg alone with s8 fixed: fit A,N
    tg=S.Bt@g0; M=np.vstack([tg,sn*one]).T; sol=np.linalg.solve(M.T@Cgi@M, M.T@Cgi@S.d['gg']); Agg=np.sqrt(sol[0])
    return S,r,err,Akg,eAkg,Agg,sol[1]
if __name__=='__main__':
  for cat in ['desi_bgs_bin0','desi_lrg_bin0','desi_elg_bin0','desi_qso_bin0','lsst_gal_bin0']:
    for ref in ['real','cross']:
        S,r,err,Akg,eAkg,Agg,Ngg=run(cat,ref=ref)
        print(f'{cat:14s} ref={ref:5s} lmax={S.lmax:4d} | joint: s8={r.x[1]:.4f}±{err[1]:.4f} ({(r.x[1]-S.s8)/err[1]:+.1f}σ) A={r.x[0]:.3f} A_sn={r.x[2]:+.2f}±{err[2]:.2f} χ2={r.fun:.0f}/{2*len(S.leff)-3} | s8 fixed: A_gg={Agg:.3f} (A_sn {Ngg:+.2f}) A_kg={Akg:.3f}±{eAkg:.3f} ratio={Agg/Akg:.3f}')
        sys.stdout.flush()
