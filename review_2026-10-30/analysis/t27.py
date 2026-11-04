from t26 import *
for cat in ['desi_bgs_bin0','desi_lrg_bin0','desi_elg_bin0','desi_qso_bin0','lsst_gal_bin0']:
    S=Setup(box,cat,rmin=40,lmin=30,nbins=25); h=S.cd['h']; s0=S.s8
    gc=cat.rsplit('_bin',1)[0]; I=D+'colore_setup/inputs/'; S.gf=GalaxyField(S.cc)
    zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
    zz=np.linspace(max(0.02,S.B['zn'].min()),S.B['zn'].max(),40)
    bx=np.array([S.gf.b_cross(x,np.interp(x,zi,bi),np.interp(x,zt,tt)) for x in zz]); S.B['bz']=interp1d(zz,bx,bounds_error=False,fill_value=np.nan)
    T=build(S,False); cl=lambda kp,X: T['sw']@(T[X]*np.exp(-(T['k']/kp)**2))
    g0=cl(h/2,'gg')*S.pw**2; k0=cl(h/2,'kg')*S.pw; kk0=cl(0.2957,'kk'); sn=S.B['gg']['sn']
    ex=S.hermite_excess()  # includes pw^2
    one=S.Bt@np.ones_like(g0)
    Cgg,Ckg,X=cov_theory(S,g0+ex,k0,kk0); Ci=np.linalg.inv(np.block([[Cgg,X],[X.T,Ckg]]))
    m=lambda p: np.concatenate([S.Bt@(g0*p[0]**2*(p[1]/s0)**2)+p[2]*sn*one, S.Bt@(k0*p[0]*(p[1]/s0)**2)])
    fake=np.concatenate([S.Bt@(g0+ex),S.Bt@k0])   # noiseless "truth" from the Hermite model
    f=lambda p:(fake-m(p))@Ci@(fake-m(p)); r=minimize(f,[1,s0,0],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=40000))
    J=np.array([(m(r.x+e)-m(r.x-e))/(2*e.max()) for e in np.diag([1e-5,1e-5,1e-3])]).T; err=np.sqrt(np.diag(np.linalg.inv(J.T@Ci@J)))
    # data fit with Hermite excess instead of constant
    d=np.concatenate([S.d['gg'],S.d['kg']])
    mh=lambda p: np.concatenate([(S.Bt@(g0*p[0]**2+ex))*(p[1]/s0)**2, S.Bt@(k0*p[0]*(p[1]/s0)**2)])
    fh=lambda p:(d-mh(p))@Ci@(d-mh(p)); rh=minimize(fh,[1,s0],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11))
    exr=(S.Bt@ex)/(S.Bt@g0)
    print(f'{cat:14s} excess/linear gg at bins 0,mid,last: {exr[0]:.3f} {exr[len(exr)//2]:.3f} {exr[-1]:.3f} | constant-model on noiseless Hermite truth: s8 bias {(r.x[1]-s0)/err[1]:+.1f}σ ({100*(r.x[1]/s0-1):+.2f}%) A_sn={r.x[2]:.2f} | data with Hermite: s8={rh.x[1]:.4f} ({(rh.x[1]-s0)/err[1]:+.1f}σ) χ2={rh.fun:.0f}')
    sys.stdout.flush()
