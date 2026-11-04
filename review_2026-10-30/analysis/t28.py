from t26 import *
B=load_box(box,'desi_qso_bin0'); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat')
zz=np.linspace(1e-4,3.8,400); _,_,W,_=kernels(cc,zz,cc.chi(cd['z_out']))
for lo,hi in [(30,169),(30,359),(30,409),(30,561)]:
    e=np.arange(lo,hi+1); t=limber(cc,e,zz,W,W,smooth_kp=0.2957); d=B['kk']['cl'][lo:hi+1]; w=2*e+1
    r=np.sum(w*d)/np.sum(w*t); er=np.sqrt(2/np.sum(w))
    print(f'kappa-kappa data/theory ell {lo}-{hi}: {r:.4f} ± {er:.4f}  -> kappa amplitude {np.sqrt(r):.4f}')
