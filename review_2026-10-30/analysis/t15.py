from pipeline import *
box='test_box_11'
B=load_box(box,'desi_qso_bin0'); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat'); h=cd['h']; cm=CambCosmo(cd,kmax=20.)
zz=np.linspace(1e-4,3.8,400)
e=np.arange(2,1600)
_,_,W,_=kernels(cc,zz,cc.chi(cd['z_out'])); t=limber(cc,e,zz,W,W,smooth_kp=h/2)
_,_,W2,_=kernels(cm,zz,cm.chi(cd['z_out'])); t2=limber(cm,e,zz,W2,W2,smooth_kp=h/2)
d=B['kk']['cl'][2:1600]; pw=pw1024[2:1600]
print('data/theory: CoLoRe no-win | CoLoRe pw^2 | CAMB pw^2')
for a,b in [(2,20),(20,50),(50,100),(100,200),(200,300),(300,500),(500,700),(700,1000),(1000,1500)]:
    m=(e>=a)&(e<b); w=2*e[m]+1; err=np.sqrt(2/np.sum(w))
    print(f'[{a},{b}) {np.sum(w*d[m])/np.sum(w*t[m]):.4f} | {np.sum(w*d[m])/np.sum(w*t[m]*pw[m]**2):.4f} | {np.sum(w*d[m])/np.sum(w*t2[m]*pw[m]**2):.4f}  ±{err:.4f}')
