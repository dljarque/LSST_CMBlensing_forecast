from pipeline import *
cd=json.load(open(D+'test_box_7/cosmo_params.json'))
for box in ['test_box_7','ensemble_1','test_box_10']:
    B=load_box(box,'desi_qso_bin0'); cdb=B['cd']; ccb=CoLoReCosmo(cdb,D+f'{box}/pk_linear.dat'); h=cdb['h']
    zz=np.linspace(1e-4,3.8,400); _,_,W,_=kernels(ccb,zz,ccb.chi(cdb['z_out']))
    e=np.arange(2,700); t=limber(ccb,e,zz,W,W,smooth_kp=h/2)
    d=B['kk']['cl'][2:700]; out=[]
    for a,b in [(2,20),(20,50),(50,100),(100,200),(200,300),(300,500),(500,700)]:
        m=(e>=a)&(e<b); w=(2*e[m]+1); r=np.sum(w*d[m])/np.sum(w*t[m]); err=np.sqrt(2/np.sum(w))  # Gaussian frac error
        out.append(f'[{a},{b}) {r:.4f}±{err:.4f}')
    print(box,'(fid s8=%.4f)'%cdb['sigma_8'],' data/theory(CoLoRe, no pixwin, R=2):',' '.join(out))
