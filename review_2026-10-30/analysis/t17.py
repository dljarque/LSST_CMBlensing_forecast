from pipeline import *
import sys
def kkfit(box,lmax,R2,pix,lmin=20,nb=30):
    B=load_box(box,'desi_qso_bin0'); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat'); h=cd['h']
    e=np.arange(lmin,lmax+1); Bt=make_hybrid_bin_matrix(e,nbins=nb,ell_min=lmin,ell_max=lmax)
    Bd=make_hybrid_bin_matrix(B['kk']['ell'],nbins=nb,ell_min=lmin,ell_max=lmax); d=Bd@B['kk']['cl']
    zz=np.linspace(1e-4,3.8,400); _,_,W,_=kernels(cc,zz,cc.chi(cd['z_out']))
    t=limber(cc,e,zz,W,W,smooth_kp=h/np.sqrt(R2))*(pw1024[e]**2 if pix else 1)
    full=B['kk']['cl'].copy(); full[lmin:lmax+1]=t
    cov=Bd@np.diag(2*full**2/(2*B['kk']['ell']+1))@Bd.T; Ci=np.linalg.inv(cov); tb=Bt@t
    a=(tb@Ci@d)/(tb@Ci@tb); ea=1/np.sqrt(tb@Ci@tb); c2=(d-a*tb)@Ci@(d-a*tb)
    return np.sqrt(a),0.5*ea/np.sqrt(a),c2
for box in ['test_box_11','test_box_7','ensemble_1']:
    cc=CoLoReCosmo(load_box(box,'desi_qso_bin0')['cd'],D+f'{box}/pk_linear.dat'); L,dx=colore_grid(cc)
    for nm,R2,pix in [('R=2, pw^2 (current)',4.,True),('R^2+dx^2/4, no pixwin',4+dx**2/4,False)]:
        line=f'{box:11s} {nm:22s}'
        for lmax in [300,500,700,1000,1500]:
            s,es,c2=kkfit(box,lmax,R2,pix); line+=f' | {lmax}: {s:.4f}±{es:.4f} χ2={c2:.0f}'
        print(line); sys.stdout.flush()
