from pipeline import *
import sys
from scipy.optimize import minimize_scalar
cd=json.load(open(D+'test_box_7/cosmo_params.json')); cc=CoLoReCosmo(cd,D+'test_box_7/pk_linear.dat'); h=cd['h']
L,dx=colore_grid(cc)
# ---- fraction of kk Limber integrand from k > kX (smoothed, as in theory) ----
z=np.linspace(1e-4,3.8,4000); cs=cc.chi(cd['z_out']); _,_,Wk,_=kernels(cc,z,cs)
ells=np.array([50,100,200,300,500,700,1000,1500,2000])
chi=cc.chi(z); H=cc.H(z)
print('Fraction of C_ell^kk integrand from k > k_X  (CoLoRe linear P, smoothing kp=h/R=%.3f/Mpc)'%(h/2))
print('ell   '+' '.join('%7s'%('k>%.1f'%x) for x in [0.1,0.2,0.3,0.5,1.0])+'   z_peak  frac(z<0.3)')
for l in ells:
    k=(l+0.5)/chi; P=cc.P(z,k)*np.exp(-(k/(h/2))**2); I=C/H*Wk**2/chi**2*P; tot=simpson(I,x=z)
    fr=[simpson(np.where(k>x,I,0),x=z)/tot for x in [0.1,0.2,0.3,0.5,1.0]]
    print('%5d '%l+' '.join('%7.3f'%f for f in fr)+'   %.2f    %.3f'%(z[np.argmax(I)],simpson(np.where(z<0.3,I,0),x=z)/tot))
print('\nkappa from linear (Gaussian) potential in CoLoRe => no lognormal; dx=%.3f Mpc/h'%dx)
def kkfit(box, model, lmax, lmin=20, nb=30):
    B=load_box(box,'desi_qso_bin0'); cdb=B['cd']; ccb=CoLoReCosmo(cdb,D+f'{box}/pk_linear.dat')
    e=np.arange(lmin,lmax+1); Bt=make_hybrid_bin_matrix(e,nbins=nb,ell_min=lmin,ell_max=lmax)
    Bd=make_hybrid_bin_matrix(B['kk']['ell'],nbins=nb,ell_min=lmin,ell_max=lmax); d=Bd@B['kk']['cl']
    zz=np.linspace(1e-4,3.8,400)
    if model=='notebook':
        cm=CambCosmo(cdb,1.0); _,_,W,_=kernels(cm,zz,cm.chi(cdb['z_out'])); t=limber(cm,e,zz,W,W,smooth_kp=0.327)*pw1024[e]**2
    else:
        _,_,W,_=kernels(ccb,zz,ccb.chi(cdb['z_out']))
        R2={'nopix_R':4.,'pix_R':4.,'nopix_Rngp':4.+dx**2/12,'nopix_Rcol':4.+1.88*dx**2/6}[model]
        t=limber(ccb,e,zz,W,W,smooth_kp=h/np.sqrt(R2))*(pw1024[e]**2 if model.startswith('pix') else 1)
    tb=Bt@t; cov=Bd@np.diag(2*full_kk(B,t,lmin,lmax)**2/(2*B['kk']['ell']+1))@Bd.T; Ci=np.linalg.inv(cov)
    a=(tb@Ci@d)/(tb@Ci@tb); ea=1/np.sqrt(tb@Ci@tb); chi2=(d-a*tb)@Ci@(d-a*tb)
    return np.sqrt(a), 0.5*ea/np.sqrt(a), chi2, nb-1
def full_kk(B,t,lmin,lmax):
    o=B['kk']['cl'].copy(); o[lmin:lmax+1]=t; return o
for box in ['test_box_7','ensemble_1']:
    print(f'\n== kappa-kappa, {box}: sigma8/sigma8_fid from amplitude fit (theory covariance), 30 bins from ell=20')
    for model in ['notebook','pix_R','nopix_R','nopix_Rngp','nopix_Rcol']:
        line=f'{model:11s}'
        for lmax in [300,500,700,1000,1500,2000]:
            s,es,c2,nd=kkfit(box,model,lmax); line+=f' | lmax{lmax}: {s:.4f}±{es:.4f} ({(s-1)/es:+.1f}σ) χ2={c2:.0f}'
        print(line); sys.stdout.flush()
