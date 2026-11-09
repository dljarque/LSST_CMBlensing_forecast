from pipeline import *
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
for box,cat,gc in [('test_box_7','desi_qso_bin0','desi_qso'),('ensemble_1','desi_qso_bin0','desi_qso'),('test_box_7','lsst_gal_bin0','lsst_gal')]:
    f=np.load(D+f'{box}/cls/nz_{cat}.npz'); z,n=f['z_nz'],f['nz']; zm,dz,sz=float(f['z_mean']),float(f['dz']),float(f['sigma_z'])
    ff=np.load(D+f'{box}/cls/nz_{gc}.npz'); zf,nf=ff['z_nz'],ff['nz']
    zz=np.linspace(z.min()-0.1,z.max()+0.1,4000)
    S=norm.cdf((zm+dz/2-zz)/sz)-norm.cdf((zm-dz/2-zz)/sz)
    nfs=np.interp(zz,zf,gaussian_filter1d(nf,2)); nmod=nfs*S; nmod/=simpson(nmod,x=zz)
    nraw=np.nan_to_num(PchipInterpolator(z,n/simpson(n,x=z),extrapolate=False)(zz))
    i2r=simpson(nraw**2,x=zz); i2m=simpson(nmod**2,x=zz)
    print(f'{box} {cat}: N_bins_hist={len(z)}, int n^2 dz raw-hist/smooth-model = {i2r/i2m:.4f};  mean z raw {simpson(zz*nraw,x=zz):.4f} model {simpson(zz*nmod,x=zz):.4f}')
