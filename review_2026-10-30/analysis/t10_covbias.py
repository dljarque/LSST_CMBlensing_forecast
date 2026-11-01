import numpy as np
from utils.binning import make_hybrid_bin_matrix
rng=np.random.default_rng(1)
for lmax in [561,1000]:
    e=np.arange(20,lmax+1); Bm=make_hybrid_bin_matrix(e,nbins=30,ell_min=20,ell_max=lmax)
    C=1e-5*(e/100.)**-1.0+7.4e-6
    As_d=[];As_t=[]
    for i in range(3000):
        Ch=C*rng.chisquare(2*e+1)/(2*e+1); d=Bm@Ch; t=Bm@C
        vd=Bm**2@(2*Ch**2/(2*e+1)); vt=Bm**2@(2*C**2/(2*e+1))
        As_d.append(np.sum(t*d/vd)/np.sum(t*t/vd)); As_t.append(np.sum(t*d/vt)/np.sum(t*t/vt))
    As_d=np.array(As_d);As_t=np.array(As_t)
    print(f'lmax={lmax}: amplitude bias data-cov: {As_d.mean()-1:+.5f} (scatter {As_d.std():.4f}) -> {(As_d.mean()-1)/As_d.std():+.2f}σ ;  theory-cov: {As_t.mean()-1:+.5f}')
