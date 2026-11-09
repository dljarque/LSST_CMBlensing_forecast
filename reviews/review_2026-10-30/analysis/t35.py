from t30 import *
from lognormal import GalaxyField
boxes5=['test_box_11','ensemble_1','ensemble_2','ensemble_3','ensemble_4']
cat='desi_qso_bin0'; tg,tk=gal_theory('test_box_11',cat); chi_h=cc.chi(2.0)*h
for rmin in [60,50,40]:
    lmax=int(2*np.pi/rmin*chi_h-0.5); mask=edges[1:]<=lmax
    Ds=[data(b,cat) for b in boxes5]
    gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
    x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=5)
    print(f'QSO 5 good boxes AVERAGED rmin={rmin} lmax={lmax}: s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd}')
# Laura: configuration-space non-linear correction at large r
gf=GalaxyField(cc); I=D+'colore_setup/inputs/'
for gc,z0 in [('desi_qso',2.0),('lsst_gal',1.0)]:
    zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
    b,t=np.interp(z0,zi,bi),np.interp(z0,zt,tt); c,s=gf.coeffs(z0,b,t)
    rho=cc.Dz(z0)**2*gf.xi0/s**2
    from math import factorial
    lin=(c[1]/c[0])**2*rho; full=sum((c[n]/c[0])**2*rho**n/factorial(n) for n in range(1,gf.nherm))
    out=[]
    for r in [10,20,30,50,80]:
        i=np.argmin(abs(gf.r-r)); out.append(f'r={r}: {full[i]/lin[i]-1:+.3f}')
    print(gc,'z=',z0,'xi_gg(full)/xi_gg(linear b_x^2 xi_G) - 1 :',' '.join(out))
    # k-space ratio at k=0.05,0.1
    k,Pe=gf.P_excess(z0,b,t); kh=k/cc.h; Plin=(c[1]/(c[0]*s))**2*cc.Dz(z0)**2*np.exp(cc._lP(np.log(kh)))*np.exp(-kh**2*gf.Rxi2)
    print('   P_excess/P_lin at k=0.02,0.05,0.1 h/Mpc:',' '.join(f'{(Pe/cc.h**-3)[np.argmin(abs(kh-x))]/Plin[np.argmin(abs(kh-x))]:.3f}' for x in [0.02,0.05,0.1]))
