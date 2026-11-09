from pipeline import *
from scipy.special import spherical_jn
for box,cat in [('test_box_7','desi_qso_bin0'),('test_box_7','lsst_gal_bin0')]:
    B=load_box(box,cat); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat'); h=cd['h']
    z=np.linspace(B["zn"].min(),B["zn"].max(),500); chi=cc.chi(z)
    W=np.nan_to_num(B['nz'](z))*np.nan_to_num(B['bz'](z))  # per dz
    Dz=cc.Dz(z); k=np.logspace(-3.3,-0.3,1500)
    P0=lambda kk: cc.P(0.,kk)*np.exp(-(kk/(h/2))**2)
    # kappa kernel for kg
    cs=cc.chi(cd['z_out']); zk=np.linspace(1e-4,3.8,1200); chik=cc.chi(zk)
    Wk=1.5*(100*h/C)**2*cd['omega_M']*chik*(cs-chik)/cs*(1+zk)*C/cc.H(zk)   # per dz
    print(box,cat)
    for l in [10,20,30,50,80,120]:
        Ig=simpson((W*Dz)[None,:]*spherical_jn(l,k[:,None]*chi[None,:]),x=z,axis=1)
        Ik=simpson((Wk*cc.Dz(zk))[None,:]*spherical_jn(l,k[:,None]*chik[None,:]),x=zk,axis=1)
        Cex=2/np.pi*simpson(k**2*P0(k)*Ig**2,x=k); Ckgex=2/np.pi*simpson(k**2*P0(k)*Ig*Ik,x=k)
        kl=(l+0.5)/chi; Cl=simpson((cc.H(z)/C)*W**2*Dz**2/chi**2*P0(kl),x=z)
        _,_,Wk2,Wg2=kernels(cc,zk,cs,B['nz'],B['bz']); Ckl=limber(cc,np.array([l]),zk,Wk2,np.nan_to_num(Wg2),smooth_kp=h/2)[0]
        print(f'   ell={l:4d}  gg: Limber/exact-1 = {Cl/Cex-1:+.4f}    kg: Limber/exact-1 = {Ckl/Ckgex-1:+.4f}')
