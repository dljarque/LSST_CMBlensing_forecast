from lognormal import *
for box,cat in [('test_box_7','desi_qso_bin0'),('test_box_7','lsst_gal_bin0')]:
    cd=json.load(open(D+f'{box}/cosmo_params.json')); cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat')
    gf=GalaxyField(cc)
    print(box,cat,'L=%.1f Mpc/h dx=%.3f Mpc/h  reff=%.3f  sigma_G(z=0)=%.3f'%(*colore_grid(cc),np.sqrt(gf.reff2),np.sqrt(gf.sig2_0)))
    gc=cat.rsplit('_bin',1)[0]; I=D+'colore_setup/inputs/'
    zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
    zr,br=np.loadtxt(I+f'bias/bz_{gc}.txt',unpack=True)
    for z in ([1.8,1.9,2.0,2.1,2.2] if 'qso' in cat else [0.8,0.9,1.0,1.1,1.2]):
        b,t=np.interp(z,zi,bi),np.interp(z,zt,tt); c,s=gf.coeffs(z,b,t)
        print('  z=%.2f b_in=%.3f t=%+.3f sigma_G=%.3f  b_x(exact)=%.4f  b_real(file)=%.4f  ratio=%.4f   (c2/c0)^2/2=%.3f'%(z,b,t,s,c[1]/(c[0]*s),np.interp(z,zr,br),c[1]/(c[0]*s)/np.interp(z,zr,br),(c[2]/c[0])**2/2))
