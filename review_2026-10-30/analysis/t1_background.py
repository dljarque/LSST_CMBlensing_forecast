from common import *
cd=json.load(open(D+'test_box_7/cosmo_params.json'))
cc=CoLoReCosmo(cd,D+'test_box_7/pk_linear.dat'); cm=CambCosmo(cd,kmax=20)
print('CoLoRe-style sigma8 of file before renorm: %.6f  (target %.6f) -> P rescale %.5f'%(cc.sig8_file,cc.s8,cc.norm))
print('CAMB(mnu=0.06) sigma8_0 = %.6f'%cm.r.get_sigma8_0())
pars=cm.r.Params; print('CAMB Omega_m incl nu = %.5f ; omnuh2=%.3e; Om_nu=%.5f'%((pars.omch2+pars.ombh2+pars.omnuh2)/cd['h']**2,pars.omnuh2,pars.omnuh2/cd['h']**2))
for z in [0.5,1,2,3,3.8,1090]:
    print('z=%6.1f chi ratio CoLoRe/CAMB = %.5f   H ratio %.5f'%(z,cc.chi(z)/cm.chi(z),cc.H(z)/cm.H(z)) if z<1100 else '')
for z in [0.5,1,2,3]:
    for k in [0.01,0.05,0.2]:
        print('z=%.1f k=%.2f/Mpc  P_CoLoRe/P_CAMB = %.5f'%(z,k,cc.P(z,k)/cm.P(z,k)))
