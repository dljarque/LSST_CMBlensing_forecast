from common import *
box,cat='test_box_7','desi_qso_bin0'
B=load_box(box,cat); cd=B['cd']
cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat')
cm1=CambCosmo(cd,kmax=1.0); cm20=CambCosmo(cd,kmax=20.)
for k in [0.5,0.9,1.0,1.2,2,5,20,100]:
    print('k=%6.1f  P_kmax1/P_kmax20 at z=0: %.4f   (P_kmax20 %.3e)'%(k,cm1.P(0.,k)/cm20.P(0.,k),cm20.P(0.,k)))
zg=np.linspace(1e-4,3.8,250); zf=np.linspace(1e-4,3.8,3000)
def cls(cos,z,ells,kp=0.327,kcut=None,chistar=None,Pfunc=None,Om=None):
    cs=cos.chi(cd['z_out']) if chistar is None else chistar
    _,_,Wk,Wg=kernels(cos,z,cs,B['nz'],B['bz'],Om=Om)
    return [limber(cos,ells,z,a,b,smooth_kp=kp,kcut=kcut,Pfunc=Pfunc) for a,b in [(Wg,Wg),(Wk,Wg),(Wk,Wk)]]
ells=np.array([20,50,100,200,300,561,1000,1500,2000,3000.])
ref=cls(cm1,zg,ells)
print('\nell:',ells)
def show(name,x):
    for i,xx in enumerate(['gg','kg','kk']): print('%-42s %s'%(name+' '+xx,' '.join('%+.4f'%v for v in x[i]/ref[i]-1)))
show('kmax=20 (no mask) vs notebook',cls(cm20,zg,ells))
show('kmax=1 + zero P(k>1/Mpc) vs notebook',cls(cm1,zg,ells,kcut=1.0))
show('nz=3000 z-grid vs 250',cls(cm1,zf,ells))
show('CoLoRe bg+growth+P(file) vs notebook',cls(cc,zf,ells))
show('kp=h/R=0.3366 vs 0.327 (camb)',cls(cm20,zg,ells,kp=cd['h']/2))
show('no smoothing vs kp=.327',cls(cm20,zg,ells,kp=None))
pw=hp.pixwin(1024,lmax=3071)
print('pixwin^2(nside1024):',' '.join('%.4f'%pw[int(l)]**2 for l in ells))
