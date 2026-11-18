# Exact (spherical-Bessel) vs Limber C_ell for kk, gg, kg, linear P with scale-independent growth (mnu=0),
# same kernels and CoLoRe windows (v2 kp) as nbpipe; kappa kernel with CoLoRe's chi_* (no radiation).
from nbpipe import *
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
from common import CoLoReCosmo
cc=CoLoReCosmo(cd,D+'ensemble_1/pk_linear.dat'); cs=float(cc.chi(cd['z_out']))
zg=np.linspace(1e-4,z_max,4000); chig=R.comoving_radial_distance(zg); Hg=R.hubble_parameter(zg)
Dg=np.sqrt(PK.P(zg,0.05)/PK.P(0,0.05))
kk_=np.logspace(-4,0.3,2500); P0=PK.P(0,kk_)
def Wk_(chi_,zz): return 1.5*(H0/C_KMS)**2*om*(cs-chi_)/cs*chi_*(1+zz)
def Wg_(cat,zz):
    nz=load('ensemble_1',cat)['nz']; zb,b=np.loadtxt(D+f'colore_setup/inputs/bias/bz_{cat}.txt',unpack=True)
    return np.nan_to_num(PchipInterpolator(zb,b,extrapolate=False)(zz))*np.nan_to_num(nz(zz))*R.hubble_parameter(zz)/C_KMS
def exact(WA,WB,kp,L):
    # C = 2/pi int dk k^2 P0 S(k) IA IB,  I = int dchi W(chi) D(chi) j_l(k chi)
    f=np.exp(-(kk_/kp)**2)*P0*kk_**2; out=[]
    for l in L:
        J=spherical_jn(l,np.outer(kk_,chig))
        IA=simpson(J*(WA*Dg)[None,:],x=chig,axis=1); IB=IA if WB is WA else simpson(J*(WB*Dg)[None,:],x=chig,axis=1)
        out.append(2/np.pi*simpson(f*IA*IB,x=kk_))
    return np.array(out)
def limber(WA,WB,kp,L):
    out=[]
    for l in L:
        kl=(l+.5)/chig; P=np.interp(kl,kk_,P0)*np.exp(-(kl/kp)**2)*Dg**2
        out.append(simpson(WA*WB/chig**2*P,x=chig))
    return np.array(out)
if __name__=='__main__':
    L=[30,40,50,70,100,150,200,300]
    Wkg=Wk_(chig,zg); res={'L':L}
    r=exact(Wkg,Wkg,KP['kk'],L)/limber(Wkg,Wkg,KP['kk'],L); res['kk']=r.tolist()
    print('kk exact/Limber:',' '.join(f'{l}:{x:.4f}' for l,x in zip(L,r)),flush=True)
    for nm,cat in SAMPLES.items():
        Wgg=Wg_(cat,zg)
        rg=exact(Wgg,Wgg,KP['gg'],L)/limber(Wgg,Wgg,KP['gg'],L); rk=exact(Wkg,Wgg,KP['kg'],L)/limber(Wkg,Wgg,KP['kg'],L)
        res[nm]={'gg':rg.tolist(),'kg':rk.tolist()}
        print(f'{nm} gg:',' '.join(f'{l}:{x:.4f}' for l,x in zip(L,rg)),'| kg:',' '.join(f'{l}:{x:.4f}' for l,x in zip(L,rk)),flush=True)
    json.dump(res,open('t46_limber.json','w'))
    