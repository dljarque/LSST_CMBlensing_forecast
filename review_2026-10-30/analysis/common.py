import os
import numpy as np, json, camb, healpy as hp
from scipy.integrate import simpson, quad, cumulative_trapezoid, solve_ivp
from scipy.interpolate import PchipInterpolator, interp1d
import sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))
from utils.binning import make_hybrid_bin_matrix
C=299792.458
D=os.environ.get('CMBLYA_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../data/'))

def load_box(box, cat):
    cd=json.load(open(D+f'{box}/cosmo_params.json'))
    gc=cat.rsplit('_bin',1)[0]
    zb,b=np.loadtxt(D+f'colore_setup/inputs/bias/bz_{gc}.txt',unpack=True)
    nzf=np.load(D+f'{box}/cls/nz_{cat}.npz')
    out=dict(cd=cd, bz=PchipInterpolator(zb,b,extrapolate=False))
    zn,nn=nzf['z_nz'],nzf['nz']; out['nz']=PchipInterpolator(zn,nn/simpson(nn,x=zn),extrapolate=False)
    out['zn']=zn; out['nn']=nn; out['zmean_nom']=float(nzf['z_mean'])
    for xx,f in [('gg',f'cl_gg_{cat}'),('kg',f'cl_kg_{cat}'),('kk','cl_kk')]:
        d=np.load(D+f'{box}/cls/{f}.npz'); out[xx]=dict(ell=d['ells'],cl=d['cls'],err=d['err'],nside=int(d['nside']),
            sn=float(d['shot_noise']) if 'shot_noise' in d.files else 0.)
    return out

# ---------------- CoLoRe-consistent background (flat LCDM, no radiation, no nu) -------------
class CoLoReCosmo:
    def __init__(self, cd, pkfile):
        self.h=cd['h']; self.Om=cd['omega_M']; self.OL=cd['omega_L']; self.s8=cd['sigma_8']
        z=np.linspace(0,1200,1); 
        za=np.concatenate([np.linspace(0,10,20001),np.logspace(1.0001,np.log10(1200),2000)])
        E=np.sqrt(self.Om*(1+za)**3+self.OL)
        chi=cumulative_trapezoid(C/(100*self.h)/E,za,initial=0)   # Mpc
        self._chi=interp1d(za,chi,kind='cubic'); self._z=interp1d(chi,za,kind='cubic')
        # growth factor (Heath integral, exact for LCDM w=-1)
        def Dun(z):
            a=1/(1+z); f=lambda aa: 1/(aa*np.sqrt(self.Om/aa**3+self.OL))**3
            return 2.5*self.Om*np.sqrt(self.Om/a**3+self.OL)*quad(f,0,a,epsabs=0,epsrel=1e-10)[0]
        zg=np.linspace(0,10,1001); Dg=np.array([Dun(x) for x in zg]); self._D=interp1d(zg,Dg/Dg[0],kind='cubic')
        kh,Ph=np.loadtxt(pkfile,unpack=True)
        # CoLoRe: re-spline onto log-equispaced grid, normalise to sigma8 with top-hat R=8 Mpc/h over table range
        lk=np.log10(kh); lkg=np.linspace(lk[0],lk[-1],len(lk)); Pg=interp1d(lk,Ph,kind='cubic')(lkg)
        kk=10**lkg; x=kk*8; W=3*(np.sin(x)-x*np.cos(x))/x**3
        s2=simpson(kk**3*Pg*W**2/(2*np.pi**2),x=np.log(kk))
        self.sig8_file=np.sqrt(s2); self.norm=self.s8**2/s2
        self.kh=kk; self.Ph=Pg*self.norm
        self._lP=interp1d(np.log(kk),np.log(self.Ph),kind='cubic',bounds_error=False,fill_value='extrapolate')
    def chi(self,z): return self._chi(z)
    def H(self,z): return 100*self.h*np.sqrt(self.Om*(1+z)**3+self.OL)
    def Dz(self,z): return self._D(z)
    def P(self,z,k):   # k in 1/Mpc, P in Mpc^3 ; linear, no smoothing
        kh=k/self.h; return np.exp(self._lP(np.log(kh)))/self.h**3*self.Dz(z)**2

# ---------------- notebook CAMB theory ---------------------------
def camb_results(cd, kmax=1.0, mnu=None, zmax=3.8, nz=250):
    h=cd['h']; p=camb.CAMBparams()
    kw=dict(H0=100*h,ombh2=cd['omega_B']*h**2,omch2=(cd['omega_M']-cd['omega_B'])*h**2)
    if mnu is None: mnu=cd.get('mnu',None)
    if mnu is not None: kw['mnu']=mnu
    p.set_cosmology(**kw); p.set_dark_energy(w=cd['w']); p.InitPower.set_params(ns=cd['ns'],As=cd['As'])
    zpk=(np.exp(np.log(zmax+1)*np.linspace(0,1,nz))-1)[::-1]
    p.set_matter_power(redshifts=zpk,kmax=kmax,nonlinear=False)
    return camb.get_results(p)

class CambCosmo:
    def __init__(self, cd, kmax=1.0, mnu=None, extrap_kmax=None):
        self.r=camb_results(cd,kmax,mnu); self.h=cd['h']; self.Om=cd['omega_M']
        self.pk=self.r.get_matter_power_interpolator(nonlinear=False,hubble_units=False,k_hunit=False,
                    var1='delta_tot',var2='delta_tot',extrap_kmax=extrap_kmax)
        self.kmax=kmax
    def chi(self,z): return self.r.comoving_radial_distance(z)
    def H(self,z): return self.r.hubble_parameter(z)
    def P(self,z,k): return self.pk.P(z,k,grid=False)

# ---------------- Limber ------------------------------------------
def kernels(cosmo, z, chi_star, nz=None, bz=None, Om=None):
    chi=cosmo.chi(z); H=cosmo.H(z); a=1/(1+z); Om=cosmo.Om if Om is None else Om
    Wk=1.5*(100*cosmo.h/C)**2*Om*chi*(chi_star-chi)/chi_star/a
    Wg=None
    if nz is not None:
        n=np.nan_to_num(nz(z)); b=np.nan_to_num(bz(z)); Wg=b*n*H/C
    return chi,H,Wk,Wg

def limber(cosmo, ells, z, WA, WB, Pfunc=None, smooth_kp=None, kcut=None):
    chi=cosmo.chi(z); H=cosmo.H(z)
    k=(ells[None,:]+0.5)/chi[:,None]
    zz=np.broadcast_to(z[:,None],k.shape)
    P=cosmo.P(zz.ravel(),k.ravel()).reshape(k.shape) if Pfunc is None else Pfunc(zz,k)
    if kcut is not None: P=np.where(k<kcut,P,0.)
    if smooth_kp is not None: P=P*np.exp(-(k/smooth_kp)**2)
    I=(C/H)[:,None]*WA[:,None]*WB[:,None]/chi[:,None]**2*P
    return simpson(I,x=z,axis=0)
