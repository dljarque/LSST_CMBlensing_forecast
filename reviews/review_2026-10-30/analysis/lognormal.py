import numpy as np, mcfit, warnings
from numpy.polynomial.hermite_e import hermegauss, hermeval
from scipy.interpolate import interp1d
from common import *
warnings.filterwarnings('ignore')

def colore_grid(cc, ngrid=4096, zmax=3.8):
    rmax=cc.chi(zmax)*cc.h          # Mpc/h
    L=2*rmax*(1+2./ngrid); return L, L/ngrid

class GalaxyField:
    """Exact two-point of CoLoRe _BIAS_MODEL_3 + threshold on the lognormal field (Hermite expansion)."""
    def __init__(self, cc, R=2.0, ngrid=4096, Rg2_extra=None, nherm=12):
        self.cc=cc; L,dx=colore_grid(cc,ngrid); self.dx=dx
        self.reff2=R**2+(0.45*dx)**2                       # CoLoRe sigma2_gauss smoothing (Mpc/h)^2
        self.Rxi2=self.reff2 if Rg2_extra is None else R**2+Rg2_extra
        kh=np.logspace(-5,1.6,8192); P=np.exp(cc._lP(np.log(kh)))   # (Mpc/h)^3, z=0, sigma8-normalised
        self.kh=kh; self.P0=P
        self.sig2_0=np.trapezoid(kh**3*P*np.exp(-kh**2*self.reff2)/(2*np.pi**2),np.log(kh))
        Ps=P*np.exp(-kh**2*self.Rxi2); Ps=np.maximum(Ps,1e-300)
        r,xi=mcfit.P2xi(kh,lowring=True)(Ps,extrap=False)
        self.r=r; self.xi0=xi
        self.nherm=nherm
        self.x=np.linspace(-12,12,240001); self.w=np.exp(-self.x**2/2)/np.sqrt(2*np.pi)*(self.x[1]-self.x[0])
    def coeffs(self, z, b, t):
        s=self.cc.Dz(z)*np.sqrt(self.sig2_0)
        dLN=np.exp(s*self.x-0.5*s**2)-1
        f=np.where((dLN<=-1)|(dLN<t),0.,np.maximum(1+b*dLN,0.))
        c=np.array([np.sum(self.w*f*hermeval(self.x,np.eye(n+1)[n])) for n in range(self.nherm)])
        return c, s
    def b_cross(self, z, b, t):
        c,s=self.coeffs(z,b,t); return c[1]/(c[0]*s)
    def Pn(self, z, b, t, nmin=1):
        """returns kh, P_gg(kh) [(Mpc/h)^3] summing Hermite orders n>=nmin (n=1 is linear b_x^2 P_G)"""
        c,s=self.coeffs(z,b,t); rho=self.cc.Dz(z)**2*self.xi0/s**2
        xi=sum((c[n]/c[0])**2*rho**n/np.math.factorial(n) if hasattr(np,'math') else (c[n]/c[0])**2*rho**n/float(np.prod(np.arange(1,n+1))) for n in range(nmin,self.nherm))
        k,P=mcfit.xi2P(self.r,lowring=True)(xi,extrap=False)
        return k,P
