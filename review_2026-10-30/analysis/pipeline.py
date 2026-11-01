from lognormal import *
pw1024=np.load('pw1024.npy')
from scipy.optimize import minimize

class Setup:
    def __init__(self, box, cat, rmin=40., lmin=20, nbins=30, lmax_override=None, nz=250, conv='2pi', chi_at='mean'):
        self.box,self.cat=box,cat; B=load_box(box,cat); self.B=B; cd=B['cd']; self.cd=cd
        self.s8=cd['sigma_8']; self.cm=CambCosmo(cd,kmax=1.0); self.cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat')
        h=cd['h']; kmaxh=(2*np.pi if conv=='2pi' else np.pi)/rmin
        zeff={'mean':B['zmean_nom']}[chi_at] if chi_at=='mean' else chi_at
        self.lmax=int(kmaxh*self.cm.chi(zeff)*h-0.5) if lmax_override is None else lmax_override
        self.lmin=lmin; self.nbins=nbins
        self.ells=np.arange(lmin,self.lmax+1)
        self.Bt=make_hybrid_bin_matrix(self.ells,nbins=nbins,ell_min=lmin,ell_max=self.lmax)
        self.leff=self.Bt@self.ells
        self.zg=np.linspace(1e-4,cd['z_max'],nz)
        self.Bd={xx:make_hybrid_bin_matrix(B[xx]['ell'],nbins=nbins,ell_min=lmin,ell_max=self.lmax) for xx in ['gg','kg']}
        self.d={xx:self.Bd[xx]@B[xx]['cl'] for xx in ['gg','kg']}
        self.pw=pw1024[lmin:self.lmax+1]

    # ---------------- theory -----------------
    def theory(self, kind='notebook', kp=0.327, bias='real', pix='both', zfine=False):
        cos=self.cm if kind=='notebook' else self.cc
        z=np.linspace(1e-4,self.cd['z_max'],3000) if zfine else self.zg
        cs=cos.chi(self.cd['z_out'])
        bz=self.B['bz']
        if bias=='cross':
            gc=self.cat.rsplit('_bin',1)[0]; I=D+'colore_setup/inputs/'
            zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
            self.gf=getattr(self,'gf',None) or GalaxyField(self.cc)
            zz=np.linspace(max(0.01,self.B['zn'].min()-0.05),self.B['zn'].max()+0.05,40)
            bx=np.array([self.gf.b_cross(x,np.interp(x,zi,bi),np.interp(x,zt,tt)) for x in zz])
            bz=interp1d(zz,bx,bounds_error=False,fill_value=np.nan)
        _,_,Wk,Wg=kernels(cos,z,cs,self.B['nz'],bz)
        gg=limber(cos,self.ells,z,Wg,Wg,smooth_kp=kp); kg=limber(cos,self.ells,z,Wk,Wg,smooth_kp=kp)
        wgg = self.pw**2 if pix in ('both','gg') else 1.
        wkg = self.pw**2 if pix=='both' else (self.pw if pix=='gg' else 1.)
        return gg*wgg, kg*wkg

    def hermite_excess(self, nmin=2, kp_extra=None):
        """C_ell^gg from Hermite orders >= nmin (non-linear lognormal+clipping+threshold terms)"""
        gc=self.cat.rsplit('_bin',1)[0]; I=D+'colore_setup/inputs/'
        zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
        self.gf=getattr(self,'gf',None) or GalaxyField(self.cc)
        h=self.cc.h
        zz=np.linspace(max(0.01,self.B['zn'].min()),self.B['zn'].max(),30)
        Pk=[]; 
        for x in zz:
            k,P=self.gf.Pn(x,np.interp(x,zi,bi),np.interp(x,zt,tt),nmin=nmin); Pk.append(np.interp(np.log(k),np.log(k),P))
        Pk=np.array(Pk); kk=k*h; Pk=Pk/h**3
        lP=[interp1d(np.log(kk),P,bounds_error=False,fill_value=0.) for P in Pk]
        def Pf(Z,K):
            out=np.zeros_like(K)
            iz=np.clip(np.searchsorted(zz,Z[:,0])-1,0,len(zz)-2)
            for i in range(K.shape[0]):
                if Z[i,0]<zz[0] or Z[i,0]>zz[-1]: continue
                j=iz[i]; w=(Z[i,0]-zz[j])/(zz[j+1]-zz[j])
                out[i]=(1-w)*lP[j](np.log(K[i]))+w*lP[j+1](np.log(K[i]))
            return out
        z=np.linspace(zz[0],zz[-1],400); cs=self.cc.chi(self.cd['z_out'])
        n=np.nan_to_num(self.B['nz'](z)); Wn=n*self.cc.H(z)/C
        return limber(self.cc,self.ells,z,Wn,Wn,Pfunc=Pf)*self.pw**2

    # ---------------- covariance --------------
    def cov(self, mode='notebook', tgg=None, tkg=None):
        B=self.B; L=B['gg']['ell']; n=len(B['kg']['ell']); sn=B['gg']['sn']
        if mode=='notebook':
            vgg=B['gg']['err']**2; vkg=B['kg']['err']**2
            xc=2/(2*L[:n]+1)*(B['gg']['cl'][:n]+sn)*B['kg']['cl']
        else:  # theory-based Knox, correct factors; tgg,tkg,tkk full-length arrays
            ckk=B['kk']['cl'] if mode=='knox_data' else self._tkk_full
            gg=(B['gg']['cl'] if mode=='knox_data' else tgg)+sn; kg=B['kg']['cl'] if mode=='knox_data' else tkg
            vgg=2*gg**2/(2*L+1); vkg=(kg**2+ckk*gg)/(2*L+1); xc=2*kg*gg/(2*L+1)
        Cgg=self.Bd['gg']@np.diag(vgg)@self.Bd['gg'].T; Ckg=self.Bd['kg']@np.diag(vkg)@self.Bd['kg'].T
        X=self.Bd['gg']@np.diag(xc)@self.Bd['kg'].T
        return Cgg,Ckg,X

    def fit(self, tgg, tkg, cov, which='joint'):
        Cgg,Ckg,X=cov; bg=self.Bt@tgg; bk=self.Bt@tkg
        if which=='joint': C=np.block([[Cgg,X],[X.T,Ckg]]); d=np.concatenate([self.d['gg'],self.d['kg']])
        elif which=='gg': C=Cgg; d=self.d['gg']
        else: C=Ckg; d=self.d['kg']
        Ci=np.linalg.inv(C)
        def model(p):
            A,s=p; g=(s/self.s8)**2
            if which=='joint': return np.concatenate([bg*A**2*g,bk*A*g])
            return bg*A**2*g if which=='gg' else bk*A*g
        chi2=lambda p: (d-model(p))@Ci@(d-model(p))
        r=minimize(chi2,[1.,self.s8],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-8))
        # Fisher errors
        eps=1e-4; J=np.array([(model(r.x+e)-model(r.x-e))/(2*eps) for e in np.eye(2)*eps]).T
        F=J.T@Ci@J; err=np.sqrt(np.diag(np.linalg.inv(F)))
        return dict(A=r.x[0],s8=r.x[1],eA=err[0],es8=err[1],chi2=r.fun,ndof=len(d)-2)

def full(arr, S):
    out=np.zeros(len(S.B['gg']['ell'])); out[S.lmin:S.lmax+1]=arr; return out

def kk_theory(S, kind='colore', kp=0.327, ells=None):
    cos=S.cm if kind=='notebook' else S.cc
    ells=S.ells if ells is None else ells
    z=np.linspace(1e-4,S.cd['z_max'],250); cs=cos.chi(S.cd['z_out'])
    _,_,Wk,_=kernels(cos,z,cs); return limber(cos,ells,z,Wk,Wk,smooth_kp=kp)

def cov_theory(S, tgg, tkg, tkk):
    B=S.B; L=B['gg']['ell']; sn=B['gg']['sn']
    gg=full(tgg,S)+sn; kg=full(tkg,S); kk=full(tkk,S)
    vgg=2*gg**2/(2*L+1); vkg=(kg**2+kk*gg)/(2*L+1); xc=2*kg*gg/(2*L+1)
    return (S.Bd['gg']@np.diag(vgg)@S.Bd['gg'].T, S.Bd['kg']@np.diag(vkg)@S.Bd['kg'].T, S.Bd['gg']@np.diag(xc)@S.Bd['kg'].T)
