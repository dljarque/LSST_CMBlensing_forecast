from pipeline import *
from scipy.optimize import minimize
import sys
box,cat='test_box_11','desi_qso_bin0'
def build(S, herm):
    cos=S.cc; z=np.linspace(1e-4,3.8,600); cs=cos.chi(S.cd['z_out'])
    _,_,Wk,Wg=kernels(cos,z,cs,S.B['nz'],S.B['bz'])
    chi=cos.chi(z); H=cos.H(z); k=(S.ells[None,:]+.5)/chi[:,None]; zz=np.broadcast_to(z[:,None],k.shape)
    P=cos.P(zz.ravel(),k.ravel()).reshape(k.shape); w=simpson(np.eye(len(z)),x=z,axis=0) if False else None
    base=lambda A,B: (C/H)[:,None]*A[:,None]*B[:,None]/chi[:,None]**2*P
    sw=np.array([simpson(np.eye(len(z))[i],x=z) for i in range(len(z))])
    ex=S.hermite_excess()/S.pw**2 if herm else 0*S.ells
    return dict(k=k,sw=sw,gg=base(Wg,Wg),kg=base(Wk,Wg),kk=base(Wk,Wk),ex=ex)
if __name__=="__main__":
  for rmin in [40,20]:
      S=Setup(box,cat,rmin=rmin); h=S.cd['h']; s0=S.s8
      for herm in [False,True]:
          T=build(S,herm)
          cl=lambda kp,X: T['sw']@(T[X]*np.exp(-(T['k']/kp)**2))
          g0=cl(h/2,'gg'); k0=cl(h/2,'kg'); kk0=cl(h/2,'kk')
          cov=cov_theory(S,(g0+T['ex'])*S.pw**2,k0*S.pw**2,kk0*S.pw**2); Cgg,Ckg,X=cov
          Ci=np.linalg.inv(np.block([[Cgg,X],[X.T,Ckg]])); d=np.concatenate([S.d['gg'],S.d['kg']])
          def model(p,free):
              A,s=p[0],p[1]; kg_=p[2] if free else h/2; kk_=p[3] if free else h/2; f=(s/s0)**2
              return np.concatenate([S.Bt@((cl(kg_,'gg')*A**2+T['ex'])*f*S.pw**2), S.Bt@(cl(kk_,'kg')*A*f*S.pw**2)])
          for free in [False,True]:
              x0=[1,s0,h/2,h/2] if free else [1,s0]
              fn=lambda p:(d-model(p,free))@Ci@(d-model(p,free)) if (not free or (0.05<p[2]<5 and 0.05<p[3]<5)) else 1e12
              r=minimize(fn,x0,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-8,maxiter=40000))
              eps=1e-4; J=np.array([(model(r.x+e,free)-model(r.x-e,free))/(2*eps) for e in np.eye(len(r.x))*eps]).T
              F=J.T@Ci@J; Cv=np.linalg.inv(F); err=np.sqrt(np.diag(Cv)); corr=Cv/np.outer(err,err)
              extra=f" kp_gg={r.x[2]:.3f}±{err[2]:.3f} kp_kg={r.x[3]:.3f}±{err[3]:.3f}  corr(s8,kp_gg)={corr[1,2]:+.2f} corr(s8,kp_kg)={corr[1,3]:+.2f} corr(A,kp_gg)={corr[0,2]:+.2f}" if free else ''
              print(f"rmin={rmin} {'Hermite' if herm else 'linear '} kp {'FREE ' if free else 'h/R  '}: s8={r.x[1]:.4f}±{err[1]:.4f} ({(r.x[1]-s0)/err[1]:+.1f}σ) A={r.x[0]:.3f} χ2={r.fun:.1f}{extra}")
              sys.stdout.flush()
