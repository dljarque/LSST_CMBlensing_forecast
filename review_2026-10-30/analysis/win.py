import numpy as np
rng=np.random.default_rng(1); N=400000
u=rng.normal(size=(N,3)); u/=np.linalg.norm(u,axis=1)[:,None]
a=rng.normal(size=(N,3)); e1=np.cross(u,a); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(u,e1)
psi=rng.uniform(0,2*np.pi,N); n=np.cos(psi)[:,None]*e1+np.sin(psi)[:,None]*e2   # unit k, perpendicular to ray (Limber)
sinc=lambda x: np.sinc(x/np.pi)
P=np.eye(3)[None]-u[:,:,None]*u[:,None,:]           # transverse projector = CoLoRe 'rot'
print(' k*dx | W_gal (cube)   -ln/(kdx)^2 | FD only  -ln/(kdx)^2 | FD+NGP kappa -ln/(kdx)^2 | cross gk -ln/(kdx)^2')
for kdx in [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8]:
    q=kdx*n                                            # k_i dx
    Wc=np.prod(sinc(q/2),axis=1)                       # top-hat cell (amplitude): galaxies placed uniformly; NGP read
    # finite-difference second derivatives (beaming.c get_element), in units of 1/dx^2
    D=-q[:,:,None]*q[:,None,:]*sinc(q)[:,:,None]*sinc(q)[:,None,:]   # off-diagonal: 0.25*(4 neighbours)
    for i in range(3): D[:,i,i]=-(q[:,i]**2)*sinc(q[:,i]/2)**2       # diagonal: phi(+)+phi(-)-2phi
    Kfd=np.einsum('nij,nij->n',P,D); Kc=-kdx**2
    Wfd=Kfd/Kc
    g=np.mean(Wc**2); f=np.mean(Wfd**2); kk=np.mean((Wfd*Wc)**2); x=np.mean(Wfd*Wc*Wc)
    print(f' {kdx:.2f} | {g:.5f}  {-np.log(g)/kdx**2:.4f} | {f:.5f} {-np.log(f)/kdx**2:.4f} | {kk:.5f} {-np.log(kk)/kdx**2:.4f} | {x:.5f} {-np.log(x)/kdx**2:.4f}')
