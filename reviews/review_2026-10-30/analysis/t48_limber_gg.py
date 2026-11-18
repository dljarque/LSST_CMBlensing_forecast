# Exact/Limber ratio for C_ell^gg, all samples, dense ell (30-300); one Bessel matrix per ell shared by all samples.
import sys
from t46_limber import *
L=[30,34,38,42,46,50,55,60,66,72,80,88,96,105,115,125,140,160,180,200,230,260,300]
W={nm:Wg_(cat,zg) for nm,cat in SAMPLES.items()}
f=np.exp(-(kk_/KP['gg'])**2)*P0*kk_**2
res={nm:[] for nm in W}
for l in L:
    J=spherical_jn(l,np.outer(kk_,chig))
    for nm,Wg in W.items():
        I=simpson(J*(Wg*Dg)[None,:],x=chig,axis=1); ex=2/np.pi*simpson(f*I*I,x=kk_)
        res[nm].append(ex/limber(Wg,Wg,KP['gg'],[l])[0])
    del J
    print(l,' '.join(f'{nm}:{res[nm][-1]:.4f}' for nm in W),flush=True)
json.dump({'L':L,**res},open('t48_limber_gg.json','w'))
