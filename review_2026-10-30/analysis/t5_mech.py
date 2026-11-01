from pipeline import *
import sys
from scipy.optimize import minimize_scalar
def fitA(S,t,which,cov):
    Cgg,Ckg,X=cov; C=Cgg if which=='gg' else Ckg; Ci=np.linalg.inv(C); d=S.d[which]; m=S.Bt@t
    if which=='kg': A=(m@Ci@d)/(m@Ci@m); e=1/np.sqrt(m@Ci@m); return A,e
    f=lambda A:(d-A**2*m)@Ci@(d-A**2*m); A=minimize_scalar(f,bounds=(0.5,1.5),method='bounded').x
    e=1/np.sqrt(4*A**2*(m@Ci@m)); return A,e
if __name__=="__main__":
 for box,cat in [('test_box_7','desi_qso_bin0'),('ensemble_1','desi_qso_bin0'),('test_box_7','lsst_gal_bin0')]:
     S=Setup(box,cat); cov=S.cov('notebook')
     tg,tk=S.theory('notebook')
     print(f'\n== {box} {cat} (lmax={S.lmax}), sigma8 fixed to fiducial')
     for nm,(g,k) in [('notebook (CAMB, b_real)',(tg,tk)),('CoLoRe-cosmo, b_real',S.theory('colore')),('CoLoRe-cosmo, b_x exact',S.theory('colore',bias='cross'))]:
         Ag=fitA(S,g,'gg',cov); Ak=fitA(S,k,'kg',cov)
         print(f'  {nm:28s} A_gg={Ag[0]:.4f}±{Ag[1]:.4f}   A_kg={Ak[0]:.4f}±{Ak[1]:.4f}   A_gg^2/A_kg^2={Ag[0]**2/Ak[0]**2:.4f}')
     g,k=S.theory('colore',bias='cross'); ex=S.hermite_excess()
     print('  Hermite n>=2 excess / linear(b_x) gg at ell=',[int(x) for x in S.ells[[0,30,100,200,300,-1]]], np.round((ex/g)[[0,30,100,200,300,-1]],4))
     Ag=fitA(S,g+ex,'gg',cov); print(f'  b_x + Hermite excess:        A_gg={Ag[0]:.4f}±{Ag[1]:.4f}')
     for nm,(gg,kg) in [('notebook',(tg,tk)),('CoLoRe-cosmo,b_x,linear',(g,k)),('CoLoRe-cosmo,b_x,+Hermite',(g+ex,k))]:
         r=S.fit(gg,kg,cov); print(f'  JOINT {nm:28s} s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.2f}σ) A={r["A"]:.4f}±{r["eA"]:.4f} chi2={r["chi2"]:.1f}/{r["ndof"]}')
     sys.stdout.flush()
