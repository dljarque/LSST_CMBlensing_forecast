from pipeline import *
import sys
from t5_mech import fitA
box,cat='test_box_11','desi_qso_bin0'
S=Setup(box,cat,rmin=40); h=S.cd['h']; kp=h/2
S.cm=CambCosmo(S.cd,kmax=20.)
print('sigma8 fid',S.s8,' CoLoRe-style sigma8 of file %.5f'%S.cc.sig8_file, ' CAMB(mnu=%s) s8 %.5f'%(S.cd.get('mnu'),S.cm.r.get_sigma8_0()))
tg,tk=S.theory('notebook',kp=kp); cov=S.cov('notebook')
r=S.fit(tg,tk,cov); print(f'NOTEBOOK replica (data cov, fixed err_kg): s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.2f}σ) A={r["A"]:.4f} chi2={r["chi2"]:.1f}/{r["ndof"]}')
for nm,(g,k) in [('notebook b_real',(tg,tk)),('CoLoRe b_real',S.theory('colore',kp=kp)),('CoLoRe b_x',S.theory('colore',bias='cross',kp=kp))]:
    Ag=fitA(S,g,'gg',cov); Ak=fitA(S,k,'kg',cov); print(f'  s8 fixed, {nm:16s}: A_gg={Ag[0]:.4f}±{Ag[1]:.4f} A_kg={Ak[0]:.4f}±{Ak[1]:.4f} ratio^2={Ag[0]**2/Ak[0]**2:.4f}')
sys.stdout.flush()
for rmin in [100,60,40,30,20]:
    S=Setup(box,cat,rmin=rmin); S.cm=CambCosmo(S.cd,kmax=20.)
    tg,tk=S.theory('notebook',kp=kp); gx,kx=S.theory('colore',bias='cross',kp=kp); ex=S.hermite_excess()
    tkk=kk_theory(S,'colore',kp=kp)*S.pw**2
    out=f'rmin={rmin:3d} lmax={S.lmax:4d}'
    for nm,(g,k) in [('NB',(tg,tk)),('Herm',(gx+ex,kx)),('Herm*0.987k',(gx+ex,kx*0.987))]:
        r=S.fit(g,k,cov_theory(S,g,k,tkk)); out+=f' | {nm}: {r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}'
    print(out); sys.stdout.flush()
