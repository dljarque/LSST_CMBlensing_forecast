from pipeline import *
import sys
for box in ['test_box_7','ensemble_1']:
    S=Setup(box,'desi_qso_bin0',rmin=40)
    tkk=kk_theory(S,'colore')*S.pw**2
    for kp in [0.30,0.327,0.3366]:
        gx,kx=S.theory('colore',bias='cross',kp=kp); ex=S.hermite_excess()
        r=S.fit(gx+ex,kx,cov_theory(S,gx+ex,kx,tkk))
        print(f'{box} Hermite kp={kp}: s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}')
    gx,kx=S.theory('colore',bias='cross'); ex=S.hermite_excess()
    for fk in [0.987]:
        r=S.fit(gx+ex,kx*fk,cov_theory(S,gx+ex,kx,tkk))
        print(f'{box} Hermite, kappa amplitude x{fk} (from kk deficit): s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}')
    # pixwin on kg: only galaxy pixwin (kappa point-sampled)
    gx1,kx1=S.theory('colore',bias='cross',pix='gg')
    r=S.fit(gx1+ex,kx1,cov_theory(S,gx1+ex,kx1,tkk)); print(f'{box} Hermite, kg with single pixwin: s8={r["s8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}')
    # sigma_G sensitivity
    for fac in [0.95,1.05]:
        S.gf=GalaxyField(S.cc); S.gf.sig2_0*=fac**2
        gx,kx=S.theory('colore',bias='cross'); ex=S.hermite_excess()
        r=S.fit(gx+ex,kx,cov_theory(S,gx+ex,kx,tkk)); print(f'{box} Hermite sigma_G x{fac}: s8={r["s8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}')
    del S.gf
    sys.stdout.flush()
