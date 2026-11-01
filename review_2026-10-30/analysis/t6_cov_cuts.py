from pipeline import *
import sys
res={}
for box,cat in [('test_box_7','desi_qso_bin0'),('ensemble_1','desi_qso_bin0')]:
    print(f'\n===== {box} {cat}')
    for rmin in [100,80,60,40,30,20,15]:
        S=Setup(box,cat,rmin=rmin)
        tg,tk=S.theory('notebook'); gx,kx=S.theory('colore',bias='cross'); ex=S.hermite_excess()
        tkk=kk_theory(S,'colore')*S.pw**2
        line=f'rmin={rmin:3d} lmax={S.lmax:4d} |'
        for nm,(g,k),cv in [('NB/NBcov',(tg,tk),S.cov('notebook')),('NB/thcov',(tg,tk),cov_theory(S,tg,tk,tkk)),('Herm/thcov',(gx+ex,kx),cov_theory(S,gx+ex,kx,tkk))]:
            r=S.fit(g,k,cv); line+=f' {nm}: s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}/{r["ndof"]} |'
        print(line); sys.stdout.flush()
