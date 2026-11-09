from pipeline import *
for box in ['test_box_7','ensemble_1']:
  for lmin in [20,30,50]:
    S=Setup(box,'desi_qso_bin0',rmin=40,lmin=lmin,nbins=25 if lmin>20 else 30)
    tkk=kk_theory(S,'colore')*S.pw**2; gx,kx=S.theory('colore',bias='cross'); ex=S.hermite_excess()
    r=S.fit(gx+ex,kx,cov_theory(S,gx+ex,kx,tkk)); tg,tk=S.theory('notebook'); r0=S.fit(tg,tk,S.cov('notebook'))
    print(f'{box} lmin={lmin}: Hermite s8={r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) χ2={r["chi2"]:.0f}/{r["ndof"]} | notebook s8={r0["s8"]:.4f} ({(r0["s8"]-S.s8)/r0["es8"]:+.1f}σ)')
