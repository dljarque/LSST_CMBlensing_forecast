from pipeline import *
import sys
for box in ['test_box_11']:
  for rmin in [100,60,40,30,20]:
    S=Setup(box,'desi_qso_bin0',rmin=rmin); h=S.cd['h']; L,dx=colore_grid(S.cc)
    kgg=h/np.sqrt(4+dx**2/12); kkg=h/np.sqrt(4+0.1547*dx**2); kkk=h/np.sqrt(4+0.2261*dx**2)
    g,_=S.theory('colore',bias='cross',kp=kgg,pix='both'); _,k=S.theory('colore',bias='cross',kp=kkg,pix='gg')
    ex=S.hermite_excess(); tkk=kk_theory(S,'colore',kp=kkk)
    # previous setup for comparison
    g0,k0=S.theory('colore',bias='cross',kp=h/2)
    r0=S.fit(g0+ex,k0,cov_theory(S,g0+ex,k0,kk_theory(S,'colore',kp=h/2)*S.pw**2))
    r=S.fit(g+ex,k,cov_theory(S,g+ex,k,tkk))
    print(f'{box:11s} rmin={rmin:3d} lmax={S.lmax:4d} | Herm, kp=h/R, pw^2 both: {r0["s8"]:.4f} ({(r0["s8"]-S.s8)/r0["es8"]:+.1f}σ) | Herm + cell windows, pw_g only on kg: {r["s8"]:.4f}±{r["es8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f} χ2={r["chi2"]:.0f}/{r["ndof"]}')
    sys.stdout.flush()
