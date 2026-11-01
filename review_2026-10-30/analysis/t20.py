from pipeline import *
S=Setup('test_box_11','desi_qso_bin0',rmin=40); h=S.cd['h']; L,dx=colore_grid(S.cc)
kgg=h/np.sqrt(4+dx**2/12); kkg=h/np.sqrt(4+dx**2/6); kkk=h/np.sqrt(4+dx**2/4)
ex=S.hermite_excess(); tkk=kk_theory(S,'colore',kp=kkk)
for b in ['real','cross']:
    g,_=S.theory('colore',bias=b,kp=kgg,pix='both'); _,k=S.theory('colore',bias=b,kp=kkg,pix='gg')
    r=S.fit(g+ex,k,cov_theory(S,g+ex,k,tkk)); print(b,'excess scaled by A^2 with the linear part:',f'{r["s8"]:.4f} ({(r["s8"]-S.s8)/r["es8"]:+.1f}σ) A={r["A"]:.3f}')
    # excess not scaled by A
    Cgg,Ckg,X=cov_theory(S,g+ex,k,tkk); C=np.block([[Cgg,X],[X.T,Ckg]]); Ci=np.linalg.inv(C); d=np.concatenate([S.d['gg'],S.d['kg']])
    from scipy.optimize import minimize
    bg,be,bk=S.Bt@g,S.Bt@ex,S.Bt@k
    m=lambda p: np.concatenate([(p[1]/S.s8)**2*(p[0]**2*bg+be),(p[1]/S.s8)**2*p[0]*bk])
    rr=minimize(lambda p:(d-m(p))@Ci@(d-m(p)),[1,S.s8],method='Nelder-Mead'); print(b,'excess fixed (not scaled by A):',f'{rr.x[1]:.4f} ({(rr.x[1]-S.s8)/r["es8"]:+.1f}σ) A={rr.x[0]:.3f}')
