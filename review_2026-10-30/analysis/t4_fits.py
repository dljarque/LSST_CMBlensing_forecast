from pipeline import *
import sys
for box,cat in [('test_box_7','desi_qso_bin0'),('ensemble_1','desi_qso_bin0'),('test_box_10','desi_qso_bin0'),('test_box_7','lsst_gal_bin0')]:
    S=Setup(box,cat)
    tg,tk=S.theory('notebook'); cov=S.cov('notebook')
    r=S.fit(tg,tk,cov); 
    print(f'{box:12s} {cat:14s} lmax={S.lmax}  NOTEBOOK REPLICA: s8={r["s8"]:.4f}±{r["es8"]:.4f} (fid {S.s8:.4f}, shift {(r["s8"]-S.s8)/r["es8"]:+.2f}σ)  A={r["A"]:.4f}±{r["eA"]:.4f}  chi2={r["chi2"]:.1f}/{r["ndof"]}')
    sys.stdout.flush()
