import t30
from t30 import *
for kp in [0.29574,0.2933]:
    t30.kk_th=lim(Wk,Wk,kp)
    allkk=[data(b,'desi_qso_bin0')['kk'] for b in boxes]
    per=[t30.fit_kk(d,edges[1:]<=1000.01) for d in allkk]
    dm=np.mean(allkk,axis=0)
    for lm in [300,500,1000.01]:
        s,e,c2,nd=t30.fit_kk(dm,edges[1:]<=lm,N=4)
        print(f'kp_kk={kp}: mean-of-4 s8/fid lmax={int(lm)}: {s:.4f}±{e:.4f} χ2={c2:.0f}/{nd}')
    print('   per box lmax1000:',[round(p[0],4) for p in per])
