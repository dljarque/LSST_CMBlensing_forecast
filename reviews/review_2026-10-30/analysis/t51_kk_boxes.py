# kk-only sigma8 fits per box and for the ten-box average vs ell_max (notebook pipeline, v2 kk window, CoLoRe chi*).
from nbpipe import *
N=len(BOXES); KK={b:load(b,'desi_qso')['kk'] for b in BOXES}; KKm=np.mean([KK[b] for b in BOXES],0)
out={}
for lmx in [100,200,300,500,1000]:
    m=mask_of(lmx); rows=[]
    for b in BOXES:
        a,ea,c2,c2f,nd=fit_kk(KK[b],m); rows.append([a*s0,ea*s0,c2,int(nd)])
    a,ea,c2,c2f,nd=fit_kk(KKm,m,N=N)
    s=np.array([r[0] for r in rows]); e=np.array([r[1] for r in rows])
    out[lmx]=dict(nb=int(m.sum()),boxes=rows,avg=[a*s0,ea*s0,c2,c2f,int(nd)],scatter=float(s.std(ddof=1)))
    print(f"lmax={lmx} nb={m.sum()}: avg s8={a*s0:.4f}±{ea*s0:.4f} ({100*(a-1):+.2f}%, {(a-1)/ea:+.1f}σ10) χ2 {c2:.1f}/{nd} fid {c2f:.1f} | mean per-box {s.mean():.4f}, scatter {s.std(ddof=1):.4f} vs Knox {e.mean():.4f} (ratio {s.std(ddof=1)/e.mean():.2f}) | boxes within 1σ: {np.sum(np.abs(s-s0)<=e)}/10, below fid: {np.sum(s<s0)}/10 | per-box pulls {' '.join(f'{(x-s0)/y:+.1f}' for x,y in zip(s,e))}")
json.dump(out,open('t51_kk_boxes.json','w'))
