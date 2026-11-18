# kk data/theory, ten ensemble boxes, notebook pipeline with v2 kp_kk (19/84); background from NB_BG (colore default).
from nbpipe import *
N=len(BOXES); KK=np.array([load(b,'desi_qso')['kk'] for b in BOXES]); KKm=KK.mean(0)
print(f"[{BG}] dx={dx:.4f} Mpc/h chi*={chi_star:.1f} Mpc  kp[h/Mpc]: gg {KP['gg']/h:.4f} kg {KP['kg']/h:.4f} kk {KP['kk']/h:.4f}")
out={'bg':BG}
ck=CKK; Rr=(B@KK.T).T/(B@ck)
w=(B@ck)**2/np.diag(B@np.diag(2*ck**2/(2*ells+1))@B.T)
for lo,hi in [(30,100),(100,200),(200,300),(300,500),(500,1000)]:
    m=(lb>=lo)&(lb<hi); v=(Rr[:,m]*w[m]).sum(1)/w[m].sum()
    print(f'  {lo}-{hi}: {v.mean():.4f} ± {v.std(ddof=1)/np.sqrt(N):.4f}  box range {v.min():.3f}-{v.max():.3f}  below 1: {(v<1).sum()}/10')
    out[f'{lo}-{hi}']=[v.mean(),v.std(ddof=1)/np.sqrt(N),v.tolist()]
for lmx in [300,500,1000]:
    a,ea,c2,c2f,nd=fit_kk(KKm,mask_of(lmx),N=N)
    print(f'  kk-only lmax={lmx}: s8/fid={a:.4f}±{ea:.4f}  chi2_min={c2:.1f}/{nd}  chi2_fid={c2f:.1f}/{nd+1}  dchi2={c2f-c2:.1f}')
    out[f'fit{lmx}']=[a,ea,c2,c2f,int(nd)]
out['bins']=[lb.tolist(),Rr.mean(0).tolist(),(Rr.std(0,ddof=1)/np.sqrt(N)).tolist()]
json.dump(out,open(f't45_kk_{BG}.json','w'))
