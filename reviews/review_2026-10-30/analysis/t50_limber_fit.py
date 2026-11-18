# Refit the ten-box averages with the gg theory multiplied by the exact/Limber ratio (t48); kg, kk unchanged (Limber good to <=0.1%).
from nbpipe import *
R48=json.load(open('t48_limber_gg.json')); L=np.array(R48['L'])
out={}
for nm,cat in SAMPLES.items():
    r=np.array(R48[nm]); corr=np.interp(ells,L,r); hi=ells>L[-1]; corr[hi]=1-(1-r[-1])*L[-1]/ells[hi]
    Ds=[load(b,cat) for b in BOXES]; cg,ck=theory(Ds[0]['nz'],cat)
    gm=np.mean([d['gg'] for d in Ds],0); km=np.mean([d['kg'] for d in Ds],0); snm=np.mean([d['sn'] for d in Ds])
    S=json.load(open(f't47_{nm}_colore.json'))['cuts']; out[nm]={}
    for rmin in [150,120,100,80,70,60,50,40,30]:
        if str(rmin) not in S or S[str(rmin)]['nb']<7: continue
        m=mask_of(lmax_of(rmin,Ds[0]['zmean']))
        x,e,c2,nd,_=fit_joint(gm,km,cg*corr,ck,snm,m,N=10)
        a=S[str(rmin)]['avg']; out[nm][rmin]=[x[1],e[1],c2,nd]
        print(f"{nm} r={rmin}: Limber {100*(a/s0-1):+.2f}% ({(a-s0)/e[1]:+.2f}σ10) -> non-Limber gg {100*(x[1]/s0-1):+.2f}% ({(x[1]-s0)/e[1]:+.2f}σ10, single {(x[1]-s0)/e[1]/np.sqrt(10):+.2f}) shift {100*(x[1]-a)/s0:+.2f}% χ2 {c2:.1f}/{nd}",flush=True)
json.dump(out,open('t50_limber_fit.json','w'))
