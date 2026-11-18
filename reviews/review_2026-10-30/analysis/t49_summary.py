# Summary numbers for the doc from t47_*_colore.json (+ camb averages) and chi2 at fiducial sigma8 for QSO/LSST.
from nbpipe import *
from scipy.optimize import minimize as mn
S={}; 
for nm in SAMPLES: S[nm]=json.load(open(f't47_{nm}_colore.json')); 
C={nm:json.load(open(f't47_{nm}_camb.json')) for nm in SAMPLES}
print('== safe cuts (ten-box: |avg-fid|<=e10 ; single: |avg-fid|<=e10*sqrt(10)); smallest within, all smaller outside')
def safe(nm,fac,src=S):
    cuts=src[nm]['cuts']; rs=sorted(int(r) for r in cuts if cuts[r]['nb']>=7)
    ok={r:abs(cuts[str(r)]['avg']-s0)<=fac*cuts[str(r)]['eavg'] for r in rs}
    cand=[r for r in rs if ok[r] and all(not ok[q] for q in rs if q<r)]
    return min(cand) if cand else None
for nm in SAMPLES:
    print(nm,'ten-box',safe(nm,1),'single',safe(nm,np.sqrt(10)),'| notebook-as-is (camb):',safe(nm,1,C),safe(nm,np.sqrt(10),C))
print('== scatter/Knox at 100/60/40; paired 60->40; avg vs mean per-box')
for nm in SAMPLES:
    c=S[nm]['cuts']; rr=[c[r]['scatter']/np.mean(c[r]['e']) for r in ('100','60','40')]
    d=100*(np.array(c['40']['s8'])-np.array(c['60']['s8']))/s0
    dm=max(abs(c[r]['avg']-np.mean(c[r]['s8'])) for r in c)
    nl=100*(c['40']['nl']-c['60']['nl'])/s0
    print(f"{nm}: ratio {rr[0]:.2f} {rr[1]:.2f} {rr[2]:.2f} | paired 60->40 {d.mean():+.2f} ± {d.std(ddof=1)/np.sqrt(10):.2f}% (noiseless {nl:+.2f}%) | max|avg-mean per-box| {dm:.4f}")
print('== closure: observed avg vs noiseless (%), and camb-colore shift')
for nm in SAMPLES:
    c=S[nm]['cuts']
    for r in ('80','60','50','40'):
        o=100*(c[r]['avg']/s0-1); n=100*(c[r]['nl']/s0-1); e=100*c[r]['eavg']/s0; cb=100*(C[nm]['cuts'][r]['avg']-c[r]['avg'])/s0
        print(f"{nm} r={r}: obs {o:+.2f} ± {e:.2f}  noiseless {n:+.2f}  obs-nl {o-n:+.2f} ({(o-n)/e:+.1f}σ10) | notebook-as-is shift {cb:+.2f}%")
print('== chi2 best vs fiducial, ten-box average')
for nm,rmins in [('QSO',[60,40]),('LSST',[60,40]),('LRG',[80,40]),('ELG',[50,40]),('BGS',[120,40])]:
    cat=SAMPLES[nm]; Ds=[load(b,cat) for b in BOXES]; cg,ck=theory(Ds[0]['nz'],cat)
    gm=np.mean([d['gg'] for d in Ds],0); km=np.mean([d['kg'] for d in Ds],0); snm=np.mean([d['sn'] for d in Ds])
    for rmin in rmins:
        m=mask_of(lmax_of(rmin,Ds[0]['zmean'])); x,e,c2,nd,f=fit_joint(gm,km,cg,ck,snm,m,N=10)
        g=lambda p: f([p[0],s0,p[1]]); r=mn(g,[x[0],x[2]],method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=20000))
        print(f"{nm} r={rmin}: best s8 {x[1]:.4f} chi2 {c2:.1f}/{nd} | fid chi2 {r.fun:.1f}/{nd+1} | dchi2 {r.fun-c2:.1f} (sqrt {np.sqrt(max(r.fun-c2,0)):.1f}) pull {(x[1]-s0)/e[1]:+.2f}")
