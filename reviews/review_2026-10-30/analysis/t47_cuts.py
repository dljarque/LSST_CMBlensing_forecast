# r_min scan per sample on ensemble_1-10 with the notebook pipeline (nbpipe: v2 windows, notebook binning/cuts/cov).
# Background from NB_BG (default 'colore'). No kappa transfer function. Noiseless Hermite test with the same windows.
# usage: NB_BG=colore python t47_cuts.py LRG [noiseless]
from nbpipe import *
name=sys.argv[1]; cat=SAMPLES[name]; do_nl=len(sys.argv)>2; N=len(BOXES)
RMINS=[200,150,120,100,80,70,60,50,40,30]
Ds=[load(b,cat) for b in BOXES]; TH=[theory(d['nz'],cat) for d in Ds]
gm=np.mean([d['gg'] for d in Ds],0); km=np.mean([d['kg'] for d in Ds],0); snm=np.mean([d['sn'] for d in Ds])
cg,ck=TH[0]; zmean=Ds[0]['zmean']
if do_nl:
    from pipeline import Setup
    from lognormal import GalaxyField
    from scipy.interpolate import interp1d
    S=Setup('ensemble_1',cat+'_bin0',rmin=40,lmin=LMIN,nbins=NB,lmax_override=LMAX)
    I=D+'colore_setup/inputs/'; gf=GalaxyField(CC)
    zi,bi=np.loadtxt(I+f'bias-input/bz-input_{cat}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{cat}.txt',unpack=True)
    zz=np.linspace(max(0.02,S.B['zn'].min()),S.B['zn'].max(),40)
    bx=np.array([gf.b_cross(x,np.interp(x,zi,bi),np.interp(x,zt,tt)) for x in zz]); bxf=interp1d(zz,bx,bounds_error=False,fill_value=0.)
    S.gf=gf; ex=S.hermite_excess()                      # includes pw^2
    nz=Ds[0]['nz']; Wx=np.nan_to_num(bxf(z))*np.nan_to_num(nz(z))*Hz/C_KMS
    cgx=lim(Wx,Wx,KP['gg']); ckx=lim(Wk,Wx,KP['kg'])
    tdg=cgx*pw**2+ex; tdk=ckx*pw
res={'name':name,'bg':BG,'zmean':zmean,'cuts':{}}
for rmin in RMINS:
    lmax=lmax_of(rmin,zmean); m=mask_of(lmax)
    if m.sum()<4: continue
    x,e,c2,nd,_=fit_joint(gm,km,cg,ck,snm,m,N=N)
    r=dict(lmax=lmax,nb=int(m.sum()),avg=x[1],eavg=e[1],A=x[0],Asn=x[2],chi2=c2,nd=nd)
    if BG=='colore':
        s8=[];es=[]
        for d,(a,b) in zip(Ds,TH):
            xx,ee,_,_,_=fit_joint(d['gg'],d['kg'],a,b,d['sn'],m); s8.append(xx[1]); es.append(ee[1])
        r.update(s8=s8,e=es,scatter=float(np.std(s8,ddof=1)))
    if do_nl:
        xn,en,_,_,_=fit_joint(tdg,tdk,cgx,ckx,snm,m,N=N); r.update(nl=xn[1])
    res['cuts'][rmin]=r
    s=f"{name} [{BG}] rmin={rmin} lmax={lmax} nb={m.sum()}: avg {x[1]:.4f}±{e[1]:.4f} ({100*(x[1]/s0-1):+.2f}%, {(x[1]-s0)/e[1]:+.2f}σ10, {(x[1]-s0)/(e[1]*np.sqrt(N)):+.2f}σ1) A={x[0]:.3f} Asn={x[2]:.3f} χ2={c2:.1f}/{nd}"
    if 's8' in r: s+=f" | per-box within 1σ1: {np.sum(np.abs(np.array(r['s8'])-s0)<=np.array(r['e']))}/10, scatter/Knox {r['scatter']/np.mean(r['e']):.2f}"
    if do_nl: s+=f" | noiseless {100*(xn[1]/s0-1):+.2f}% ({(xn[1]-s0)/e[1]:+.2f}σ10)"
    print(s,flush=True)
json.dump(res,open(f't47_{name}_{BG}.json','w'),default=float)
