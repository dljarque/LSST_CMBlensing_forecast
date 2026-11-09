# Ten-box (ensemble_1-10) QSO + LSST fits; test_box_11 excluded.
# t30.py takes its reference cosmology from test_box_11; re-run its setup with ensemble_1 (identical in all ten boxes).
import t30
t30._REF='ensemble_1'
exec("""
cd=json.load(open(D+_REF+'/cosmo_params.json')); cc=CoLoReCosmo(cd,D+_REF+'/pk_linear.dat'); h=cd['h']; s0=cd['sigma_8']
z=np.linspace(1e-4,3.8,600); chi=cc.chi(z); H=cc.H(z); cs=cc.chi(cd['z_out'])
k=(ells[None,:]+.5)/chi[:,None]; P=cc.P(np.broadcast_to(z[:,None],k.shape).ravel(),k.ravel()).reshape(k.shape)
_,_,Wk,_=kernels(cc,z,cs)
kk_th=lim(Wk,Wk,KP_KK)
""",t30.__dict__)
from t30 import *

# t43: all five samples on ensemble_1-10; sigma8 vs r_min (data, with/without kappa T(ell)) + noiseless Hermite test with the ensemble n(z)
import sys
from pipeline import Setup
from lognormal import GalaxyField
boxes=[f'ensemble_{i}' for i in range(1,11)]; N=len(boxes)
SAMPLES={'BGS':('desi_bgs_bin0',0.4),'LRG':('desi_lrg_bin0',1.0),'ELG':('desi_elg_bin0',1.2),'QSO':('desi_qso_bin0',2.0),'LSST':('lsst_gal_bin0',1.0)}
RMINS=[150,120,100,80,60,50,40,30]
KK={b:data(b,'desi_qso_bin0')['kk'] for b in boxes}
lb=Bt@ells
def T_of(excl):
    m=np.mean([KK[b] for b in boxes if b!=excl],axis=0)
    r=(Bt@m)/(Bt@kk_th); x=np.log(lb); w=1/((Bt**2)@(2*kk_th**2/(2*ells+1)))*(Bt@kk_th)**2
    c=np.polyfit(x,r,3,w=np.sqrt(w)); return np.polyval(c,np.log(ells))
T_all=T_of(None)
def run(name):
    cat,z0=SAMPLES[name]; gc=cat.rsplit('_bin',1)[0]
    tg,tk=gal_theory('ensemble_1',cat)
    Ds=[data(b,cat) for b in boxes]; TH=[gal_theory(b,cat) for b in boxes]
    gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
    # noiseless truth: linear with mock bias b_x(z) + Hermite n>=2 excess (ensemble_1 n(z)); template = linear with b_x
    S=Setup('ensemble_1',cat,rmin=40,lmin=LMIN,nbins=NB,lmax_override=LMAX)
    I=D+'colore_setup/inputs/'; gf=GalaxyField(cc)
    zi,bi=np.loadtxt(I+f'bias-input/bz-input_{gc}.txt',unpack=True); zt,tt=np.loadtxt(I+f'thresholds/tz_{gc}.txt',unpack=True)
    zz=np.linspace(max(0.02,S.B['zn'].min()),S.B['zn'].max(),40)
    bx=np.array([gf.b_cross(x,np.interp(x,zi,bi),np.interp(x,zt,tt)) for x in zz]); bxf=interp1d(zz,bx,bounds_error=False,fill_value=np.nan)
    S.gf=gf; ex=S.hermite_excess()
    _,_,_,Wgx=kernels(cc,z,cs,S.B['nz'],bxf)
    tgx=lim(Wgx,Wgx,h/2)*pw**2; tkx=lim(Wk,Wgx,h/2)*pw
    res={'z0':z0,'ex_ratio':((Bt@ex)/(Bt@tgx)).tolist(),'lb':lb.tolist(),'cuts':{}}
    for rmin in RMINS:
        lmax=int(2*np.pi/rmin*cc.chi(z0)*h-0.5); mask=edges[1:]<=lmax
        if mask.sum()<3: continue
        s8=[];e=[]
        for dd,(tgb,tkb) in zip(Ds,TH):
            x,err,_,_=fit_joint(dd['gg'],dd['kg'],tgb,tkb,dd['sn'],mask); s8.append(x[1]); e.append(err[1])
        xa,ea,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=N)
        xT,eT,c2T,_=fit_joint(gm,km,tg,tk*np.sqrt(T_all),snm,mask,N=N)
        xn,en,c2n,_=fit_joint(tgx+ex,tkx,tgx,tkx,snm,mask,N=N)
        r=dict(lmax=lmax,nbins=int(mask.sum()),s8=s8,e=e,avg=xa[1],eavg=ea[1],A=xa[0],Asn=xa[2],chi2=c2,nd=nd,avgT=xT[1],chi2T=c2T,noiseless=xn[1],enl=en[1],scatter=float(np.std(s8,ddof=1)))
        res['cuts'][rmin]=r
        print(f'{name} rmin={rmin} lmax={lmax} nb={mask.sum()}: avg {xa[1]:.4f}±{ea[1]:.4f} ({100*(xa[1]/s0-1):+.2f}%, {(xa[1]-s0)/ea[1]:+.1f}σ) withT {xT[1]:.4f} ({(xT[1]-s0)/ea[1]:+.1f}σ) | noiseless {100*(xn[1]/s0-1):+.2f}% ({(xn[1]-s0)/ea[1]:+.1f}σ10) | A={xa[0]:.3f} Asn={xa[2]:.3f} χ2={c2:.0f}/{nd} | scatter/Knox {np.std(s8,ddof=1)/np.mean(e):.2f}',flush=True)
    json.dump(res,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),f't43_{name}.json'),'w'),default=float)
if __name__=='__main__':
    run(sys.argv[1])
