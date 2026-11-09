import t30
from t30 import *
from lognormal import GalaxyField, colore_grid
# A. CoLoRe P_ll vs my theory
L,dx=colore_grid(cc); rsm2=4+1.88*dx**2/6
for zz in [0.0,1.0,2.0,3.0]:
    kh,Ptt,Ptl,Pll=np.loadtxt(D+f'ensemble_4/results/out_pk_srcs_pop4_z{zz:.3f}.txt',unpack=True)
    mine=cc.P(zz,kh*h)*h**3*np.exp(-rsm2*kh**2)
    r=Pll/mine; sel=[np.argmin(abs(kh-x)) for x in [0.001,0.01,0.05,0.1,0.2]]
    print(f'z={zz}: P_ll(CoLoRe)/P_mine at k=0.001,0.01,0.05,0.1,0.2 h/Mpc:',' '.join(f'{r[i]:.5f}' for i in sel))
# B. kk 6 boxes
boxes6=['test_box_11','ensemble_1','ensemble_2','ensemble_3','ensemble_4','ensemble_5']
allkk=[data(b,'desi_qso_bin0')['kk'] for b in boxes6]
for b,dk in zip(boxes6,allkk):
    row=[]
    for a_,b_ in [(30,100),(100,200),(200,300)]:
        m=(ells>=a_)&(ells<b_); w=2*ells[m]+1; row.append(f'{np.sum(w*dk[m])/np.sum(w*kk_th[m]):.4f}')
    print(b,' '.join(row), 's8/fid(lmax300)=%.4f'%t30.fit_kk(dk,edges[1:]<=300)[0])
dm=np.mean(allkk,axis=0); row=[]
for a_,b_ in [(30,100),(100,200),(200,300),(300,500),(500,1000)]:
    m=(ells>=a_)&(ells<b_); w=2*ells[m]+1; row.append(f'{np.sum(w*dm[m])/np.sum(w*kk_th[m]):.4f}±{np.sqrt(2/np.sum(w)/6):.4f}')
print('mean6',' '.join(row))
for lm in [300,500,1000.01]:
    s,e,c2,nd=t30.fit_kk(dm,edges[1:]<=lm,N=6); print(f'mean6 kk s8/fid lmax={int(lm)}: {s:.4f}±{e:.4f} χ2={c2:.0f}/{nd}')
sys.stdout.flush()
# C. QSO joint, 6 boxes
cat='desi_qso_bin0'; tg,tk=gal_theory('test_box_11',cat); chi_h=cc.chi(2.0)*h
for rmin in [60,40]:
    lmax=int(2*np.pi/rmin*chi_h-0.5); mask=edges[1:]<=lmax; Ds=[]; s8s=[]; es=[]
    for b in boxes6:
        dd=data(b,cat); Ds.append(dd); tgb,tkb=gal_theory(b,cat)
        x,err,c2,nd=fit_joint(dd['gg'],dd['kg'],tgb,tkb,dd['sn'],mask); s8s.append(x[1]); es.append(err[1])
        print(f'  {b:12s} rmin={rmin}: s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ)')
    gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
    x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=6)
    print(f'QSO rmin={rmin} lmax={lmax}: per-box scatter {np.std(s8s,ddof=1):.4f} vs Knox {np.mean(es):.4f}; AVERAGED(6): s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd}')
    sys.stdout.flush()
