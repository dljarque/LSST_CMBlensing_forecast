from t30 import *
boxes6=['test_box_11','ensemble_1','ensemble_2','ensemble_3','ensemble_4','ensemble_5']
KK={b:data(b,'desi_qso_bin0')['kk'] for b in boxes6}
lb=Bt@ells
def T_of(excl):
    m=np.mean([KK[b] for b in boxes6 if b!=excl],axis=0)
    r=(Bt@m)/(Bt@kk_th); x=np.log(lb); w=1/((Bt**2)@(2*kk_th**2/(2*ells+1)))*(Bt@kk_th)**2
    c=np.polyfit(x,r,3,w=np.sqrt(w)); return np.polyval(c,np.log(ells))
T_all=T_of(None)
print('smooth T(ell) all 6 at ell=30,50,100,200,300,500,1000:',' '.join(f'{np.interp(l,ells,T_all):.4f}' for l in [30,50,100,200,300,500,1000]))
for cat,rmin,z0 in [('desi_qso_bin0',60,2.0),('lsst_gal_bin0',40,1.0)]:
    lmax=int(2*np.pi/rmin*cc.chi(z0)*h-0.5); mask=edges[1:]<=lmax
    res0=[];res1=[]
    for b in boxes6:
        dd=data(b,cat); tg,tk=gal_theory(b,cat); T=T_of(b)
        x0,e0,_,_=fit_joint(dd['gg'],dd['kg'],tg,tk,dd['sn'],mask)
        x1,e1,_,_=fit_joint(dd['gg'],dd['kg'],tg,tk*np.sqrt(T),dd['sn'],mask)
        res0.append(x0[1]);res1.append(x1[1])
        print(f'{cat} {b}: s8 no-T={x0[1]:.4f}  with leave-one-out T={x1[1]:.4f} ±{e1[1]:.4f}')
    res0=np.array(res0);res1=np.array(res1); e=e1[1]/np.sqrt(6)
    print(f'{cat} rmin={rmin}: mean s8 no-T={res0.mean():.4f} ({(res0.mean()-s0)/e:+.1f}σ); with T={res1.mean():.4f} ({(res1.mean()-s0)/e:+.1f}σ); σ_mean≈{e:.4f}')
