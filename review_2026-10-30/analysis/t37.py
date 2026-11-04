from t30 import *
boxes6=['test_box_11','ensemble_1','ensemble_2','ensemble_3','ensemble_4','ensemble_5']
cat='desi_qso_bin0'; tg,tk=gal_theory('test_box_11',cat); chi_h=cc.chi(2.0)*h
print('ens5 Ngal=%.4e'%(4*np.pi/data('ensemble_5',cat)['sn']), 'ens4 Ngal=%.4e'%(4*np.pi/data('ensemble_4',cat)['sn']))
for rmin in [60,50,40]:
    lmax=int(2*np.pi/rmin*chi_h-0.5); mask=edges[1:]<=lmax; Ds=[]; s8s=[]; es=[]
    for b in boxes6:
        dd=data(b,cat); Ds.append(dd); tgb,tkb=gal_theory(b,cat)
        x,err,c2,nd=fit_joint(dd['gg'],dd['kg'],tgb,tkb,dd['sn'],mask); s8s.append(x[1]); es.append(err[1])
        print(f'  {b:12s} rmin={rmin}: s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd}')
    gm=np.mean([d['gg'] for d in Ds],axis=0); km=np.mean([d['kg'] for d in Ds],axis=0); snm=np.mean([d['sn'] for d in Ds])
    x,err,c2,nd=fit_joint(gm,km,tg,tk,snm,mask,N=6)
    print(f'QSO rmin={rmin} lmax={lmax}: scatter {np.std(s8s,ddof=1):.4f} vs Knox {np.mean(es):.4f}; AVERAGED(6): s8={x[1]:.4f}±{err[1]:.4f} ({(x[1]-s0)/err[1]:+.1f}σ) A={x[0]:.3f} A_sn={x[2]:.3f} χ2={c2:.0f}/{nd}')
