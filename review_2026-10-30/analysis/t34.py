from t30 import *
cat='desi_qso_bin0'; tg,tk=gal_theory('test_box_11',cat)
for b in ['test_box_11','ensemble_3','ensemble_4','ensemble_5']:
    dd=data(b,cat); f=np.load(D+f'{b}/cls/nz_{cat}.npz'); z_,n_=f['z_nz'],f['nz']
    print(b,'sn=%.4e'%dd['sn'],'Ngal=%.0f'%(4*np.pi/dd['sn']),'nz mean=%.4f std=%.4f range %.2f-%.2f nbins %d'%(np.sum(z_*n_)/np.sum(n_),np.sqrt(np.sum((z_-np.sum(z_*n_)/np.sum(n_))**2*n_)/np.sum(n_)),z_.min(),z_.max(),len(z_)), {k:(float(f[k]) if f[k].ndim==0 else '') for k in f.files if f[k].ndim==0})
    for xx,t in [('gg',tg),('kg',tk)]:
        row=[]
        for a_,b_ in [(30,60),(60,120),(120,240),(240,360),(360,560)]:
            m=(ells>=a_)&(ells<b_); w=2*ells[m]+1; row.append(f'{np.sum(w*dd[xx][m])/np.sum(w*t[m]):.3f}')
        print('   ',xx,' '.join(row))
