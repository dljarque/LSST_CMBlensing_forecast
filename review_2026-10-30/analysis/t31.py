from t30 import *
for cat in ['lsst_gal_bin0','desi_qso_bin0']:
    tg,tk=gal_theory('test_box_11',cat); Ds=[data(b,cat) for b in boxes]
    for xx,t in [('gg',tg),('kg',tk)]:
        for a_,b_ in [(30,60),(60,120),(120,240),(240,360)]:
            m=(ells>=a_)&(ells<b_); w=2*ells[m]+1
            vals=np.array([np.sum(w*d[xx][m])/np.sum(w*t[m]) for d in Ds])
            sn=Ds[0]['sn']
            if xx=='gg': knox=np.sqrt(2/np.sum(w))*np.sum(w*(t[m]+sn))/np.sum(w*t[m])
            else: knox=np.sqrt(np.sum(w*(t[m]**2+kk_th[m]*(tg[m]+sn)))/np.sum(w)**2*1)/ (np.sum(w*t[m])/np.sum(w))
            print(f'{cat} {xx} ell[{a_},{b_}) data/theory per box: {np.round(vals,4)} scatter={vals.std(ddof=1):.4f} Knox frac err={knox:.4f}')
