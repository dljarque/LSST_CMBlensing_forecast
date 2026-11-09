from pipeline import *
for box in ['test_box_7','ensemble_1','test_box_11']:
    B=load_box(box,'desi_qso_bin0'); cd=B['cd']; cc=CoLoReCosmo(cd,D+f'{box}/pk_linear.dat'); h=cd['h']; L,dx=colore_grid(cc)
    zz=np.linspace(1e-4,3.8,400); e=np.arange(2,1600); _,_,W,_=kernels(cc,zz,cc.chi(cd['z_out'])); d=B['kk']['cl'][2:1600]
    res={}
    for nm,R2 in [('R2+dx2/4',4+dx**2/4),('R2+dx2/3',4+dx**2/3)]:
        t=limber(cc,e,zz,W,W,smooth_kp=h/np.sqrt(R2)); row=[]
        for a,b in [(20,100),(100,200),(200,300),(300,500),(500,700),(700,1000),(1000,1500)]:
            m=(e>=a)&(e<b); w=2*e[m]+1; row.append('%.4f'%(np.sum(w*d[m])/np.sum(w*t[m])))
        print(box,nm,' '.join(row))
