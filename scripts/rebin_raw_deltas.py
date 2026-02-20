import argparse
import numpy as np
from astropy.io import fits
from multiprocessing import Pool


# hard-coded constants (great coding here)
lya_rf=1215.67


def rebin_deltas(hdus, iq, lmin, lmax, dlam):
    info ={}
    header = hdus[iq].header
    info['zq'] = header['Z']
    info['ra_rad'] = header['RA']
    info['dec_rad'] = header['DEC']
    info['targetid'] = header['THING_ID']
    # read forest before rebinning
    data = hdus[iq].data
    lam = np.power(10.0, data['LOGLAM'])
    delta = data['DELTA']
    weight = data['WEIGHT']
    # rebin quantities
    Nlam=int((lmax-lmin)/dlam)
    info['sum_w']=np.zeros(Nlam)
    info['mean_z']=np.zeros(Nlam)
    info['mean_delta']=np.zeros(Nlam)
    for il in range(Nlam):
        ilmin=lmin+il*dlam
        ilmax=ilmin+dlam
        mask=(lam>=ilmin) & (lam<ilmax)
        info['sum_w'][il]=np.sum(weight[mask])
        if info['sum_w'][il]>0:
            info['mean_z'][il]=np.sum(weight[mask]*(lam[mask]/lya_rf-1))/info['sum_w'][il]
            info['mean_delta'][il]=np.sum(weight[mask]*delta[mask])/info['sum_w'][il]
    return info


def rebin_deltas_in_healpix(fname, lmin, lmax, dlam):
    hdus = fits.open(fname)
    Nq = len(hdus)
    dataset = []
    for iq in range(1, Nq):
        if (iq%200==0):
            print('read quasar', iq)
        dataset.append(rebin_deltas(hdus, iq, lmin, lmax, dlam))
    return dataset


def write_rebinned_deltas(fname, dataset):
    zq = np.array([data['zq'] for data in dataset])
    ra_rad = np.array([data['ra_rad'] for data in dataset])
    dec_rad = np.array([data['dec_rad'] for data in dataset])
    targetid = np.array([data['targetid'] for data in dataset])
    Nq=len(dataset)
    Nlam=len(dataset[0]['sum_w'])
    sum_w=np.zeros([Nq, Nlam])
    mean_z=np.zeros([Nq, Nlam])
    mean_delta=np.zeros([Nq, Nlam])
    for iq in range(Nq):
        sum_w[iq] = dataset[iq]['sum_w']
        mean_z[iq] = dataset[iq]['mean_z']
        mean_delta[iq] = dataset[iq]['mean_delta']
    np.savez(fname, zq=zq, ra_rad=ra_rad, dec_rad=dec_rad, targetid=targetid,
             sum_w=sum_w, mean_z=mean_z, mean_delta=mean_delta)


# folder with the raw deltas
delta_dir='/global/cfs/cdirs/desicollab/users/lauracdp/cmb_lensing_box/mock_analysis/raw_master/deltas_lya/Delta'
# folder with rebinned deltas
out_dir = '/global/cfs/cdirs/desicollab/users/font/cmb_lensing_fun/rebinned_deltas'

# define range such that no pixel falls through the cracks
lmin=3599.9
lmax=5600
dlam=100

parser=argparse.ArgumentParser(description="Rebin raw deltas")
parser.add_argument('--imin', action="store", default=0, type=int)
parser.add_argument('--npix', action="store", default=10, type=int)
parser.add_argument('--nproc', action="store", default=128, type=int)
info = parser.parse_args()
print('got imin =', info.imin)
print('got npix =', info.npix)
print('got nproc =', info.nproc)

def rebin_and_write(ipix):
    in_fname=delta_dir+'/delta-{}.fits.gz'.format(ipix)
    print('read pixel', ipix, 'from file', in_fname)
    dataset = rebin_deltas_in_healpix(in_fname, lmin, lmax, dlam)
    out_fname = out_dir + '/rebinned_deltas_{}.npz'.format(ipix)
    print('write info to file', out_fname)
    write_rebinned_deltas(out_fname, dataset)

ipixs = range(info.imin, info.imin + info.npix)
with Pool(info.nproc) as p:
    p.map(rebin_and_write, ipixs)

