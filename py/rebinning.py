import numpy as np
import healpy as hp
from astropy.io import fits

def log_bin_cls(ell, cls, nls=None, nbins=None, ell_min=None, ell_max=None):
    '''
    Rebin the Cls into logarithmic bins, weighting by 2l+1 to preserve power.
    '''

    # set default values for ell_min and ell_max if not provided
    if ell_min is None:
        ell_min = ell[2]   # mask out the monopole and dipole

    if ell_max is None:
        ell_max = ell.max()

    # set nls to zero if not provided
    if nls is None:
        nls = np.zeros_like(cls)

    # log bins
    bins = np.logspace(
        np.log10(ell_min),
        np.log10(ell_max),
        nbins + 1
    )

    # arrays to hold the binned values
    ell_binned = np.zeros(nbins)
    cls_binned = np.zeros(nbins)
    err_tt_binned = np.zeros(nbins)
    
    # weights for averaging, proportional to the number of modes in each ell bin
    weights = 2*ell + 1

    # variance of the Cls in the bin, assuming Gaussianity and no CMB lensing noise
    var = 2*(cls + nls)**2 / (2*ell + 1)

    for i in range(nbins):

        # create a mask for the current bin
        mask = (ell >= bins[i]) & (ell < bins[i+1])
 
        # if there are any ells in the bin, compute the weighted average and error
        if np.any(mask):

            # apply the weights to the cls and ells in the bin
            w = weights[mask]

            # compute the weighted average of the cls and ells in the bin, and the rebinned variance
            cls_binned[i] = np.sum(w * cls[mask]) / np.sum(w)  # weighted average of the cls in the bin
            ell_binned[i] = np.sum(w * ell[mask]) / np.sum(w)  # weighted average of the ells in the bin
            var_binned = np.sum((w**2) * var[mask]) / (np.sum(w)**2)  # rebinned variance of the Cls, assuming Gaussianity and no CMB lensing noise
            err_tt_binned[i] = np.sqrt(var_binned)  

        # if there are no ells in the bin, set the binned values to nan
        else:
            cls_binned[i] = np.nan
            ell_binned[i] = np.nan
            err_tt_binned[i] = np.nan
        
    # return the binned ells, Cls, and errors
    return ell_binned, cls_binned, err_tt_binned


def save_rebinned_cls_kk(kappa_file, nbins, out_dir):
    '''
    Read the kappa map from the specified CoLoRe box, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the kappa map and compute the power spectrum
    m_k = hp.read_map(kappa_file, hdu=1)

    # Compute the power spectrum of the kappa map
    cl_kk = hp.anafast(m_k)  # get the power spectrum Cls
    ell = np.arange(len(cl_kk))  # corresponding multipole moments
    err_tt = np.sqrt(2*(cl_kk)**2 / (2*ell + 1))  # error on the Cls, assuming Gaussianity and no CMB lensing noise
    ell_colore_log, cl_kk_log, err_tt_log = log_bin_cls(ell, cl_kk, nbins=nbins)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_kk', ells=ell, cls=cl_kk, err=err_tt)
    np.savez(f'{out_dir}/cl_kk_log', ells=ell_colore_log, cls=cl_kk_log, err=err_tt_log)


############################################################################################################################


def merge_quasar_catalog(quasar_folder):
    '''
    Merge the quasar catalogs from the different HEALPix tiles into a single catalog of RA and DEC.
    The catalogs are assumed to be named out_srcs_s1_i.fits, where i is the tile number.
    '''

    ra_list = []
    dec_list = []
    
    i = 0
    # Loop over the HEALPix tiles until we run out of files
    while True:
        try:
            cat = fits.getdata(f'{quasar_folder}out_srcs_s1_{i}.fits', ext=1)
        except FileNotFoundError:
            break
        ra_list.append(cat['RA'])
        dec_list.append(cat['DEC'])
        i += 1
        ra_list.append(cat['RA'])
        dec_list.append(cat['DEC'])

    RA = np.concatenate(ra_list)
    DEC = np.concatenate(dec_list)

    print(f'Merged {i} catalogs from {quasar_folder} into a single catalog with {len(RA)} quasars.')

    return RA, DEC


def compute_cls_qq(quasars_folder, nside=256):

    ra, dec = merge_quasar_catalog(quasars_folder)

    # Choose resolution
    nside = nside  # resoluition parameter for HEALPix map, nside=256 corresponds to ~14 arcmin pixels

    # Convert coordinates
    theta = np.radians(90.0 - dec)   # colatitude
    phi   = np.radians(ra)  # longitude

    # Pixel indices
    pix = hp.ang2pix(nside, theta, phi)

    # Count quasars per pixel
    npix = hp.nside2npix(nside)
    counts = np.bincount(pix, minlength=npix)
    
    # Create a mask for pixels with counts > 0
    mask = counts > 0
    mean_counts = np.mean(counts)
    print(f"Mean counts per pixel: {mean_counts:.2f}")

    # Check for empty pixels and warn if there are any
    if np.sum(mask) == 0:
        raise ValueError("No pixels with counts > 0. Check the quasar catalog and HEALPix parameters.")
    if False in mask:
        nfalse = np.sum(~mask)
        print(f"Warning: {nfalse} pixels have 0 counts ({nfalse/npix*100:.2f}%). Check the nside resolution.")

    # Compute the overdensity map (delta) for pixels with counts > 0, and set delta=0 for empty pixels to avoid division by zero
    delta = np.zeros_like(counts, dtype=float)
    delta[mask] = (counts[mask] - mean_counts) / mean_counts
    lmax_map = hp.npix2nside(npix) * 3 - 1  # maximum multipole for the map, set by the HEALPix resolution
    cl = hp.anafast(delta, lmax=lmax_map)  # compute the auto power spectrum of the overdensity map

    # compute the effective sky fraction and number density for shot noise calculation
    fsky = 1.0 * np.sum(mask) / npix
    nbar_sr = len(ra) / (4*np.pi*fsky)
    shot_noise = 1.0 / nbar_sr

    # subtract shot noise from the auto-spectrum and deconvolve the pixel window function to get an estimate of the underlying quasar clustering power spectrum
    cl_signal = cl - shot_noise
    pixwin = hp.pixwin(nside, lmax=lmax_map)  # pixel window function for HEALPix map at given nside and lmax NSIDE*3 - 1
    cl_signal_deconvolved = cl_signal / 1

    return cl_signal_deconvolved, cl, fsky, shot_noise


def save_rebinned_cls_qq(quasars_folder, nbins, out_dir):
    '''
    Read the quasar maps from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the quasar power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    cl_qq, cl, fsky, shot_noise = compute_cls_qq(quasars_folder)
    ell = np.arange(len(cl_qq))  # corresponding multipole moments
    err_tt = np.sqrt(2/ ((2*ell + 1)*fsky) ) + (cl + shot_noise)**2   # error on the Cls, assuming Gaussianity
    ell_colore_log, cl_qq_log, err_tt_log = log_bin_cls(ell, cl_qq, nbins=nbins)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_qq', ells=ell, cls=cl_qq, err=err_tt)
    np.savez(f'{out_dir}/cl_qq_log', ells=ell_colore_log, cls=cl_qq_log, err=err_tt_log)