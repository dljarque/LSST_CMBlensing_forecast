import numpy as np
import healpy as hp
from astropy.io import fits

def hybrid_bin_cls(
    ell,
    cls,
    err=None,
    nbins_low=4,
    nbins=None,
    ell_switch=20,
    ell_min=None,
    ell_max=None
):
    '''
    Rebin the Cls using linear bins at low ell and logarithmic bins at high ell,
    weighting by 2l+1 to preserve power.
    '''

    # set default values for ell_min and ell_max if not provided
    if ell_min is None:
        ell_min = ell[2]   # mask out the monopole and dipole

    if ell_max is None:
        ell_max = ell.max()

    # set default variance to 0 if not provided
    err = np.zeros_like(cls) if err is None else err


    # make nbins consistent with nbins before ell_switch (nbins_low) and after ell_switch (nbins_high)
    nbins_high = nbins - nbins_low

    # linear bins at low ell
    bins_low = np.linspace(
        ell_min,
        ell_switch,
        nbins_low + 1
    )

    # logarithmic bins at high ell
    bins_high = np.logspace(
        np.log10(ell_switch),
        np.log10(ell_max),
        nbins_high + 1
    )

    # merge low and high ell bins, avoiding duplicate ell_switch edge
    bins = np.concatenate([
        bins_low[:-1],
        bins_high
    ])

    # total number of bins
    nbins = len(bins) - 1

    # arrays to hold the binned values
    ell_binned = np.zeros(nbins)
    cls_binned = np.zeros(nbins)
    err_binned = np.zeros(nbins)

    # weights for averaging, proportional to the number of modes in each ell bin
    weights = 2*ell + 1

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
            var_binned = np.sum((w**2) * err[mask]**2) / (np.sum(w)**2)  # rebinned variance of the Cls
            err_binned[i] = np.sqrt(var_binned)  # rebinned error of the Cls

        # if there are no ells in the bin, set the binned values to nan
        else:
            cls_binned[i] = np.nan
            ell_binned[i] = np.nan
            err_binned[i] = np.nan

    # return the binned ells, Cls, and errors
    return ell_binned, cls_binned, err_binned


################################################################################
#                                                                              #
#   CMB LENSING AUTO-SPECTRUM SECTION                                          #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps,   #
#   computing the CMB Lensing auto-power spectrum (shot-noise subtraction), and  #
#   saving rebinned Cls.                                                       #
#                                                                              #
################################################################################


def compute_cls_kk(kappa_file):
    '''
    Compute the CMB lensing auto-spectrum
    '''
    # Read the kappa map
    m_k = hp.read_map(kappa_file, hdu=1)
    fsky = 1  # full sky for the moment

    # Compute the power spectrum of the kappa map
    print('Computing the auto-spectrum of the lensing map...')
    cl_kk = hp.anafast(m_k)  # get the power spectrum Cls
    ell = np.arange(len(cl_kk))  # corresponding multipole moments
    var_kk = (2 / ((2*ell + 1)*fsky)) * (cl_kk)**2   # variance on the Cls, assuming Gaussianity and no CMB lensing noise

    return ell, cl_kk, np.sqrt(var_kk)


def save_rebinned_cls_kk(kappa_file, lmin, lmax, nbins, out_dir):
    '''
    Read the kappa map from the specified CoLoRe box, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the kappa map and compute the power spectrum
    ell, cl_kk, err_kk = compute_cls_kk(kappa_file)
    ell_log, cl_kk_log, err_kk_log = hybrid_bin_cls(ell, cl_kk, err_kk, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_kk', ells=ell, cls=cl_kk, err=err_kk)
    np.savez(f'{out_dir}/cl_kk_log', ells=ell_log, cls=cl_kk_log, err=err_kk_log)

    print(f'Saved cl_kk.npz and cl_kk_log.npz arrays in {out_dir} succesfully.')


################################################################################
#                                                                              #
#   SOURCES AUTO-SPECTRUM SECTION                                              #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps, #
#   computing the sources auto-power spectrum (shot-noise subtraction), and    #
#   saving rebinned Cls.                                                       #
#                                                                              #
################################################################################


def merge_sources_catalogs(sources_folder):
    '''
    Merge the sources catalogs from different files into a single catalog of RA, DEC and Z.
    The catalogs are assumed to be named out_srcs_s1_i.fits, where i is the tile number.
    '''

    ra_list = []
    dec_list = []
    z_list = []
    
    i = 0
    # Loop over the HEALPix tiles until we run out of files
    while True:
        try:
            cat = fits.getdata(f'{sources_folder}out_srcs_s1_{i}.fits', ext=1)
        except FileNotFoundError:
            break
        ra_list.append(cat['RA'])
        dec_list.append(cat['DEC'])
        z_list.append(cat['Z_COSMO'])
        i += 1

    RA = np.concatenate(ra_list)
    DEC = np.concatenate(dec_list)
    Z = np.concatenate(z_list)

    print(f'Merged {i} catalogs from {sources_folder} into a single catalog with {len(RA)} quasars.')

    return RA, DEC, Z


def compute_cls_ss(sources_folder, nside, lmax):
    '''
    Computes the sources clustering auto-spectrum, substracting shot noise
    '''

    # Merge sources catalogs
    ra, dec = merge_sources_catalogs(sources_folder)[:2]

    # Choose resolution
    nside = nside  # resolution parameter for HEALPix map, nside=4096 corresponds to ~1 arcmin pixels

    # Convert coordinates
    theta = np.radians(90.0 - dec)   # colatitude
    phi   = np.radians(ra)  # longitude

    # Pixel indices
    pix = hp.ang2pix(nside, theta, phi)

    # Count sources per pixel
    npix = hp.nside2npix(nside)
    counts = np.bincount(pix, minlength=npix)
    
    # Mean counts per pixel for pixel
    mean_counts = np.mean(counts)
    print(f"Mean counts per pixel for the input nside: {mean_counts:.2f}")

    # Compute the overdensity map (delta) for pixels with counts > 0, and set delta=0 for empty pixels to avoid division by zero
    delta = (counts - mean_counts) / mean_counts

    lmax_map = nside * 3 - 1  # maximum multipole for the map, set by the HEALPix resolution
    if lmax > lmax_map:
        raise ValueError('lmax is too big to compute the Cls. Maximum multipole allowed: lmax = (3 * NSIDE) - 1')
    elif lmax is None:
        print('Computing the auto-spectrum of the overdensity map...')
        cl_ss = hp.anafast(delta, lmax=lmax_map)  # compute the auto power spectrum of the overdensity map
    else:
        print('Computing the auto-spectrum of the overdensity map...')
        cl_ss = hp.anafast(delta, lmax=lmax)  # compute the auto power spectrum of the overdensity map

    # compute the effective sky fraction and number density for shot noise calculation
    fsky = 1.0  # full sky so far, otherwise 1.0 * np.sum(mask) / npix
    nbar_sr = np.sum(counts) / (4*np.pi*fsky)
    shot_noise = 1.0 / nbar_sr
    nl_ss = np.ones(len(cl_ss)) * shot_noise

    # subtract shot noise from the auto-spectrum to get an estimate of the signal
    ell = np.arange(len(cl_ss))  # corresponding multipole moments
    cl_ss_signal = cl_ss - nl_ss
    var_ss  = (2 / ((2*ell + 1)*fsky)) * (cl_ss_signal + nl_ss)**2   # variance on the Cls, assuming Gaussianity

    return ell, cl_ss_signal, np.sqrt(var_ss)


def save_rebinned_cls_ss(sources, sources_folder, nside, lmin, lmax, nbins, out_dir):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ss, err_ss = compute_cls_ss(sources_folder, nside, lmax)
    ell_log, cl_ss_log, err_ss_log = hybrid_bin_cls(ell, cl_ss, err_ss, nbins=nbins, ell_min=lmin, ell_max=lmax)

    if sources == 'galaxies':
        xx = 'gg'
    elif sources == 'quasars':
        xx = 'qq'
    else:
        xx = 'xx'

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_{xx}', ells=ell, cls=cl_ss, err=err_ss)
    np.savez(f'{out_dir}/cl_{xx}_log', ells=ell_log, cls=cl_ss_log, err=err_ss_log)
    print(f'Saved cl_{xx}.npz and cl_{xx}_log.npz arrays in {out_dir} succesfully.')
    

################################################################################
#                                                                              #
#   CMB LENSING X SOURCES CROSS-SPECTRUM SECTION                               #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps, #
#   computing the cmb lensing x sources cross-power spectrum (shot-noise subtraction), 
#   and saving rebinned Cls.                                                   #
#                                                                              #
################################################################################



def compute_cls_ks(kappa_file, sources_folder, nside, lmax):
    '''
    Computes the cmb lensing x sources cross-spectrum, substracting shot noise
    '''

    # Read the kappa map (nside=512 by colore default) and upscale it to input nside
    m_k = hp.read_map(kappa_file, hdu=1)
    m_k = hp.ud_grade(m_k, nside)

    # Merge sources catalogs
    ra, dec = merge_sources_catalogs(sources_folder)[:2]

    # Convert coordinates
    theta = np.radians(90.0 - dec)   # colatitude
    phi   = np.radians(ra)  # longitude

    # Pixel indices
    pix = hp.ang2pix(nside, theta, phi)

    # Count sources per pixel
    npix = hp.nside2npix(nside)
    counts = np.bincount(pix, minlength=npix)
    
    # Mean counts per pixel for pixel
    mean_counts = np.mean(counts)
    print(f"Mean counts per pixel for the input nside: {mean_counts:.2f}")

    # Compute the overdensity map (delta) for pixels with counts > 0, and set delta=0 for empty pixels to avoid division by zero
    delta = (counts - mean_counts) / mean_counts

    lmax_map = nside * 3 - 1  # maximum multipole for the map, set by the HEALPix resolution
    if lmax > lmax_map:
        raise ValueError('lmax is too big to compute the Cls. Maximum multipole allowed: lmax = (3 * NSIDE) - 1')
    elif lmax is None:
        print('Computing the cross-spectrum of CMB lensing x sources overdensity...')
        cl_ks = hp.anafast(m_k, delta, lmax=lmax_map)  # compute the cross power spectrum
        cl_ss = hp.anafast(delta, lmax=lmax_map)  # compute the ss auto power spectrum
        cl_kk = hp.anafast(m_k, lmax=lmax_map)  # compute the kk auto power spectrum
    else:
        print('Computing the cross-spectrum of CMB lensing x sources overdensity...')
        cl_ks = hp.anafast(m_k, delta, lmax=lmax)  # compute the auto power spectrum of the overdensity map
        cl_ss = hp.anafast(delta, lmax=lmax)  # compute the auto power spectrum of the overdensity map
        cl_kk = hp.anafast(m_k, lmax=lmax)  # compute the kk auto power spectrum


    # compute the effective sky fraction and number density for shot noise calculation
    fsky = 1  # full sky for the moment
    nbar_sr = np.sum(counts) / (4*np.pi*fsky)
    shot_noise = 1.0 / nbar_sr
    nl_ss = np.ones(len(cl_ss)) * shot_noise
    nl_kk = 0  # no CMB lensing reconstruction noise here


    # subtract shot noise from the auto-spectrum to get an estimate of the signal
    ell = np.arange(len(cl_ks))  # corresponding multipole moments
    cl_ss_signal = cl_ss - nl_ss
    var_ks = (2 / ((2*ell + 1)*fsky)) * ((cl_ks)**2 + (cl_kk + nl_kk)*(cl_ss_signal + nl_ss))

    return ell, cl_ks, np.sqrt(var_ks)


def save_rebinned_cls_ks(kappa_file, sources, sources_folder, nside, lmin, lmax, nbins, out_dir):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ks, err_ks = compute_cls_ks(kappa_file, sources_folder, nside, lmax)
    ell_log, cl_ks_log, err_ks_log = hybrid_bin_cls(ell, cl_ks, err_ks, nbins=nbins, ell_min=lmin, ell_max=lmax)

    if sources == 'galaxies':
        xx = 'kg'
    elif sources == 'quasars':
        xx = 'kq'
    else:
        xx = 'kx'

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_{xx}', ells=ell, cls=cl_ks, err=err_ks)
    np.savez(f'{out_dir}/cl_{xx}_log', ells=ell_log, cls=cl_ks_log, err=err_ks_log)
    print(f'Saved cl_{xx}.npz and cl_{xx}_log.npz arrays in {out_dir} succesfully.')