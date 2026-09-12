from glob import glob
import os

import numpy as np
import healpy as hp
from astropy.io import fits
import warnings


################################################################################
#                                                                              #
#   COMMON FUNCTIONS SECTION                                                   #
#   ------------------------------------------------------------------------   #
#   The functions below handle common tasks such as merging source catalogs,   #
#   rebinning power spectra, and saving results.                               #
#                                                                              #
################################################################################


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

    # check that nbins is provided
    if nbins is None:
            raise ValueError("nbins must be specified")

    # set default errors to zero if not provided
    err = np.zeros_like(cls) if err is None else err


    # if ell_min is above ell_switch, skip the linear part entirely
    if ell_min >= ell_switch:

        warnings.warn(
            f"ell_min ({ell_min}) >= ell_switch ({ell_switch}). "
            "Skipping linear binning and using logarithmic bins only."
        )

        bins = np.logspace(
            np.log10(ell_min),
            np.log10(ell_max),
            nbins + 1
        )

    # if ell_max is below ell_switch, skip the logarithmic part entirely
    elif ell_max <= ell_switch:

        warnings.warn(
            f"ell_max ({ell_max}) <= ell_switch ({ell_switch}). "
            "Skipping logarithmic binning and using linear bins only."
        )

        bins = np.linspace(
            ell_min,
            ell_max,
            nbins + 1
        )

    # otherwise use the standard hybrid binning
    else:

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

    # arrays to hold the binned values
    ell_binned = np.full(nbins, np.nan)
    cls_binned = np.full(nbins, np.nan)
    err_binned = np.full(nbins, np.nan)

    # weights for averaging, proportional to the number of modes in each ell bin
    weights = 2*ell + 1

    for i in range(nbins):

        # define the mask for the ells in the current bin, 
        # including the upper edge for the last bin to ensure we include all ells up to ell_max
        if i == nbins - 1:
            mask = (ell >= bins[i]) & (ell <= bins[i+1])
        else:
            mask = (ell >= bins[i]) & (ell < bins[i+1])

        # if there are any ells in the bin, compute the weighted average and error
        if np.any(mask):

            # apply the weights to the cls and ells in the bin
            w = weights[mask]

            # compute the weighted average of the cls and ells in the bin, and the rebinned variance
            ell_binned[i] = np.sum(w * ell[mask]) / np.sum(w)  # weighted average of the ells in the bin
            cls_binned[i] = np.sum(w * cls[mask]) / np.sum(w)  # weighted average of the cls in the bin
            var_binned = np.sum((w**2) * err[mask]**2) / (np.sum(w)**2)  # rebinned variance of the Cls
            err_binned[i] = np.sqrt(var_binned)  # rebinned error of the Cls

    # return the binned ells, Cls, and errors
    return ell_binned, cls_binned, err_binned


def merge_sources_catalogs(sources_prefix, out_dir, filename):
    '''
    Merge and store all CoLoRe source catalogs sharing the same prefix into a single
    catalog of RA, DEC and Z.

    Example
    -------
    sources_prefix = 'sources_folder/out_srcs_s1'

    matches:
        out_srcs_s1_0.fits
        out_srcs_s1_1.fits
        ...
    '''

    # Get all source catalogs matching the prefix
    files = sorted(glob(f'{sources_prefix}_*.fits'))

    if not files:
        raise FileNotFoundError(f'No source catalogs found matching {sources_prefix}_*.fits')

    ra_list, dec_list, z_list = [], [], []

    # Loop over the files and read the RA, DEC, and Z_COSMO columns
    for file in files:
        cat = fits.getdata(file, ext=1)
        ra_list.append(cat['RA'])
        dec_list.append(cat['DEC'])
        z_list.append(cat['Z_COSMO'])

    RA = np.concatenate(ra_list)
    DEC = np.concatenate(dec_list)
    Z = np.concatenate(z_list)

    print(f'Merged {len(files)} catalogs into a single catalog with {len(RA)} sources.')

    # Store sources catalog and bins count
    os.makedirs(out_dir, exist_ok=True)  # Create the output directory if it doesn't exist
    nz, z_edges = np.histogram(Z, bins=200)  # Binning the redshift distribution into 200 bins
    np.savez_compressed(f'{out_dir}/nz_hist_{filename}.npz', z_edges=z_edges, nz=nz)

    return RA, DEC, Z


################################################################################
#                                                                              #
#   CMB LENSING AUTO-SPECTRUM SECTION                                          #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps,   #
#   computing the CMB Lensing auto-power spectrum (shot-noise subtraction), and  #
#   saving rebinned Cls.                                                       #
#                                                                              #
################################################################################


def compute_cls_kk(kappa_map, noise_spectrum_file=None, return_raw_cls=False):
    '''
    Compute the CMB lensing auto-spectrum
    '''
    # Read the kappa map
    m_k = hp.read_map(kappa_map, hdu=1)
    fsky = 1  # full sky for the moment

    # Get the map resolution
    nside = hp.get_nside(m_k)
    lmax_map = 3 * nside - 1

    # Add a realization of the reconstruction noise if provided
    if noise_spectrum_file is not None:

        _, nl_kk = np.loadtxt(noise_spectrum_file, unpack=True)

        # Match the spectrum length to the map resolution
        if len(nl_kk) < lmax_map + 1:
            nl_kk_match = np.zeros(lmax_map + 1)
            nl_kk_match[:len(nl_kk)] = nl_kk
        else:
            nl_kk_match = nl_kk[:lmax_map + 1]

        # Generate a Gaussian noise realization and add it to the map
        noise_map = hp.synfast(nl_kk_match, nside=nside, lmax=lmax_map)
        m_k = m_k + noise_map

    else:

        nl_kk_match = np.zeros(lmax_map + 1)

    # Compute the power spectrum of the kappa map
    print('Computing the auto-spectrum of the lensing map...')
    cl_kk = hp.anafast(m_k, lmax=lmax_map)  # get the power spectrum Cls
    cl_kk_signal = cl_kk - nl_kk_match
    ell = np.arange(len(cl_kk))  # corresponding multipole moments
    var_kk = (2 / ((2*ell + 1)*fsky)) * (cl_kk_signal + nl_kk_match)**2   # variance on the Cls, assuming Gaussianity and no CMB lensing noise

    if return_raw_cls:
        return ell, cl_kk_signal, np.sqrt(var_kk), cl_kk, nside
    else:
        return ell, cl_kk_signal, np.sqrt(var_kk), nside


def save_rebinned_cls_kk(kappa_map, lmin, lmax, nbins, out_dir):
    '''
    Read the kappa map from the specified CoLoRe box, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the kappa map and compute the power spectrum
    ell, cl_kk, err_kk, nside = compute_cls_kk(kappa_map)
    if lmax is None:
        lmax = 3 * nside - 1  # maximum multipole for the map, set by the HEALPix resolution
    ell_log, cl_kk_log, err_kk_log = hybrid_bin_cls(ell, cl_kk, err_kk, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Create the output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_kk', ells=ell, cls=cl_kk, err=err_kk, nside=nside)
    np.savez(f'{out_dir}/cl_kk_log', ells=ell_log, cls=cl_kk_log, err=err_kk_log, nside=nside)

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



def compute_cls_ss(sources_prefix, nside, lmax, zmin_bin, zmax_bin, out_dir, filename):
    '''
    Computes the sources clustering auto-spectrum, substracting shot noise
    '''

    # Merge sources catalogs
    ra, dec, z = merge_sources_catalogs(sources_prefix, out_dir, filename)

    # Create a bin mask to filter out sources outside the bin
    mask = (z > zmin_bin) & (z < zmax_bin)
    ra = ra[mask]
    dec = dec[mask]
    z = z[mask]

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

    return ell, cl_ss_signal, np.sqrt(var_ss), shot_noise


def save_rebinned_cls_ss(sources_prefix, nside, lmin, lmax, zmin_bin, zmax_bin, nbins, out_dir, filename):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ss, err_ss, shot_noise = compute_cls_ss(sources_prefix, nside, lmax, zmin_bin, zmax_bin, out_dir, filename)
    ell_log, cl_ss_log, err_ss_log = hybrid_bin_cls(ell, cl_ss, err_ss, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_gg_{filename}', ells=ell, cls=cl_ss, err=err_ss, nside=nside, zmin=zmin_bin, zmax=zmax_bin, shot_noise=shot_noise)
    np.savez(f'{out_dir}/cl_gg_{filename}_log', ells=ell_log, cls=cl_ss_log, err=err_ss_log, nside=nside, zmin=zmin_bin, zmax=zmax_bin, shot_noise=shot_noise)
    print(f'Saved cl_gg_{filename}.npz and cl_gg_{filename}_log.npz arrays in {out_dir} succesfully.')
    

################################################################################
#                                                                              #
#   CMB LENSING X SOURCES CROSS-SPECTRUM SECTION                               #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps, #
#   computing the cmb lensing x sources cross-power spectrum (shot-noise subtraction), 
#   and saving rebinned Cls.                                                   #
#                                                                              #
################################################################################


def compute_cls_ks(kappa_map, sources_prefix, lmax, zmin_bin, zmax_bin, out_dir, filename):
    '''
    Computes the cmb lensing x sources cross-spectrum, substracting shot noise
    '''

    # Read the kappa map (nside=512 by colore default) and upscale it to input nside
    m_k = hp.read_map(kappa_map, hdu=1)

    # Get the map resolution
    nside = hp.get_nside(m_k)
    lmax_map = 3 * nside - 1  # maximum multipole for the map, set by the HEALPix resolution

    # Merge sources catalogs
    ra, dec, z = merge_sources_catalogs(sources_prefix, out_dir, filename)

    # Create a bin mask to filter out sources outside the bin
    mask = (z > zmin_bin) & (z <= zmax_bin)
    ra = ra[mask]
    dec = dec[mask]
    z = z[mask]

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

    if lmax is None:
        print('Computing the cross-spectrum of CMB lensing x sources overdensity...')
        cl_ks = hp.anafast(m_k, delta, lmax=lmax_map)  # compute the cross power spectrum
        cl_ss = hp.anafast(delta, lmax=lmax_map)  # compute the ss auto power spectrum
        cl_kk = hp.anafast(m_k, lmax=lmax_map)  # compute the kk auto power spectrum
    else:
        if lmax > lmax_map:
            raise ValueError('lmax is too big to compute the Cls. Maximum multipole allowed: lmax = (3 * NSIDE) - 1')
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

    return ell, cl_ks, np.sqrt(var_ks), nside, shot_noise


def save_rebinned_cls_ks(kappa_map, sources_prefix, lmin, lmax, zmin_bin, zmax_bin, nbins, out_dir, filename):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ks, err_ks, nside, shot_noise = compute_cls_ks(kappa_map, sources_prefix, lmax, zmin_bin, zmax_bin, out_dir, filename)
    ell_log, cl_ks_log, err_ks_log = hybrid_bin_cls(ell, cl_ks, err_ks, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_kg_{filename}', ells=ell, cls=cl_ks, err=err_ks, nside=nside, zmin=zmin_bin, zmax=zmax_bin, shot_noise=shot_noise)
    np.savez(f'{out_dir}/cl_kg_{filename}_log', ells=ell_log, cls=cl_ks_log, err=err_ks_log, nside=nside, zmin=zmin_bin, zmax=zmax_bin, shot_noise=shot_noise)
    print(f'Saved cl_kg_{filename}.npz and cl_kg_{filename}_log.npz arrays in {out_dir} succesfully.')