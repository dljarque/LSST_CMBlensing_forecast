from glob import glob
import os

import numpy as np
import healpy as hp
from astropy.io import fits

from utils.rebinning import hybrid_bin_cls


################################################################################
#                                                                              #
#   COMMON FUNCTIONS SECTION                                                   #
#   ------------------------------------------------------------------------   #
#   The functions below handle common tasks such as merging source catalogs,   #
#   and saving results.                                           #
#                                                                              #
################################################################################


def merge_sources_catalogs(sources_prefix, out_dir, filename, zmean, dz, sigma_z):
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

    Z_true = np.concatenate(z_list)

    print(f'Merged {len(files)} catalogs into a single catalog with {len(Z_true)} sources.')

    # Check if all tomographic parameters are provided
    run_tomography = (zmean is not None) and (sigma_z is not None) and (dz is not None)

    if run_tomography:
        # Create the tomographic redshift bin by adding Gaussian noise to true redshifts and filtering out sources outside the bin
        print(f'Running tomographic binning with z_mean={zmean}, sigma_z={sigma_z}, dz={dz}')

        Z_min_bin = zmean - dz / 2
        Z_max_bin = zmean + dz / 2
        Z_photo = np.random.normal(loc=Z_true, scale=sigma_z, size=len(Z_true))  # add Gaussian noise to the redshifts
        mask = (Z_photo > Z_min_bin) & (Z_photo < Z_max_bin)  # filter out sources for the bin

        RA = np.concatenate(ra_list)[mask]
        DEC = np.concatenate(dec_list)[mask]
        Z = Z_true[mask]
    else:
        # If parameters are missing, bypass noise addition and keep the entire distribution
        print('No tomographic binning applied. Using all sources.')
        RA = np.concatenate(ra_list)
        DEC = np.concatenate(dec_list)
        Z = Z_true
        
    # Compute the redshift distribution n(z) for the selected sources in the bin
    Nz, z_edges = np.histogram(Z, bins='fd')  # Binning the redshift distribution with optimal nbins
    z_bin_width = z_edges[1] - z_edges[0]  # Bin width
    z_hist = 0.5 * (z_edges[1:] + z_edges[:-1])  # Bin centers
    print(len(Nz), len(z_hist), z_bin_width)

    # Store sources catalog and bins count
    os.makedirs(out_dir, exist_ok=True)  # Create the output directory if it doesn't exist

    # Dynamically choose what variables to save based on the tomography flag    
    if run_tomography:
        np.savez_compressed(f'{out_dir}/nz_{filename}.npz', z_nz=z_hist, nz=Nz/z_bin_width, z_max=Z.max(), z_min=Z.min(), z_mean=zmean, dz=dz, sigma_z=sigma_z)
    else:
        np.savez_compressed(f'{out_dir}/nz_{filename}.npz', z_nz=z_hist, nz=Nz/z_bin_width, z_max=Z.max(), z_min=Z.min())

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


def compute_cls_ss(sources_prefix, nside, lmax, zmean, dz, sigma_z, out_dir, filename):
    '''
    Computes the sources clustering auto-spectrum, substracting shot noise
    '''

    # Merge sources catalogs
    ra, dec, _ = merge_sources_catalogs(sources_prefix, out_dir, filename, zmean=zmean, dz=dz, sigma_z=sigma_z)

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
    print(f"Mean counts per pixel for the input nside: {mean_counts:.3f}")

    # Compute the overdensity map (delta) for pixels with counts > 0, and set delta=0 for empty pixels to avoid division by zero
    delta = counts / mean_counts - 1.0

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

    # subtract shot noise from the auto-spectrum to get an estimate of the signal
    ell = np.arange(len(cl_ss))  # corresponding multipole moments
    cl_ss_signal = cl_ss - shot_noise
    var_ss  = (2 / ((2*ell + 1)*fsky)) * (cl_ss)**2   # variance on the Cls, assuming Gaussianity

    return ell, cl_ss_signal, np.sqrt(var_ss), shot_noise


def save_rebinned_cls_ss(sources_prefix, nside, lmin, lmax, zmean, dz, sigma_z, nbins, out_dir, filename):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''
    import time
    t0 = time.perf_counter()

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ss, err_ss, shot_noise = compute_cls_ss(sources_prefix, nside, lmax, zmean, dz, sigma_z, out_dir, filename)
    ell_log, cl_ss_log, err_ss_log = hybrid_bin_cls(ell, cl_ss, err_ss, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_gg_{filename}', ells=ell, cls=cl_ss, err=err_ss, nside=nside, shot_noise=shot_noise)
    np.savez(f'{out_dir}/cl_gg_{filename}_log', ells=ell_log, cls=cl_ss_log, err=err_ss_log, nside=nside, shot_noise=shot_noise)
    print(f'Saved cl_gg_{filename}.npz, cl_gg_{filename}_log.npz and nz_{filename}.npz arrays in {out_dir} succesfully.')
    print(f'Time taken: {time.perf_counter() - t0:.2f} seconds')

################################################################################
#                                                                              #
#   CMB LENSING X SOURCES CROSS-SPECTRUM SECTION                               #
#   ------------------------------------------------------------------------   #
#   The functions below handle merging source catalogs, building HEALPix maps, #
#   computing the cmb lensing x sources cross-power spectrum (shot-noise subtraction), 
#   and saving rebinned Cls.                                                   #
#                                                                              #
################################################################################


def compute_cls_ks(kappa_map, sources_prefix, lmax, zmean, dz, sigma_z, out_dir, filename):
    '''
    Computes the cmb lensing x sources cross-spectrum, substracting shot noise
    '''

    # Read the kappa map (nside=512 by colore default) and upscale it to input nside
    m_k = hp.read_map(kappa_map, hdu=1)

    # Get the map resolution
    nside = hp.get_nside(m_k)
    lmax_map = 3 * nside - 1  # maximum multipole for the map, set by the HEALPix resolution

    # Merge sources catalogs
    ra, dec, _ = merge_sources_catalogs(sources_prefix, out_dir, filename, zmean=zmean, dz=dz, sigma_z=sigma_z)

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
    print(f"Mean counts per pixel for the input nside: {mean_counts:.3f}")

    # Compute the overdensity map (delta) for pixels with counts > 0, and set delta=0 for empty pixels to avoid division by zero
    delta = counts / mean_counts - 1.0

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
    nl_kk = 0  # no CMB lensing reconstruction noise here


    # subtract shot noise from the auto-spectrum to get an estimate of the signal
    ell = np.arange(len(cl_ks))  # corresponding multipole moments
    var_ks = (2 / ((2*ell + 1)*fsky)) * ((cl_ks)**2 + (cl_kk + nl_kk)*(cl_ss))

    return ell, cl_ks, np.sqrt(var_ks), nside, shot_noise


def save_rebinned_cls_ks(kappa_map, sources_prefix, lmin, lmax, zmean, dz, sigma_z, nbins, out_dir, filename):
    '''
    Read the sources catalogs from the specified folder, compute the power spectrum, rebin it into logarithmic bins, and save the rebinned Cls and their errors to .npy files.
    '''

    # Read the sources power spectrum and compute the error bars assuming Gaussianity and no CMB lensing noise
    ell, cl_ks, err_ks, nside, shot_noise = compute_cls_ks(kappa_map, sources_prefix, lmax, zmean, dz, sigma_z, out_dir, filename)
    ell_log, cl_ks_log, err_ks_log = hybrid_bin_cls(ell, cl_ks, err_ks, nbins=nbins, ell_min=lmin, ell_max=lmax)

    # Save the rebinned Cls and their errors to .npy files
    np.savez(f'{out_dir}/cl_kg_{filename}', ells=ell, cls=cl_ks, err=err_ks, nside=nside, shot_noise=shot_noise)
    np.savez(f'{out_dir}/cl_kg_{filename}_log', ells=ell_log, cls=cl_ks_log, err=err_ks_log, nside=nside, shot_noise=shot_noise)
    print(f'Saved cl_kg_{filename}.npz, cl_kg_{filename}_log.npz and nz_{filename}.npz arrays in {out_dir} succesfully.')