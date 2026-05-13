import numpy as np
import healpy as hp

def log_bin_cls(ell, cls, nls=None, nbins=20, ell_min=None, ell_max=None):
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


def save_rebinned_cls(kappa_file, nbins, out_dir):
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