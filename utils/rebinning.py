import numpy as np
import warnings


################################################################################
#                                                                              #
#   COMMON FUNCTIONS SECTION                                                   #
#   ------------------------------------------------------------------------   #
#   The functions below handle common tasks such as rebinning power spectra,   #
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


def make_hybrid_bin_matrix(
    ell,
    nbins,
    nbins_low=4,
    ell_switch=20,
    ell_min=None,
    ell_max=None,
):
    """
    Build the linear matrix B such that

        cls_binned = B @ cls

    reproduces the weighted binning of hybrid_bin_cls.
    """

    if ell_min is None:
        ell_min = ell[2]

    if ell_max is None:
        ell_max = ell.max()

    if ell_min >= ell_switch:
        bins = np.logspace(
            np.log10(ell_min),
            np.log10(ell_max),
            nbins + 1
        )

    elif ell_max <= ell_switch:
        bins = np.linspace(
            ell_min,
            ell_max,
            nbins + 1
        )

    else:
        nbins_high = nbins - nbins_low

        bins_low = np.linspace(
            ell_min,
            ell_switch,
            nbins_low + 1
        )

        bins_high = np.logspace(
            np.log10(ell_switch),
            np.log10(ell_max),
            nbins_high + 1
        )

        bins = np.concatenate([
            bins_low[:-1],
            bins_high
        ])

    B = np.zeros((nbins, len(ell)))

    weights = 2 * ell + 1

    for i in range(nbins):

        if i == nbins - 1:
            mask = (ell >= bins[i]) & (ell <= bins[i + 1])
        else:
            mask = (ell >= bins[i]) & (ell < bins[i + 1])

        if np.any(mask):
            B[i, mask] = weights[mask] / np.sum(weights[mask])

    return B