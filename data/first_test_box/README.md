# CoLoRe simulation boxes

## Description

First CoLoRe box with CMB lensing output, run by Laura. The source plane redshift `z_out` is equal to the box maximum redshift `z_max`. This issue is solved in the second box.

## Products

CMB lensing convergence $\kappa$ values of the test box from the [CoLoRe simulation](https://github.com/damonge/CoLoRe). 

- File [out_kappa.fits](cmb_lensing/out_kappa.fits) contains the convergence $\kappa$ values.
- File [cl_kk.npz](cmb_lensing/cl_kk.npz) contains the CMB lensing autospectrum $C_{\ell}^{\kappa \kappa}$
- File [cl_kk_log.npz](cmb_lensing/cl_kk_log.npz) contains the log-rebinned CMB lensing autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

## CoLoRe config file

**This box was run using the NORMALIZED bias model of CoLoRe.**

The file [param.cfg](param.cfg) contains the input parameters to reproduce the colore box. The config file refers to three different files:
- File [pk_linear.dat](pk_linear.dat) contains the linear matter power spectrum P(k) at z=0.
- File [nz_desi_qso.txt](../colore_files/nz/nz_desi_qso.txt) contains the quasars redshift distribution n(z).
- File [bz_desi_qso.txt](../colore_files/bz/bz_desi_qso.txt) contains the quasar bias parameter b_Q(z) evolution with redshift.

## Cosmological parameters

The file [cosmo_params.json](cosmo_params.json) contains the relevant cosmological parameters matching the box along with the box max redshift `z_max` and source plane redshift `z_out`.

