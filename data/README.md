# CoLoRe simulation CMB lensing values

CMB lensing convergence $\kappa$ values of the test box from the [CoLoRe simulation](https://github.com/damonge/CoLoRe). 

- File [out_kappa.fits](out_kappa.fits) contains the convergence $\kappa$ values.
- File [cl_kk.npz](cl_kk.npoz) contains the CMB lensing autospectrum $C_{\ell}^{\kappa \kappa}$
- File [cl_kk_log.npz](cl_kk_log.npz) contains the log rebinned CMB lensing autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

## CoLoRe config file

The file [param.cfg](param.cfg) contains the input parameters to reproduce the colore box.

## Cosmological parameters

The file [cosmology.json](cosmology.json) contains the relevant cosmological parameters matching the box along with the box max redshift `z_max` and source plane redshift `z_out`.

