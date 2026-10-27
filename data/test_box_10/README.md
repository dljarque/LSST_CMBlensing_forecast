# CoLoRe simulation boxes

## Description

Seventh CoLoRe box, with CMB lensing, DESI galaxies and LSST galaxies output, run by Dani. This box mimics the seventh box, same seed, same NSIDE, [pk_linear.dat](pk_linear.dat) extended from 500 points to 1000, but still same pk, changing R_smoothing to 1 Mpc/h. This box was run to cross-check that parameters that shouldn't matter in the param.cfg file for lognormal sims actually don't change the results.

## Products

CMB lensing convergence $\kappa$ values, Galaxy auto-spectrum, and CMB x Galaxy cross-spectrum of the test box from the [CoLoRe simulation](https://github.com/damonge/CoLoRe). 

- File [out_kappa.fits](cls/out_kappa_z000.fits) contains the convergence $\kappa$ values.

- File [cl_kk.npz](cls/cl_kk.npz) contains the CMB lensing autospectrum $C_{\ell}^{\kappa \kappa}$

- File [cl_kk_log.npz](cls/cl_kk_log.npz) contains the log-rebinned CMB lensing autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

- Files [cl_gg_survey_galaxy.npz](cls/cl_gg_desi_qso.npz) contain the galaxy clustering autospectrum $C_{\ell}^{gg}$ of a particular survey galaxy sample.

- Files [cl_gg_survey_galaxy_log.npz](cls/cl_gg_desi_qso.npz) contain the log-rebinned galaxy clustering autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

- Files [cl_kg_survey_galaxy.npz](cls/cl_kg_desi_qso.npz) contain the CMB lensing x Galaxy cross-spectrum $C_{\ell}^{kq}$ of a particular survey galaxy sample.

- Files [cl_kg_survey_galaxy_log.npz](cls/cl_kg_desi_qso_log.npz) contain the log-rebinned CMB lensing x galaxy cross-spectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

## CoLoRe config file

**This box was run using the LINEAR bias model of CoLoRe.**

The file [param.cfg](param.cfg) contains the input parameters to reproduce the colore box. The config file refers to five different files:

- File [pk_linear.dat](pk_linear.dat) contains the linear matter power spectrum $P(k)$ at $z=0$.
- Files [nz_survey_galaxy.txt](../colore_files/nz/nz_desi_qso.txt) contain the redshift distributions per square degree $dn(z)/d\Omega$ of the different galaxy samples used for this box.
- Files [bz_survey_galaxy.txt](../colore_files/bz/bz_desi_qso.txt) contain the bias parameter $b(z)$ evolution with redshift of the different galaxy samples used for this box.
- Files [bz-tunning_survey_galaxy.txt](../colore_files/bz-tunning/bz-tunning_desi_qso.txt) contains the tunning bias parameter $b_{\rm tunning}(z)$ evolution with redshift of the different galaxy samples used for this box.
- Files [Tz_survey_galaxy.txt](../colore_files/Tz/Tz_desi_qso.txt) contain the threshold parameter $T_{\delta}(z)$ evolution with redshift of the different galaxy samples used for this box.

## Cosmological parameters

The file [cosmo_params.json](cosmo_params.json) contains the relevant cosmological parameters matching the box along with the box max redshift `z_max` and source plane redshift `z_out`.
