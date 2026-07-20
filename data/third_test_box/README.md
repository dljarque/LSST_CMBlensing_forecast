# CoLoRe simulation boxes

CMB lensing convergence $\kappa$ values, and QSO catalogs of the test box from the [CoLoRe simulation](https://github.com/damonge/CoLoRe). 

- File [out_kappa.fits](cmb_lensing/out_kappa_z000.fits) contains the convergence $\kappa$ values.
- File [cl_kk.npz](cmb_lensing/cl_kk.npz) contains the CMB lensing autospectrum $C_{\ell}^{\kappa \kappa}$
- File [cl_kk_log.npz](cmb_lensing/cl_kk_log.npz) contains the log-rebinned CMB lensing autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

- Files [out_srcs_s1_i.fits](desi_qso/out_srcs_s1_0.fits) contains the ra, dec coordinates and redshift of the quasars.
- File [cl_qq.npz](desi_qso/cl_qq.npz) contains the QSO clustering autospectrum $C_{\ell}^{qq}$
- File [cl_qq_log.npz](desi_qso/cl_qq_log.npz) contains the log-rebinned QSO clustering autospectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

- File [cl_kq.npz](desi_qso/cl_kq.npz) contains the CMB lensing x QSO cross-spectrum $C_{\ell}^{kq}$.
- File [cl_kq_log.npz](desi_qso/cl_kq_log.npz) contains the log-rebinned CMB lensing x QSO cross-spectrum using the [rebin_cls.py](../../scripts/rebin_cls.py) script.

## CoLoRe config file

**This box was run using the NORMALIZED bias model of CoLoRe.**

The file [param.cfg](param.cfg) contains the input parameters to reproduce the colore box. The config file refers to three different files:
- File [pk_linear.dat](pk_linear.dat) contains the linear matter power spectrum $P(k)$ at $z=0$.
- File [nz_desi_qso.txt](../colore_files/nz/nz_desi_qso.txt) contains the quasars redshift distribution $n(z)$.
- File [bz_desi_qso.txt](../colore_files/bz/bz_desi_qso.txt) contains the quasar bias parameter $b_Q(z)$ evolution with redshift.

## Cosmological parameters

The file [cosmo_params.json](cosmo_params.json) contains the relevant cosmological parameters matching the box along with the box max redshift `z_max` and source plane redshift `z_out`.

## Fitting results

The file [fit_results.txt](fit_results.txt) contains the best-fit values of two parameters implemented in the theory model to agree with CoLoRe data. The fitting is done for the CMB lensing auto, QSO auto, CMB x QSO both for the full QSO sample and for the tomographic bins.
- $k_p$ parameter accounting for the finite cell size in CoLoRe
$$P_{\rm mm, lin}(k) \longrightarrow P_{\rm mm, lin}(k) \times e^{-(k/k_p)^2}$$
- $A_{\rm bias}$ parameter accounting for a rescaling of the input QSO bias to match CoLoRe's output bias. The bias model that was used in CoLoRe to run this box had this problem. Otherwise, theoretically, input and output bias $b(z)$ should be the same.
$$b_{\rm QSO}(z) \longrightarrow  A_{\rm bias} \times b_{\rm QSO}(z)$$

