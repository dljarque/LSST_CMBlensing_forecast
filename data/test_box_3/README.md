# CoLoRe box `test_box_3`

## Description

Third CoLoRe box, with CMB lensing and QSO output, run by Dani.

**Bias model:** this box was run with the **NORMALIZED** bias model of CoLoRe.

Source populations in `results/`: `srcs1` = DESI QSO (`out_srcs_s1_*`).

## Contents

| Path | Content |
| --- | --- |
| [`param.cfg`](param.cfg) | CoLoRe configuration used for this run |
| [`cosmo_params.json`](cosmo_params.json) | cosmology of the box, `z_max` and the source-plane redshift `z_out` |
| [`pk_linear.dat`](pk_linear.dat) | input linear matter P(k) at z = 0 [h/Mpc, (Mpc/h)³] |
| `results/` | raw CoLoRe output: κ map [`out_kappa_z000.fits`](results/out_kappa_z000.fits) (HEALPix, N_side = 512), source catalogs `out_srcs_s<N>_<i>.fits` (RA, DEC, Z_COSMO, …) and, when kept, CoLoRe's P(k) predictions and run logs |
| `cls/` | measured spectra and n(z), see below |
| `cmb_lensing/` | Simons Observatory κ reconstruction-noise curves (`nlkk_*.dat`) used for the noisy-κ tests |
| `plots/` | figures made from this box |
| [`fit_results.txt`](fit_results.txt) | best-fit $k_p$ and $A_{\rm bias}$ for κκ, QSO auto and κ × QSO (see below) |

### Spectra (`cls/`)

The spectra are produced by [`scripts/compute_colore_cls.py`](../../scripts/compute_colore_cls.py) (see `scripts/run_compute_cls*.sh`). Each `*.npz` holds the per-multipole spectrum; the `*_log.npz` versions are log-rebinned with `utils/binning.py`. The files in this box were made with earlier versions of that script, so some names differ from the current convention.

| File | Content |
| --- | --- |
| [`cl_kk.npz`](cls/cl_kk.npz) | CMB lensing auto-spectrum $C_\ell^{\kappa\kappa}$ per multipole (`ells, cls, err, nside`) |
| [`cl_kk_log.npz`](cls/cl_kk_log.npz) | same, log-rebinned |
| `cl_qq.npz`, `cl_qq_bin0.npz`, `cl_qq_bin1.npz` (+ `_log`) | QSO auto-spectrum, full sample and two true-z bins (old naming: q = QSO) |
| `cl_kq.npz`, `cl_kq_bin0.npz`, `cl_kq_bin1.npz` (+ `_log`) | κ × QSO cross-spectrum, same samples |
| `cl_gg_desi_qso.npz`, `cl_kg_desi_qso.npz` (+ `_log`) | QSO auto and κ × QSO with the later naming |
| `nz_hist_desi_qso.npz` | measured QSO n(z) |
| [`README.md`](cls/README.md) | map resolutions and the redshift bins used |

## CoLoRe input files

The [param.cfg](param.cfg) records the input files with the names they had when the box was run (`data/colore_files/nz/nz_desi_qso.txt`, `data/colore_files/bz-tunning/bz-tunning_desi_qso.txt`). They now live in [`../colore_setup/inputs/`](../colore_setup/inputs/): `nz/nz_desi_qso.txt` and `bias-input/bz-input_desi_qso.txt`. Their contents may have changed since this box was run. With the NORMALIZED bias model there is no threshold file.

## Cosmological parameters

[`cosmo_params.json`](cosmo_params.json) gives the cosmology used to generate the box (Ω_M, Ω_Λ, Ω_B, h, A_s, n_s, w, σ8, Σm_ν), the box maximum redshift `z_max` and the source-plane redshift `z_out`. The analysis reads it to build the theory.

## Fitting results

The file [fit_results.txt](fit_results.txt) contains the best-fit values of two parameters implemented in the theory model to agree with CoLoRe data. The fitting is done for the CMB lensing auto, QSO auto, CMB x QSO both for the full QSO sample and for the tomographic bins.
- $k_p$ parameter accounting for the finite cell size in CoLoRe
$$P_{\rm mm, lin}(k) \longrightarrow P_{\rm mm, lin}(k) \times e^{-(k/k_p)^2}$$
- $A_{\rm bias}$ parameter accounting for a rescaling of the input QSO bias to match CoLoRe's output bias. The bias model that was used in CoLoRe to run this box had this problem. Otherwise, theoretically, input and output bias $b(z)$ should be the same.
$$b_{\rm QSO}(z) \longrightarrow  A_{\rm bias} \times b_{\rm QSO}(z)$$
