# CoLoRe box `test_box_4`

## Description

Fourth CoLoRe box, with CMB lensing, DESI galaxies and LSST galaxies output, run by Dani. First time using the Linear bias model of CoLoRe.

**Bias model:** this box was run with the **LINEAR** bias model of CoLoRe.

Source populations in `results/`: `srcs1` = LSST galaxies, `srcs2` = DESI BGS, `srcs3` = DESI LRG, `srcs4` = DESI ELG, `srcs5` = DESI QSO (so `out_srcs_s5_*` are the QSO catalogs).

## Contents

| Path | Content |
| --- | --- |
| [`param.cfg`](param.cfg) | CoLoRe configuration used for this run |
| [`cosmo_params.json`](cosmo_params.json) | cosmology of the box, `z_max` and the source-plane redshift `z_out` |
| [`pk_linear.dat`](pk_linear.dat) | input linear matter P(k) at z = 0 [h/Mpc, (Mpc/h)³] |
| `results/` | raw CoLoRe output: κ map [`out_kappa_z000.fits`](results/out_kappa_z000.fits) (HEALPix, N_side = 1024), source catalogs `out_srcs_s<N>_<i>.fits` (RA, DEC, Z_COSMO, …) and, when kept, CoLoRe's P(k) predictions and run logs |
| `cls/` | measured spectra and n(z), see below |
| `chains/` | MCMC chains from the fits to this box |
| `chains_old/`, `plots_old/` | chains and plots from earlier versions of the analysis |
| `plots/` | figures made from this box |

### Spectra (`cls/`)

The spectra are produced by [`scripts/compute_colore_cls.py`](../../scripts/compute_colore_cls.py) (see `scripts/run_compute_cls*.sh`). Each `*.npz` holds the per-multipole spectrum; the `*_log.npz` versions are log-rebinned with `utils/binning.py`. The files in this box were made with earlier versions of that script, so some names differ from the current convention.

| File | Content |
| --- | --- |
| [`cl_kk.npz`](cls/cl_kk.npz) | CMB lensing auto-spectrum $C_\ell^{\kappa\kappa}$ per multipole (`ells, cls, err, nside`) |
| [`cl_kk_log.npz`](cls/cl_kk_log.npz) | same, log-rebinned |
| `cl_gg_desi_bgs.npz`, `cl_gg_desi_bgs_bin0.npz` (+ `_log`) | DESI BGS auto-spectrum $C_\ell^{gg}$, full sample and tomographic bin 0 (shot noise subtracted, `shot_noise` stored) |
| `cl_kg_desi_bgs.npz`, `cl_kg_desi_bgs_bin0.npz` (+ `_log`) | κ × DESI BGS cross-spectrum $C_\ell^{\kappa g}$, same samples |
| `nz_desi_bgs.npz`, `nz_desi_bgs_bin0.npz` | measured n(z) of the full sample and of bin 0 |
| `cl_gg_desi_lrg.npz`, `cl_gg_desi_lrg_bin0.npz` (+ `_log`) | DESI LRG auto-spectrum $C_\ell^{gg}$, full sample and tomographic bin 0 (shot noise subtracted, `shot_noise` stored) |
| `cl_kg_desi_lrg.npz`, `cl_kg_desi_lrg_bin0.npz` (+ `_log`) | κ × DESI LRG cross-spectrum $C_\ell^{\kappa g}$, same samples |
| `nz_desi_lrg.npz`, `nz_desi_lrg_bin0.npz` | measured n(z) of the full sample and of bin 0 |
| `cl_gg_desi_elg.npz`, `cl_gg_desi_elg_bin0.npz` (+ `_log`) | DESI ELG auto-spectrum $C_\ell^{gg}$, full sample and tomographic bin 0 (shot noise subtracted, `shot_noise` stored) |
| `cl_kg_desi_elg.npz`, `cl_kg_desi_elg_bin0.npz` (+ `_log`) | κ × DESI ELG cross-spectrum $C_\ell^{\kappa g}$, same samples |
| `nz_desi_elg.npz`, `nz_desi_elg_bin0.npz` | measured n(z) of the full sample and of bin 0 |
| `cl_gg_desi_qso.npz`, `cl_gg_desi_qso_bin0.npz` (+ `_log`) | DESI QSO auto-spectrum $C_\ell^{gg}$, full sample and tomographic bin 0 (shot noise subtracted, `shot_noise` stored) |
| `cl_kg_desi_qso.npz`, `cl_kg_desi_qso_bin0.npz` (+ `_log`) | κ × DESI QSO cross-spectrum $C_\ell^{\kappa g}$, same samples |
| `nz_desi_qso.npz`, `nz_desi_qso_bin0.npz` | measured n(z) of the full sample and of bin 0 |
| `cl_gg_lsst_gal.npz`, `cl_gg_lsst_gal_bin0.npz` (+ `_log`) | LSST galaxies auto-spectrum $C_\ell^{gg}$, full sample and tomographic bin 0 (shot noise subtracted, `shot_noise` stored) |
| `cl_kg_lsst_gal.npz`, `cl_kg_lsst_gal_bin0.npz` (+ `_log`) | κ × LSST galaxies cross-spectrum $C_\ell^{\kappa g}$, same samples |
| `nz_lsst_gal.npz`, `nz_lsst_gal_bin0.npz` | measured n(z) of the full sample and of bin 0 |
| [`README.md`](cls/README.md) | map resolutions used |

## CoLoRe input files

The [param.cfg](param.cfg) records the input files with the names they had when the box was run (`data/colore_files/...`). They now live in [`../colore_setup/inputs/`](../colore_setup/inputs/) under new names, and their contents may have changed since this box was run:

| Name in `param.cfg` | Current location |
| --- | --- |
| `colore_files/nz/nz_<sample>.txt` | [`inputs/nz/nz_<sample>.txt`](../colore_setup/inputs/nz/) |
| `colore_files/bz-tunning/bz-tunning_<sample>.txt` | [`inputs/bias-input/bz-input_<sample>.txt`](../colore_setup/inputs/bias-input/) (CoLoRe input bias) |
| `colore_files/Tz/Tz_<sample>.txt` | [`inputs/thresholds/tz_<sample>.txt`](../colore_setup/inputs/thresholds/) |

The real bias b(z) used in the theory is in [`inputs/bias/bz_<sample>.txt`](../colore_setup/inputs/bias/).

## Cosmological parameters

[`cosmo_params.json`](cosmo_params.json) gives the cosmology used to generate the box (Ω_M, Ω_Λ, Ω_B, h, A_s, n_s, w, σ8, Σm_ν), the box maximum redshift `z_max` and the source-plane redshift `z_out`. The analysis reads it to build the theory.
