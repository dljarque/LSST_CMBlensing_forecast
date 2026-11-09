# CoLoRe box `test_box_1`

## Description

First CoLoRe box with CMB lensing output, run by Laura. The source plane redshift `z_out` is equal to the box maximum redshift `z_max`. This issue is solved in the second box.

**Bias model:** this box was run with the **NORMALIZED** bias model of CoLoRe.

## Contents

| Path | Content |
| --- | --- |
| [`param.cfg`](param.cfg) | CoLoRe configuration used for this run |
| [`cosmo_params.json`](cosmo_params.json) | cosmology of the box, `z_max` and the source-plane redshift `z_out` |
| `results/` | raw CoLoRe output: κ map [`out_kappa_z000.fits`](results/out_kappa_z000.fits) (HEALPix, N_side = 512) |
| `cls/` | measured spectra and n(z), see below |
| `plots/` | figures made from this box |

### Spectra (`cls/`)

The spectra are produced by [`scripts/compute_colore_cls.py`](../../scripts/compute_colore_cls.py) (see `scripts/run_compute_cls*.sh`). Each `*.npz` holds the per-multipole spectrum; the `*_log.npz` versions are log-rebinned with `utils/binning.py`. The files in this box were made with earlier versions of that script, so some names differ from the current convention.

| File | Content |
| --- | --- |
| [`cl_kk.npz`](cls/cl_kk.npz) | CMB lensing auto-spectrum $C_\ell^{\kappa\kappa}$ per multipole (`ells, cls, err, nside`) |
| [`cl_kk_log.npz`](cls/cl_kk_log.npz) | same, log-rebinned |

## CoLoRe input files

The [param.cfg](param.cfg) is the one used by Laura Casas at NERSC. It points to her own input files: the QSO n(z) `Nz_qso_130618_2_colore1_hZs.txt` and the QSO bias `Bz_qso_G18.txt` (not stored in this repository), and to the LyaCoLoRe linear P(k) `PlanckDR12_kmax_matterpower_z0.dat`.

## Cosmological parameters

[`cosmo_params.json`](cosmo_params.json) gives the cosmology used to generate the box (Ω_M, Ω_Λ, Ω_B, h, A_s, n_s, w, σ8, Σm_ν), the box maximum redshift `z_max` and the source-plane redshift `z_out`. The analysis reads it to build the theory.
