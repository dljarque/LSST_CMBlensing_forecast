# CoLoRe setup

Everything needed to configure, build and run the CoLoRe boxes in [`../`](../).

## Input files (`inputs/`)

CoLoRe reads, for each source population, a redshift distribution, an input bias and a density threshold. The theory additionally uses a real bias. All files are two-column text (`z  value`) and, unless they carry the `_original` suffix, span z = 0–3.8 (the range simulated by CoLoRe; CoLoRe fails if an input file does not cover it).

| Folder | Files | Used by | Content |
| --- | --- | --- | --- |
| [`nz/`](inputs/nz/) | `nz_<sample>.txt` | CoLoRe | n(z) = dN/dΩ dz per deg² |
| [`bias-input/`](inputs/bias-input/) | `bz-input_<sample>.txt` | CoLoRe | input bias b_in(z) of the linear bias model |
| [`thresholds/`](inputs/thresholds/) | `tz_<sample>.txt` | CoLoRe | density threshold t(z) |
| [`bias/`](inputs/bias/) | `bz_<sample>.txt` | theory only | real bias b(z) that the tuned mock reproduces |

Samples: `desi_bgs`, `desi_lrg`, `desi_elg`, `desi_qso`, `lsst_gal`.

- `_original` suffix: file as received, before it was extended to z = 3.8 with [`plot_galaxy_distributions.ipynb`](plot_galaxy_distributions.ipynb).
- `tz_desi_bgs_new.txt`: attempt to lower the BGS threshold to reduce the mismatch between the mock bias and the real bias.
- `tz_desi_qso_old.txt`: earlier QSO tuning by Laura Casas; it goes with a different input bias and is not compatible with the current `bz-input_desi_qso.txt`. The current `tz_desi_qso.txt` comes from LyaCoLoRe.
- `bz_desi_lrg.txt` is 0.9 × the DESI FDR LRG bias, because the LRG input bias was tuned to reproduce that.

The origin of every file (DESI FDR, LyaCoLoRe, LSST Y10, tunings by Laura Casas) is listed in §3.3 of the [main README](../../README.md).

## Notebooks

- [`plot_galaxy_distributions.ipynb`](plot_galaxy_distributions.ipynb): extends the input files to z = 3.8, builds the real-bias files and plots n(z), b(z) and t(z) ([`plots/`](plots/)).
- [`tomographic_binning.ipynb`](tomographic_binning.ipynb): tests of the Gaussian photo-z binning used to build the `_bin0` samples.

## Building and running CoLoRe at NERSC (`configs/`)

- [`Makefile_linear`](configs/Makefile_linear): Makefile for the CoLoRe build with the **linear** bias model (`_BIAS_MODEL_3`), used for `test_box_4` onwards and the ensemble.
- [`Makefile_normalized`](configs/Makefile_normalized): Makefile for the build with the **normalized** bias model (`test_box_1`–`3`).
- `Install CoLoRe in NERSC (new).pdf` / `(old).pdf`: step-by-step installation of CoLoRe and its dependencies in a NERSC allocation.
- [`param_example.cfg`](configs/param_example.cfg): annotated example configuration.

The configuration of a new box is generated with [`scripts/generate_colore_config.py`](../../scripts/generate_colore_config.py) (example call in [`scripts/run_generate_colore_config.sh`](../../scripts/run_generate_colore_config.sh)).

## References (`references/`)

- `DESI_FDR_science.pdf`: DESI Experiment Part I (science, targeting and survey design); §2.4.2 gives the real bias of BGS, LRG and ELG.
- `DESI_DR1_fNL_2411.17623v2.pdf`: DESI DR1 f_NL analysis; source of the QSO inputs together with [LyaCoLoRe](https://github.com/igmhub/LyaCoLoRe/tree/colore-ql/colore-ql/colore_box/input_files).
