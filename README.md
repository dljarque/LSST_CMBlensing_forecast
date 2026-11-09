# Forecast on LSST Y10 × CMB lensing

Forecast for the cross-correlation of LSST year-10 galaxies with CMB lensing, where the LSST redshift distribution and bias, n(z) b(z), are calibrated through cross-correlations with DESI galaxies instead of the LSST auto-spectrum. The forecast is compared with DESI DR3 × CMB lensing.

The current stage of the project validates the inference pipeline on **CoLoRe lognormal simulations**. We break the degeneracy between galaxy bias and σ8 by combining the galaxy auto-spectrum C_ℓ^gg with the galaxy × CMB-lensing cross-spectrum C_ℓ^κg (with C_ℓ^κκ as a consistency check). The aim is to recover the fiducial σ8 of the mocks, and to understand every systematic offset before moving to more realistic modelling and survey settings.

> Author: Daniel L. Jarque (IFAE). The first two CoLoRe boxes and the original galaxy-bias tuning of the samples are by Laura Casas.
> This repository continues the earlier `cmb_lensing_lya` repository; the research line no longer involves the Lyman-α forest.

---

## Contents

1. [Science goal](#1-science-goal)
2. [Current status](#2-current-status)
3. [Method](#3-method)
4. [Repository structure](#4-repository-structure)
5. [Workflow: from a CoLoRe box to σ8](#5-workflow-from-a-colore-box-to-σ8)
6. [Data products and file formats](#6-data-products-and-file-formats)
7. [Conventions](#7-conventions)
8. [Simulation history](#8-simulation-history)
9. [Environment](#9-environment)
10. [Known issues and housekeeping](#10-known-issues-and-housekeeping)
11. [References](#11-references)

---

## 1. Science goal

**Final goal (the paper).** A forecast of LSST Y10 × CMB lensing, in which the clustering amplitude of LSST galaxies, b(z) n(z), is calibrated with the cross-correlation against DESI tracers (BGS, LRG, ELG, QSO) rather than with the LSST auto-spectrum. The forecast is compared against DESI DR3 × CMB lensing.

**Why it works.** In the linear, Limber regime and for a single tomographic bin,

```
C_ℓ^gg  ∝  A_bias² σ8²
C_ℓ^κg  ∝  A_bias  σ8²
```

so b and σ8 are fully degenerate in either spectrum alone, but the combination fixes both: σ8 ∝ C_ℓ^κg / √C_ℓ^gg. C_ℓ^κκ ∝ σ8² gives an independent check.

**Current stage (validation on mocks).** Before forecasting, the pipeline must recover the input cosmology of simulations whose truth is known exactly. We use full-sky CoLoRe lognormal boxes containing a CMB-lensing convergence map and five galaxy samples (LSST galaxies and the four DESI tracers), and fit (A_bias, σ8, A_sn) to the measured spectra. The paper is a forecast, not a simulation paper, so mock-specific modelling is kept to a minimum.

## 2. Current status

As of November 2026 (ten-box ensemble `ensemble_1`–`ensemble_10`):

- **The pipeline is fixed.** It uses the CoLoRe-consistent background and P(k), the exact grid-smoothing and pixel windows, 40 log bins in ℓ = 30–1000 selected by their upper edge, and a Knox covariance from theory. For an average of N boxes the covariance is divided by N.
- **Single boxes recover σ8 within ≈ 2σ** for every sample at its safe scale cut.
- **Averaging ten boxes reveals two systematics, both making σ8 low:**
  1. *Non-linear galaxy excess in C_ℓ^gg.* CoLoRe makes the galaxy density a non-linear function of a Gaussian field (lognormal + bias + threshold). This adds an almost-white excess power that a constant shot-noise term A_sn only partly absorbs. It is computed exactly with a Hermite (Mehler) expansion (`reviews/review_2026-10-30/colore_theory.py`). It is largest for high-threshold samples (LRG, BGS) and negligible for LSST.
  2. *κ deficit.* The CoLoRe κκ is ≈ 0.6% below theory at ℓ > 200 and 1.5–4% below at ℓ < 100, in every box. A leave-one-out mock transfer function T(ℓ) corrects it (κg × √T).
- **Safe scale cuts**, defined as the smallest r_min where the constant-A_sn model bias is below 0.5σ:

  | Sample | z_eff | Ten-box precision r_min [Mpc/h] | Single-box precision r_min [Mpc/h] |
  | --- | --- | --- | --- |
  | DESI BGS | 0.4 | 100 | 80 |
  | DESI LRG | 1.0 | 150 | 100 |
  | DESI ELG | 1.2 | 80 | 60 |
  | DESI QSO | 2.0 | 100 | 60 |
  | LSST gal | 1.0 | 50 | 40 |

  ℓ_max = 2π χ(z_eff) h / r_min − ½. The analysis notebook currently uses the single-box cuts.
- **Open issues:** an unexplained −0.8% shift of LSST σ8 between r_min = 60 and 40 Mpc/h; the cause of the κ deficit; Knox overestimating the σ8 error by 20–40% for most samples at small scales; and whether to model the excess (Hermite) or retune the mocks.

The full discussion, with figures, is in the meeting notes “σ8 bias in CoLoRe mocks — meeting notes (ten ensemble boxes)” and in `reviews/review_2026-10-30/`.

## 3. Method

### 3.1 Simulations: CoLoRe

[CoLoRe](https://github.com/damonge/CoLoRe) generates lognormal light-cone simulations on a 4096³ grid up to z_max = 3.8. The runs use the **linear bias model** (`_BIAS_MODEL_3`) with a density threshold. For each cell:

1. A Gaussian field δ_G is drawn from the input P(k), normalised to σ8 and smoothed with exp(−k²R²), R = `r_smooth` = 2 Mpc/h.
2. Lognormal transform: δ_LN = exp(D δ_G − D² σ_G²/2) − 1.
3. Galaxy weight: w = max(1 + b_in(z) δ_LN, 0) if δ_LN ≥ t(z), else 0. b_in is the *input* bias (`bz-input_*`) and t(z) the threshold (`tz_*`). Both were tuned so that the large-scale bias of the mock matches the *real* bias b(z) (`bz_*`).
4. Normalisation to the target n(z) and Poisson sampling.
5. κ is integrated along the line of sight from the linear Gaussian potential (`n_lensing` = 11 shells), on a HEALPix map with N_side = 1024, with the source plane at z_out = z_* ≈ 1090.

Fiducial cosmology (all ensemble boxes): Ω_M = 0.3147, Ω_Λ = 0.6853, Ω_B = 0.04904, h = 0.6731, n_s = 0.9655, A_s = 2.1 × 10⁻⁹, w = −1, Σm_ν = 0. This gives σ8 = 0.82509.

### 3.2 Galaxy samples

Five source populations are simulated in every box (`srcs1`–`srcs5`, in this order):

| Population | Sample | Catalog prefix |
| --- | --- | --- |
| srcs1 | LSST galaxies | `out_srcs_s1` |
| srcs2 | DESI BGS | `out_srcs_s2` |
| srcs3 | DESI LRG | `out_srcs_s3` |
| srcs4 | DESI ELG | `out_srcs_s4` |
| srcs5 | DESI QSO | `out_srcs_s5` |

For each sample, one tomographic bin (`<sample>_bin0`) is built from the true redshifts with a Gaussian photo-z: z_photo = z_true + N(0, σ_z), keeping galaxies with |z_photo − z_mean| < Δz/2 (their true-z n(z) is what the theory uses). The bins use Δz = 0.2 and σ_z = 0.05, centred at z_mean = 0.4 (BGS), 1.0 (LRG), 1.2 (ELG), 2.0 (QSO) and 1.0 (LSST). The photo-z draw uses a fixed seed (seed = 0), so the gg and κg runs select exactly the same galaxies.

### 3.3 Input files and their provenance

CoLoRe needs, for each sample, an n(z), an input bias b_in(z) and a threshold t(z) (`data/colore_setup/inputs/`). The theory uses a separate **real** bias b(z), the one the tuned mock is meant to reproduce.

| Sample | n(z) (`nz/`) | Real bias b(z) (`bias/bz_*`) | Input bias b_in(z) (`bias-input/bz-input_*`) | Threshold t(z) (`thresholds/tz_*`) |
| --- | --- | --- | --- | --- |
| DESI BGS | sent to Laura Casas (original source unknown) | DESI FDR science paper (§2.4.2) | tuned by Laura Casas | tuned by Laura Casas |
| DESI LRG | sent to Laura Casas (original source unknown) | **0.9 ×** the DESI FDR bias (see note) | tuned by Laura Casas to reproduce 0.9 × the FDR bias | tuned by Laura Casas |
| DESI ELG | sent to Laura Casas (original source unknown) | DESI FDR science paper (§2.4.2) | tuned by Laura Casas | tuned by Laura Casas |
| DESI QSO | [LyaCoLoRe](https://github.com/igmhub/LyaCoLoRe/tree/colore-ql/colore-ql/colore_box/input_files) (`colore-ql` branch) and the DESI DR1 f_NL paper | provided by Laura Casas | LyaCoLoRe tuning | LyaCoLoRe (`tz_desi_qso.txt`) |
| LSST gal | LSST Y10, provided to Laura Casas by Qianjun Hang | LSST Y10, provided to Laura Casas by Qianjun Hang | tuned by Laura Casas | tuned by Laura Casas |

- **DESI LRG exception.** The LRG input bias was tuned to reproduce 0.9 × the FDR bias, so `bz_desi_lrg.txt` must be the FDR bias times 0.9. The theory then matches the mock; with the plain FDR bias it does not.
- **LSST.** The n(z) and b(z) are an LSST year-10 forecast. They were obtained about a year ago, probably from a forecast or directly from LSST members.
- **Alternative threshold files.** `tz_desi_qso_old.txt` is an earlier QSO tuning by Laura Casas that goes with a different input bias from the current one. It is not compatible with `bz-input_desi_qso.txt`. `tz_desi_bgs_new.txt` is an attempt to lower the BGS threshold to reduce the mismatch between the mock bias and the real bias.
- `_original` files are the inputs as received, before they were extended to z = 3.8 (see §7).

### 3.4 Measured spectra

`utils/colore_cls.py` builds full-sky HEALPix maps and runs `healpy.anafast`:

- **C_ℓ^κκ** from `out_kappa_z000.fits`;
- **C_ℓ^gg** from the overdensity map δ = N/N̄ − 1, with the shot noise 1/n̄ (per steradian) subtracted and stored;
- **C_ℓ^κg** from the κ map × δ.

Per-ℓ Gaussian (Knox) errors are stored together with log-rebinned versions (`*_log.npz`, see `utils/binning.py`).

### 3.5 Theory

The Limber integrals are computed with the CAMB linear matter power, using the same cosmology as the box (`cosmo_params.json`):

- Galaxy kernel: W_g(z) = b(z) n(z) H(z)/c, with b(z) the **real** bias (`bz_*`) and n(z) the measured n(z) of the tomographic bin.
- Lensing kernel: W_κ(χ) = (3/2) Ω_M H₀² χ (χ_* − χ)/(χ_* a).
- Grid smoothing: P(k) × exp(−(k/k_p)²). For gg and κg, k_p = h/(2 Mpc) (R = 2 Mpc/h). For κκ, k_p = 0.2957 Mpc⁻¹, fitted on κκ with σ8 fixed.
- HEALPix pixel window: pw² for gg, pw for κg.

`reviews/review_2026-10-30/colore_theory.py` re-implements the background, growth and σ8 normalisation exactly as CoLoRe does, plus the Hermite expansion of the galaxy two-point function (the source of the gg excess).

### 3.6 Likelihood and fit

- **Model:** C^gg = A_bias² (σ8/σ8,fid)² C^gg,lin + A_sn/n̄ and C^κg = A_bias (σ8/σ8,fid)² C^κg,lin. A_bias is a free amplitude (it is not expected to equal 1). A_sn absorbs the white part of the non-linear excess.
- **Data vector:** 40 log bins, ℓ = 30–1000, keeping bins whose upper edge is ≤ ℓ_max(r_min).
- **Covariance:** theory Knox, including the gg–κg cross-covariance. It is built at the fiducial, then rebuilt once at the best-fit (A_bias, A_sn). For an average of N boxes it is divided by N.
- **Fit:** `scipy.optimize` minimisation, then MCMC with [zeus](https://github.com/minaskar/zeus). Flat priors: 0.5 < σ8, A_bias < 1.5 and −5 < A_sn < 100. Posteriors are plotted with GetDist.
- **κ deficit (optional):** κg theory × √T(ℓ), κκ × T(ℓ). T is a cubic in ln ℓ fitted to ⟨data/theory⟩ of κκ over the other boxes (leave-one-out).

### 3.7 Validation tests

- **Noiseless test:** fit a zero-noise data vector made of linear theory with the mock's own bias plus the Hermite excess. Any σ8 shift is then pure model error. This test sets the safe scale cuts.
- **Closure:** noiseless bias + κ-deficit shift is compared with the observed ten-box offset. It closes within 1σ for all four DESI samples.
- **Box-to-box scatter vs Knox**, **χ² at best fit vs at fiducial σ8**, and **paired r_min differences** (which cancel most cosmic variance).

## 4. Repository structure

```
LSST_CMBlensing_forecast/
├── README.md                 this file
├── LICENSE                   CC0 1.0
├── .gitignore
├── data/
│   ├── colore_setup/         everything needed to configure and run CoLoRe
│   │   ├── inputs/
│   │   │   ├── nz/           n(z) = dN/dΩ dz per deg² for each sample (nz_<sample>.txt)
│   │   │   ├── bias/         REAL bias b(z) used in the theory (bz_<sample>.txt)
│   │   │   ├── bias-input/   CoLoRe INPUT bias b_in(z) (bz-input_<sample>.txt)
│   │   │   └── thresholds/   CoLoRe density thresholds t(z) (tz_<sample>.txt)
│   │   ├── configs/          NERSC Makefiles (linear / normalized bias builds),
│   │   │                     param_example.cfg, CoLoRe-on-NERSC install guides (PDF)
│   │   ├── references/       DESI FDR science paper (real b(z) of BGS, LRG, ELG) and DESI DR1 f_NL paper (QSO)
│   │   ├── plots/            input n(z), b(z), t(z) figures
│   │   ├── plot_galaxy_distributions.ipynb   extends inputs to z = 3.8, builds real-bias files
│   │   └── tomographic_binning.ipynb         tests of the Gaussian photo-z binning
│   ├── ensemble_1 … ensemble_10/   production boxes (same config, seed = N × 10⁵)
│   │   ├── param.cfg               CoLoRe configuration
│   │   ├── cosmo_params.json       cosmology, z_max, z_out
│   │   ├── pk_linear.dat           input linear P(k) at z = 0 [h/Mpc, (Mpc/h)³]
│   │   ├── README.md
│   │   ├── results/                raw CoLoRe output (κ map, source catalogs, P(k) predictions, logs)
│   │   ├── cls/                    measured spectra and n(z) (see §6)
│   │   └── plots/
│   ├── test_box_1 … test_box_11/   development boxes (see §8); same layout, each with its own README
│   ├── test_box_6_legacy, test_box_6_master/   CoLoRe-branch cross-check
│   └── test_box_nside/             κκ vs N_side comparison plots
├── scripts/
│   ├── generate_colore_config.py   CAMB P(k) + param.cfg + cosmo_params.json for a box
│   ├── run_generate_colore_config.sh   example call (ensemble_1)
│   ├── compute_colore_cls.py       CLI: cmb_lensing_auto | sources_auto | cmblensing_x_sources
│   ├── run_compute_cls.sh          gg + κg for one sample and box (tomographic bin)
│   ├── run_compute_cls_kk.sh       κκ for one box
│   └── run_compute_cls_gg.sh, run_compute_cls_kg.sh   single-spectrum variants (older)
├── utils/
│   ├── colore_cls.py         catalog merging, tomographic binning, map making, anafast, saving
│   ├── binning.py            hybrid linear/log ℓ binning (hybrid_bin_cls, make_hybrid_bin_matrix)
│   └── colore_param_cfg.py   CAMB cosmology, P(k) output, CoLoRe param.cfg writer
├── notebooks/
│   ├── cosmological_analysis.ipynb   MAIN ANALYSIS: theory, covariance, fits, MCMC, plots
│   ├── plot_fitting_cls_fixedkp.ipynb, plot_fitting_cls_kp.ipynb   earlier gg/κg validation (fixed vs free k_p)
│   ├── plot_cmblensing_cls.ipynb     earlier κκ validation, k_p fit, SO noise test
│   ├── plot_quasar_auto_cross.ipynb  first QSO auto / cross validation
│   └── plots/                        kernels, κ and QSO maps
└── reviews/
    └── review_2026-10-30/            pipeline review (Oct–Nov 2026)
        ├── colore_theory.py          CoLoRe-consistent theory + Hermite galaxy field
        ├── colore_cls.patch          fixes applied to utils/colore_cls.py (photo-z seed = 0, Knox factor)
        └── analysis/                 numbered diagnostic scripts t1 … t43 (see its README)
```

### Key review scripts (`reviews/review_2026-10-30/analysis/`)

| Script | Purpose |
| --- | --- |
| `common.py`, `pipeline.py`, `lognormal.py`, `win.py` | shared loaders, CoLoRe cosmology, Limber, Hermite galaxy field, windows |
| `t30.py` | reference implementation of the fixed pipeline (binning, theory, Knox, joint fit) |
| `t27.py`, `t29.py` | noiseless test of the constant-A_sn model |
| `t36.py`–`t38.py` | six-box QSO / LSST fits and leave-one-out T(ℓ) (superseded) |
| `t41_ten_boxes.py` | ten-box QSO + LSST fits, κκ ratios, T(ℓ), scatter vs Knox |
| `t42_kk_chi2.py` | κκ data/theory per bin with and without T(ℓ); χ² at best fit vs fiducial |
| `t43_all_samples.py` | all five samples: r_min scan of the ten-box fits + noiseless test (safe cuts) |

## 5. Workflow: from a CoLoRe box to σ8

All commands run from the repository root.

**1. Prepare the inputs (once).** Edit the files in `data/colore_setup/inputs/`. Files without the `_original` suffix are extended to z = 3.8 (CoLoRe requires full coverage of the simulated range); `plot_galaxy_distributions.ipynb` does the extension.

**2. Generate the box configuration.**

```bash
bash scripts/run_generate_colore_config.sh      # edit test_box, --seed and paths first
```

This runs CAMB with the given cosmology (Σm_ν = 0), writes `pk_linear.dat`, `param.cfg` (z_out = z_* from CAMB) and `cosmo_params.json` into `data/<box>/`. `--prefix_out` and the input paths must point to where CoLoRe will run (NERSC `$PSCRATCH`).

**3. Run CoLoRe on NERSC.** Build CoLoRe with the linear-bias Makefile (`data/colore_setup/configs/Makefile_linear`; installation guide in the same folder), then run it with `param.cfg`. Copy `results/` back to `data/<box>/results/`.

**4. Measure the spectra.**

```bash
bash scripts/run_compute_cls_kk.sh   # set test_box
bash scripts/run_compute_cls.sh      # set test_box, sources (out_srcs_sN), filename, zmean, dz, sigma_z
```

or call the CLI directly:

```bash
python scripts/compute_colore_cls.py sources_auto \
    --sources_prefix ./data/ensemble_1/results/out_srcs_s5 --nside 1024 \
    --lmin 2 --nbins 30 --out_dir ./data/ensemble_1/cls --filename desi_qso_bin0 \
    --zmean 2.0 --dz 0.2 --sigma_z 0.05
python scripts/compute_colore_cls.py cmblensing_x_sources \
    --kappa_map ./data/ensemble_1/results/out_kappa_z000.fits \
    --sources_prefix ./data/ensemble_1/results/out_srcs_s5 \
    --lmin 2 --nbins 30 --out_dir ./data/ensemble_1/cls --filename desi_qso_bin0 \
    --zmean 2.0 --dz 0.2 --sigma_z 0.05
python scripts/compute_colore_cls.py cmb_lensing_auto \
    --kappa_map ./data/ensemble_1/results/out_kappa_z000.fits --lmin 2 --nbins 30 \
    --out_dir ./data/ensemble_1/cls
```

The analysis re-bins the unbinned `cl_*.npz` itself; the `*_log.npz` files are for quick plots.

**5. Fit a single box.** In `notebooks/cosmological_analysis.ipynb`, set `box` and `galaxy_catalog`, then run the sections in order: theory → data → covariance → κκ (k_p, σ8) → joint gg + κg (minimiser, then MCMC) → single-spectrum fits → plots. `r_min_by_sample` holds the scale cuts.

**6. Ensemble averages and systematics.** Use the review scripts, for example:

```bash
python reviews/review_2026-10-30/analysis/t43_all_samples.py LSST   # one sample per call
```

The scripts find `data/` and `utils/` from their own location, so they can be run from any directory. Set `LSST_FORECAST_DATA` to read the boxes from somewhere else (e.g. `$PSCRATCH`). Output JSON files are written next to the scripts.

## 6. Data products and file formats

Inside `data/<box>/cls/`:

| File | Keys | Content |
| --- | --- | --- |
| `cl_kk.npz` | `ells, cls, err, nside` | κκ per ℓ (ℓ = 0 … 3 N_side − 1) |
| `cl_gg_<sample>_bin0.npz` | `ells, cls, err, nside, shot_noise` | gg per ℓ, shot noise already subtracted; `shot_noise` = 1/n̄ [sr] |
| `cl_kg_<sample>_bin0.npz` | `ells, cls, err, nside, shot_noise` | κg per ℓ |
| `cl_*_log.npz` | `ells, cls, err, …` | hybrid linear/log rebinned versions |
| `nz_<sample>_bin0.npz` | `z_nz, nz, z_min, z_max, z_mean, dz, sigma_z` | measured n(z) of the tomographic bin |
| `nz_<sample>.npz` | `z_nz, nz, z_min, z_max, z_mean` | measured n(z) of the full sample |

Inside `data/<box>/results/` (CoLoRe output): `out_kappa_z000.fits` (κ, HEALPix, HDU 1), `out_srcs_s<N>_<i>.fits` (catalogs: RA, DEC, Z_COSMO, …), `out_pk_srcs_pop<N>_z<z>.txt` (CoLoRe's own P(k) predictions), `out_gbias.txt`, `out_params.cfg` and the job logs.

## 7. Conventions

- **Bias files.** `bz-input_*` is the CoLoRe *input* bias (small b); `bz_*` is the *real* bias of the sample. The theory deliberately uses the real bias `bz_*`; A_bias absorbs any constant mismatch with the mock's effective bias.
- **`_original` suffix:** input file as received, not yet extended to z = 3.8. Files without the suffix are the extended ones that CoLoRe reads. `tz_desi_bgs_new.txt` and `tz_desi_qso_old.txt` are alternative threshold tunings kept for reference (see §3.3).
- **Ensemble boxes:** `data/ensemble_N` with CoLoRe seed = N × 10⁵ (spaced to avoid correlated seeds). Configuration otherwise identical.
- **Neutrinos:** from `test_box_11` onwards Σm_ν = 0 in both the input P(k) and the analysis. CoLoRe does not model neutrinos.
- **Tomographic bin name:** `<survey>_<tracer>_bin<k>` (e.g. `desi_qso_bin0`, `lsst_gal_bin0`).
- **Units:** P(k) files in h/Mpc and (Mpc/h)³; theory internally in Mpc⁻¹ and Mpc³; r_min in Mpc/h.

## 8. Simulation history

| Box | Bias model | Purpose / change |
| --- | --- | --- |
| `test_box_1` | normalized | first κ box (Laura); z_out = z_max (bug) |
| `test_box_2` | normalized | z_out ≠ z_max fixed (Laura) |
| `test_box_3` | normalized | κ + QSO (first run by Dani) |
| `test_box_4` | linear | κ + DESI + LSST; first linear-bias run |
| `test_box_5` | linear | new BGS threshold; lensing N_side = 512 |
| `test_box_6_legacy` / `_master` | linear | CoLoRe master vs arbitrary-source-planes branch (master broken) |
| `test_box_7` | linear | P(k) table 500 → 1000 points; insensitivity checks |
| `test_box_8` / `test_box_9` | linear | κ N_side × 2 / × ½ |
| `test_box_10` | linear | r_smooth = 1 Mpc/h (its `param.cfg` says 2, see its README) |
| `test_box_11` | linear | Σm_ν = 0 in the input P(k) |
| `ensemble_1`–`ensemble_10` | linear | production ensemble, current input bias and thresholds |

## 9. Environment

Python ≥ 3.12 (developed with 3.14) with `numpy`, `scipy`, `matplotlib`, `healpy`, `astropy`, `camb`, `mcfit`, `zeus-mcmc` and `getdist`. The shell scripts call `/home/dljarque/miniforge3/envs/IFAE/bin/python`; change `PYTHON=` to your interpreter. Scripts and notebooks import `utils` relative to the repository root (`sys.path.append('.')` or `'../'`), so run scripts from the root and notebooks from `notebooks/`.

CoLoRe itself runs on NERSC (Perlmutter); see `data/colore_setup/configs/` for the Makefiles and the installation guide.

## 10. Known issues and housekeeping

- **Large files in git.** `.gitignore` only excludes caches, so κ maps (≈ 50 MB each) and source catalogs end up in the history. Consider ignoring `data/*/results/` (or using Git LFS) and keeping only `cls/`, configs and READMEs under version control.
- **Covariance.** The Knox covariance is not uniformly conservative (see §2); review it before quoting forecast errors.

## 11. References

- CoLoRe: D. Alonso et al., <https://github.com/damonge/CoLoRe>
- DESI Collaboration, *The DESI Experiment Part I: Science, Targeting, and Survey Design* (`data/colore_setup/references/DESI_FDR_science.pdf`; §2.4.2 for the real bias of BGS, LRG and ELG)
- DESI DR1 f_NL analysis, arXiv:2411.17623 (`data/colore_setup/references/DESI_DR1_fNL_2411.17623v2.pdf`; QSO inputs, with LyaCoLoRe)
- LyaCoLoRe input files for the QSO: <https://github.com/igmhub/LyaCoLoRe/tree/colore-ql/colore-ql/colore_box/input_files>
- CAMB: <https://camb.info>; zeus: Karamanis et al. 2021; GetDist: Lewis 2019

License: CC0 1.0 Universal (see `LICENSE`).
