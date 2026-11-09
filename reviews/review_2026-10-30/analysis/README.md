# Review scripts (2026-10-30)

Run from this folder with the repo's python env (`numpy scipy camb healpy mcfit`).
Data are read from the repository's `data/` folder (`../../../data/`, override with `LSST_FORECAST_DATA`); `utils` is imported from the repository root. Scripts can be run from any directory, but they import each other by name, so run them as `python reviews/review_2026-10-30/analysis/tXX.py` or from this folder. `pw1024.npy` is the HEALPix nside=1024 pixel window.

- `t4_fits.py`  : replica of the notebook joint fit (test_box_7, ensemble_1, test_box_10, lsst)
- `t5_mech.py`  : sigma8-fixed A_gg vs A_kg; Hermite model joint fits
- `t6_cov_cuts.py` : scale-cut stability (r_min = 100 ... 15 Mpc/h), notebook vs theory covariance
- `t7_kk.py`, `t8.py` : kappa-kappa integrand fractions, sigma8 vs lmax, data/theory per ell range
- `t9_limber.py` : Limber vs exact Bessel at low ell
- `t10_covbias.py` : data-based covariance bias toy
- `t2_theory.py`, `t1_background.py` : CAMB vs CoLoRe background, k-mask, z-grid
