# Notebook-faithful pipeline (notebooks/cosmological_analysis.ipynb, Nov 2026 version):
# CAMB linear P(k,z) (delta_nonu, mnu from cosmo_params), CAMB chi(z), H(z); Limber with Simpson weights on 1000 z;
# CoLoRe windows exp(-k^2/kp^2) with kp = h/sqrt(R^2 + c dx^2), c = 1/12 (gg), 13/84 (kg), 19/84 (kk);
# pixel window: gg pw^2 (signal and A_sn term), kg pw, kk none; 40 hybrid(log) bins l=30-1000;
# fit bins: B@ell >= 30 and highest ell in bin <= lmax_fit, lmax_fit = int(2pi/r_min * chi(z_mean_bin) h - 0.5);
# two-pass Knox covariance (build_theory_cov), divided by N for averaged spectra.
import os, sys, json, warnings
import numpy as np, camb, healpy as hp
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
warnings.filterwarnings('ignore')
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..'))
sys.path.insert(0, REPO)
from utils.binning import make_hybrid_bin_matrix
D = os.path.join(os.environ.get('LSST_FORECAST_DATA', os.path.join(REPO, 'data')), '')
C_KMS = 299792.458
LMIN, LMAX, NB = 30, 1000, 40
ells = np.arange(LMIN, LMAX + 1)
B = make_hybrid_bin_matrix(ells, nbins=NB, ell_min=LMIN, ell_max=LMAX)
ells_hi = np.array([ells[row > 0].max() for row in B]); lb = B @ ells
BOXES = [f'ensemble_{i}' for i in range(1, 11)]
SAMPLES = {'BGS': 'desi_bgs', 'LRG': 'desi_lrg', 'ELG': 'desi_elg', 'QSO': 'desi_qso', 'LSST': 'lsst_gal'}

cd = json.load(open(D + 'ensemble_1/cosmo_params.json'))
h = cd['h']; s0 = cd['sigma_8']; z_max = cd['z_max']; om = cd['omega_M']; H0 = 100 * h

def camb_results(cd):
    mnu = cd.get('mnu', 0.06)
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * cd['h'], ombh2=cd['omega_B'] * cd['h'] ** 2,
                    omch2=(cd['omega_M'] - cd['omega_B']) * cd['h'] ** 2 - mnu / 93.14, mnu=mnu, omk=cd.get('omk', 0.))
    p.set_dark_energy(w=cd['w']); p.InitPower.set_params(ns=cd['ns'], As=cd['As'])
    zp = (np.exp(np.log(z_max + 1) * np.linspace(0, 1, 256)) - 1)[::-1]
    p.WantTransfer = True; p.set_matter_power(redshifts=zp, kmax=20, nonlinear=False)
    return camb.get_results(p)
R = camb_results(cd)
PK = R.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False, var1='delta_nonu', var2='delta_nonu')
BG = os.environ.get('NB_BG', 'colore')   # 'camb' = notebook as is; 'colore' = CoLoRe background (no radiation)
z = np.linspace(1e-4, z_max, 1000)
if BG == 'colore':
    from common import CoLoReCosmo
    CC = CoLoReCosmo(cd, D + 'ensemble_1/pk_linear.dat')
    chi = CC.chi(z); Hz = CC.H(z); chi_star = float(CC.chi(cd['z_out'])); chi_of = lambda zz: float(CC.chi(zz))
else:
    chi = R.comoving_radial_distance(z); Hz = R.hubble_parameter(z); chi_star = R.comoving_radial_distance(cd['z_out'])
    chi_of = lambda zz: R.comoving_radial_distance(zz)
sw = np.array([simpson(np.eye(len(z))[i], x=z) for i in range(len(z))])
k = (ells[None, :] + .5) / chi[:, None]
if BG == 'colore':   # pk_linear.dat (sigma8-normalised, as fed to CoLoRe) x CoLoRe growth
    P = CC.P(np.broadcast_to(z[:, None], k.shape).ravel(), k.ravel()).reshape(k.shape)
else:
    P = PK.P(np.broadcast_to(z[:, None], k.shape), k, grid=False)
Wk = 1.5 * (H0 / C_KMS) ** 2 * om * (chi_star - chi) / chi_star * chi * (1 + z)
dx = chi_of(z_max) * 2 / 4096 * h * (1 + (2 / 4096 if BG == 'colore' else 0))   # Mpc/h; CoLoRe: l_box = 2 r_max (1 + 2/n_grid)
R_SM = 2.0
def kp_of(c, Rs=R_SM): return h / np.sqrt(Rs ** 2 + c * dx ** 2)   # 1/Mpc
KP = dict(gg=kp_of(1 / 12), kg=kp_of(13 / 84), kk=kp_of(19 / 84))
pw = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pw1024.npy'))[LMIN:LMAX + 1]   # = hp.pixwin(1024)
def lim(WA, WB, kp): return sw @ ((C_KMS / Hz)[:, None] * WA[:, None] * WB[:, None] / chi[:, None] ** 2 * P * np.exp(-(k / kp) ** 2))
CKK = lim(Wk, Wk, KP['kk'])

def load(box, cat):
    out = {}
    for xx, f in [('gg', f'cl_gg_{cat}_bin0'), ('kg', f'cl_kg_{cat}_bin0'), ('kk', 'cl_kk')]:
        d = np.load(D + f'{box}/cls/{f}.npz'); out[xx] = d['cls'][LMIN:LMAX + 1]
        if 'shot_noise' in d.files: out['sn'] = float(d['shot_noise'])
    n = np.load(D + f'{box}/cls/nz_{cat}_bin0.npz'); out['zmean'] = float(n['z_mean'])
    zn, nn = n['z_nz'], n['nz']; out['nz'] = PchipInterpolator(zn, nn / simpson(nn, x=zn), extrapolate=False)
    return out
def Wg_of(nz, cat):
    zb, b = np.loadtxt(D + f'colore_setup/inputs/bias/bz_{cat}.txt', unpack=True)
    bz = np.nan_to_num(PchipInterpolator(zb, b, extrapolate=False)(z)); n = np.nan_to_num(nz(z))
    return bz * n * Hz / C_KMS
def theory(nz, cat, kp=KP):
    Wg = Wg_of(nz, cat); return lim(Wg, Wg, kp['gg']), lim(Wk, Wg, kp['kg'])
def lmax_of(rmin, zmean): return int(2 * np.pi / rmin * chi_of(zmean) * h - 0.5)
def mask_of(lmax): return (lb >= LMIN) & (ells_hi <= lmax)

def cov(A, Asn, cgg, ckg, ckk, sn, m, N=1):
    G = (A ** 2 * cgg + Asn * sn + sn) * pw ** 2; K = A * ckg * pw; L = 2 * ells + 1
    vgg = 2 * G ** 2 / L; vkg = (K ** 2 + ckk * G) / L; xc = 2 * K * G / L
    M = [(B @ np.diag(v) @ B.T)[np.ix_(m, m)] / N for v in (vgg, vkg, xc)]
    return np.block([[M[0], M[2]], [M[2].T, M[1]]])
def fit_joint(dgg, dkg, cgg, ckg, sn, m, N=1, ckk=None):
    ckk = CKK if ckk is None else ckk
    d = np.concatenate([(B @ dgg)[m], (B @ dkg)[m]])
    Bg = B * pw[None, :] ** 2; Bk = B * pw[None, :]
    tg = (Bg @ cgg)[m]; tk = (Bk @ ckg)[m]; ts = (Bg @ np.full(len(ells), sn))[m]
    mod = lambda p: np.concatenate([p[0] ** 2 * (p[1] / s0) ** 2 * tg + p[2] * ts, p[0] * (p[1] / s0) ** 2 * tk])
    A0, N0 = 1., 0.
    for it in range(2):
        Ci = np.linalg.inv(cov(A0, N0, cgg, ckg, ckk, sn, m, N))
        f = lambda p: (d - mod(p)) @ Ci @ (d - mod(p))
        r = minimize(f, [1, s0, 0], method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=40000))
        r = minimize(f, r.x, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=40000))
        A0, N0 = r.x[0], r.x[2]
    J = np.array([(mod(r.x + e) - mod(r.x - e)) / (2 * e.max()) for e in np.diag([1e-5, 1e-5, 1e-3])]).T
    err = np.sqrt(np.diag(np.linalg.inv(J.T @ Ci @ J)))
    return r.x, err, r.fun, len(d) - 3, f
def fit_kk(dkk, m, N=1, ckk=None):
    ckk = CKK if ckk is None else ckk
    t = (B @ ckk)[m]; d = (B @ dkk)[m]; Ci = np.linalg.inv((B @ np.diag(2 * ckk ** 2 / (2 * ells + 1)) @ B.T / N)[np.ix_(m, m)])
    a = (t @ Ci @ d) / (t @ Ci @ t); ea = 1 / np.sqrt(t @ Ci @ t)
    return np.sqrt(a), .5 * ea / np.sqrt(a), (d - a * t) @ Ci @ (d - a * t), (d - t) @ Ci @ (d - t), m.sum() - 1
