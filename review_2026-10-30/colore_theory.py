"""
CoLoRe-consistent theory for C_ell^gg, C_ell^kg, C_ell^kk (lognormal sims, _BIAS_MODEL_3 + threshold).

Drop into utils/ and use from the notebook. Everything follows CoLoRe master source:
  * cosmo.c      : input P(k) re-splined and renormalised to param-file sigma_8 (top-hat R=8 Mpc/h);
                   flat LCDM background (omega_M, omega_L, w), no radiation, no neutrinos.
  * fourier.c    : field smoothed by exp(-R^2 k^2 / 2)  ->  power x exp(-R^2 k^2).
  * density.c    : delta_LN = exp(D delta_G - D^2 sigma2_gauss / 2) - 1,
                   sigma2_gauss computed with R_eff^2 = R^2 + (0.45 dx)^2.
  * common.h/srcs.c (_BIAS_MODEL_3): lambda ∝ max(1 + b delta_LN, 0), = 0 if delta_LN < t(z).
  * kappa.c      : kappa from the *linear Gaussian* potential  ->  C^kk and C^kg are linear.

Two-point of the galaxy field (Mehler / Hermite expansion, exact for a pointwise function of a Gaussian):
  xi_gg(r) = sum_{n>=1} (c_n/c_0)^2 rho^n / n!,  c_n = <f(y) He_n(y)>,  rho = xi_G / sigma_G^2
  xi_gG(r) = b_x xi_G,  b_x = c_1 / (c_0 sigma_G)          (Price's theorem: kg is exactly linear)
Units: internal 3D quantities in h/Mpc, (Mpc/h)^3; returned P(k) in 1/Mpc, Mpc^3 to match the notebook.
Requires: numpy, scipy, mcfit (pip install mcfit).
"""
import numpy as np
import mcfit
from math import factorial
from scipy.integrate import simpson, quad, cumulative_trapezoid
from scipy.interpolate import interp1d

C_KMS = 299792.458


class CoLoReCosmo:
    """Background, growth and normalised linear P(k) exactly as CoLoRe uses them."""

    def __init__(self, cosmo_dict, pk_filename):
        self.h = cosmo_dict['h']; self.Om = cosmo_dict['omega_M']; self.OL = cosmo_dict['omega_L']
        self.sigma8 = cosmo_dict['sigma_8']
        za = np.concatenate([np.linspace(0, 10, 20001), np.logspace(np.log10(10.0005), np.log10(1200), 2000)])
        E = np.sqrt(self.Om * (1 + za)**3 + self.OL)
        chi = cumulative_trapezoid(C_KMS / (100 * self.h) / E, za, initial=0)      # Mpc
        self._chi = interp1d(za, chi, kind='cubic')
        a_of = lambda z: 1 / (1 + z)
        integrand = lambda a: 1 / (a * np.sqrt(self.Om / a**3 + self.OL))**3
        Dun = lambda z: 2.5 * self.Om * np.sqrt(self.Om / a_of(z)**3 + self.OL) * quad(integrand, 0, a_of(z), epsrel=1e-10)[0]
        zg = np.linspace(0, 10, 1001); Dg = np.array([Dun(z) for z in zg])
        self._D = interp1d(zg, Dg / Dg[0], kind='cubic')
        kh, Ph = np.loadtxt(pk_filename, unpack=True)
        lk = np.log10(kh); lkg = np.linspace(lk[0], lk[-1], len(lk))
        Pg = interp1d(lk, Ph, kind='cubic')(lkg); kk = 10**lkg
        x = 8 * kk; W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
        s2 = simpson(kk**3 * Pg * W**2 / (2 * np.pi**2), x=np.log(kk))
        self.sigma8_file = np.sqrt(s2); self.norm = self.sigma8**2 / s2
        self._lP = interp1d(np.log(kk), np.log(Pg * self.norm), kind='cubic', bounds_error=False, fill_value='extrapolate')

    def chi(self, z): return self._chi(z)                                      # Mpc
    def H(self, z): return 100 * self.h * np.sqrt(self.Om * (1 + z)**3 + self.OL)   # km/s/Mpc
    def D(self, z): return self._D(z)
    def P0_h(self, kh): return np.exp(self._lP(np.log(kh)))                     # (Mpc/h)^3, z=0
    def P(self, z, k): return self.P0_h(k / self.h) / self.h**3 * self.D(z)**2  # Mpc^3, k in 1/Mpc


def colore_cell(cc, n_grid=4096, z_max=3.8):
    """Box side and cell size in Mpc/h, as cosmo.c computes them."""
    rmax = cc.chi(z_max) * cc.h
    L = 2 * rmax * (1 + 2. / n_grid)
    return L, L / n_grid


class GalaxyField:
    """Exact cell-level two-point of CoLoRe _BIAS_MODEL_3 + threshold galaxies."""

    def __init__(self, cc, r_smooth=2.0, n_grid=4096, z_max=3.8, nherm=12):
        self.cc = cc; _, dx = colore_cell(cc, n_grid, z_max); self.dx = dx
        self.reff2 = r_smooth**2 + (0.45 * dx)**2
        kh = np.logspace(-5, 1.6, 8192); P = cc.P0_h(kh)
        self.sig2_0 = np.trapezoid(kh**3 * P * np.exp(-kh**2 * self.reff2) / (2 * np.pi**2), np.log(kh))
        Ps = np.maximum(P * np.exp(-kh**2 * self.reff2), 1e-300)
        self.r, self.xi0 = mcfit.P2xi(kh, lowring=True)(Ps, extrap=False)
        self.nherm = nherm
        self.y = np.linspace(-12, 12, 240001)
        self.w = np.exp(-self.y**2 / 2) / np.sqrt(2 * np.pi) * (self.y[1] - self.y[0])
        self._He = [np.ones_like(self.y), self.y.copy()]
        for n in range(2, nherm): self._He.append(self.y * self._He[-1] - (n - 1) * self._He[-2])

    def coeffs(self, z, b, t):
        s = self.cc.D(z) * np.sqrt(self.sig2_0)
        dLN = np.exp(s * self.y - 0.5 * s**2) - 1
        f = np.where((dLN <= -1) | (dLN < t), 0., np.maximum(1 + b * dLN, 0.))
        return np.array([np.sum(self.w * f * He) for He in self._He]), s

    def b_cross(self, z, b, t):
        c, s = self.coeffs(z, b, t); return c[1] / (c[0] * s)

    def P_excess(self, z, b, t, nmin=2):
        """Non-linear (n>=nmin) galaxy power. Returns k [1/Mpc], P [Mpc^3]."""
        c, s = self.coeffs(z, b, t); rho = self.cc.D(z)**2 * self.xi0 / s**2
        xi = sum((c[n] / c[0])**2 * rho**n / factorial(n) for n in range(nmin, self.nherm))
        k, P = mcfit.xi2P(self.r, lowring=True)(xi, extrap=False)
        return k * self.cc.h, P / self.cc.h**3


def bx_interpolator(gf, bz_input_file, tz_file, zmin, zmax, nz=40):
    """Exact linear (cross) bias b_x(z) of the CoLoRe sample; use it instead of bz_* in W_g."""
    zi, bi = np.loadtxt(bz_input_file, unpack=True); zt, tt = np.loadtxt(tz_file, unpack=True)
    zz = np.linspace(max(zmin, 1e-3), zmax, nz)
    bx = np.array([gf.b_cross(z, np.interp(z, zi, bi), np.interp(z, zt, tt)) for z in zz])
    return interp1d(zz, bx, bounds_error=False, fill_value=np.nan)


def cls_gg_excess(gf, ells, nz_interp, bz_input_file, tz_file, zmin, zmax, nzk=30, nzint=400):
    """Limber C_ell of the Hermite n>=2 term (add to b_x^2 linear C_ell^gg, before the pixel window)."""
    zi, bi = np.loadtxt(bz_input_file, unpack=True); zt, tt = np.loadtxt(tz_file, unpack=True)
    cc = gf.cc; zk = np.linspace(zmin, zmax, nzk)
    tabs = []
    for z in zk:
        k, P = gf.P_excess(z, np.interp(z, zi, bi), np.interp(z, zt, tt))
        tabs.append(interp1d(np.log(k), P, bounds_error=False, fill_value=0.))
    z = np.linspace(zmin, zmax, nzint); chi = cc.chi(z); H = cc.H(z)
    n = np.nan_to_num(nz_interp(z)); W = n * H / C_KMS
    out = np.zeros(len(ells))
    lk = np.log((np.asarray(ells)[None, :] + 0.5) / chi[:, None])
    j = np.clip(np.searchsorted(zk, z) - 1, 0, nzk - 2); wz = (z - zk[j]) / (zk[j + 1] - zk[j])
    P = np.array([(1 - wz[i]) * tabs[j[i]](lk[i]) + wz[i] * tabs[j[i] + 1](lk[i]) for i in range(len(z))])
    return simpson((C_KMS / H)[:, None] * (W**2)[:, None] / chi[:, None]**2 * P, x=z, axis=0)
