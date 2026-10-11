#!/usr/bin/env python3
"""
SNe probe: Pantheon+ (1588 light curves, z > 0.01) vs bimetric baseline
========================================================================

Port of archive v10.6 section 3 (sne_full_cov_v5.py) onto the shared core.

Data files (downloaded by data/download_data.sh):
  data/Pantheon+SH0ES.dat
  data/Pantheon+SH0ES_STAT+SYS.cov

Conventions (verified rev 38 on genuine data):
  - diagonal: m_b_corr_err_DIAG, sigma_int = 0     -> dchi2 = 15.58 (upper bound)
  - full:     STAT+SYS covariance, NO sigma_int    -> dchi2 = 6.50, pure 0.96,
              Om_best = 0.3470, chi2_min = 1433.06, chi2_bim = 1439.56
  - sigma_int must NOT be added on top of full STAT+SYS (double counting)
  - zCMB for the distance integral, zHEL for the (1+z) factor
  - analytic marginalisation over M

Usage:  python3 -m probes.sne [--full]
"""

import sys
import os
import argparse
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid

sys.path.insert(0, '.')
from bimetric import C_OVER_H0, compute_B_params

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
DAT = os.path.join(DATA_DIR, 'Pantheon+SH0ES.dat')
COV = os.path.join(DATA_DIR, 'Pantheon+SH0ES_STAT+SYS.cov')

# Baseline model (archive v10.6)
THETA_DEG, M_FP, ALPHA, BETA = 6.1, 25.0, 0.0, 3.0
OMEGA_LAMBDA, OMEGA_R = 0.699, 9.08e-5


def make_bimetric_mu():
    """Return (Omega_m, mu_function) for the bimetric baseline."""
    th = np.radians(THETA_DEG)
    sin2, tan2, m2 = np.sin(th) ** 2, np.tan(th) ** 2, M_FP ** 2
    B = compute_B_params(THETA_DEG, M_FP, ALPHA, BETA, OMEGA_LAMBDA)

    def solve_y(z, Om):
        Ot = Om * (1 + z) ** 3 + OMEGA_R * (1 + z) ** 4
        a0 = -B['B1'] / (3 * tan2)
        a1 = B['B0'] / 3 - B['B2'] / tan2 + Ot
        a2 = B['B1'] - B['B3'] / tan2
        a3 = B['B2'] - B['B4'] / (3 * tan2)
        a4 = B['B3'] / 3
        return brentq(lambda y: a0 + a1 * y + a2 * y ** 2 + a3 * y ** 3
                      + a4 * y ** 4, 1e-6, 1.0)

    def E2_bim(z, Om, y):
        Ot = Om * (1 + z) ** 3 + OMEGA_R * (1 + z) ** 4
        ODE = OMEGA_LAMBDA - sin2 * m2 * (1 - y) * (1 + ALPHA * (1 - y)
                                                    + BETA / 3 * (1 - y) ** 2)
        return Ot + ODE

    Omega_m = brentq(lambda Om: E2_bim(0, Om, solve_y(0, Om)) - 1.0,
                     0.01, 0.99)

    z_fine = np.linspace(0, 5, 10000)
    y_fine = np.array([solve_y(z, Omega_m) for z in z_fine])
    E_fine = np.sqrt(E2_bim(z_fine, Omega_m, y_fine))
    int_fine = cumulative_trapezoid(1.0 / E_fine, z_fine, initial=0)
    return Omega_m, z_fine, int_fine


def load_data():
    if not (os.path.exists(DAT) and os.path.exists(COV)):
        sys.exit("Data files not found in data/. Run: bash data/download_data.sh")
    import pandas as pd
    # The release file has a '#'-commented header (as on GitHub, rev 51);
    # read positionally to support both the full release and the 10-column
    # reconstruction used in this repository.
    df = pd.read_csv(DAT, sep=r'\s+', comment='#', header=None,
                     usecols=[4, 6, 8, 9],
                     names=['zCMB', 'zHEL', 'm_b_corr', 'm_b_corr_err_DIAG'])
    zCMB = df['zCMB'].values
    zHEL = df['zHEL'].values
    mb = df['m_b_corr'].values                 # FIX-3: m_b_corr, not MU
    mb_err = df['m_b_corr_err_DIAG'].values
    with open(COV) as f:                       # FIX-2: no f.seek(0)
        N_cov = int(f.readline())
        vals = [float(next(f).strip()) for _ in range(N_cov * N_cov)]
    cov_full = np.array(vals).reshape(N_cov, N_cov)
    mask = zCMB > 0.01
    return (zCMB[mask], zHEL[mask], mb[mask], mb_err[mask],
            cov_full[np.ix_(mask, mask)])


def chi2_diag(r, err):
    w = 1.0 / err ** 2
    return np.sum(w * r ** 2) - np.sum(w * r) ** 2 / np.sum(w)


def chi2_full(r, Cinv, S_full):
    Cr = Cinv @ r
    return r @ Cr - np.sum(Cr) ** 2 / S_full


def main(full=False):
    zC, zH, mb_c, err_c, cov_c = load_data()
    N = len(zC)
    Omega_m, z_fine, int_bim = make_bimetric_mu()
    int_bim_zC = np.interp(zC, z_fine, int_bim)
    mu_bim = 5 * np.log10((1 + zH) * C_OVER_H0 * int_bim_zC) + 25

    def mu_lcdm(Om):
        E = np.sqrt(E2_lcdm_arr(z_fine, Om))
        integ = cumulative_trapezoid(1.0 / E, z_fine, initial=0)
        return 5 * np.log10((1 + zH) * C_OVER_H0
                            * np.interp(zC, z_fine, integ)) + 25

    def E2_lcdm_arr(z, Om):
        return Om * (1 + z) ** 3 + OMEGA_R * (1 + z) ** 4 + (1 - Om - OMEGA_R)

    if full:
        Cinv = np.linalg.inv(cov_c)
        S_full = np.sum(Cinv)
        chi2_fn = lambda mu: chi2_full(mb_c - mu, Cinv, S_full)
    else:
        chi2_fn = lambda mu: chi2_diag(mb_c - mu, err_c)

    Om_grid = np.linspace(0.05, 0.60, 500)
    chi2_grid = np.array([chi2_fn(mu_lcdm(Om)) for Om in Om_grid])
    Om_best = Om_grid[np.argmin(chi2_grid)]
    Om_fine = np.linspace(Om_best - 0.01, Om_best + 0.01, 100)
    chi2_fine = np.array([chi2_fn(mu_lcdm(Om)) for Om in Om_fine])
    i = np.argmin(chi2_fine)
    Om_best, chi2_min = Om_fine[i], chi2_fine[i]

    chi2_bim = chi2_fn(mu_bim)
    chi2_fixed = chi2_fn(mu_lcdm(Omega_m))
    print(f"SNe ({'full STAT+SYS' if full else 'diagonal'}, sigma_int=0, "
          f"N={N})")
    print(f"  Omega_m (quartic) = {Omega_m:.5f}")
    print(f"  LCDM: Om_best = {Om_best:.4f}, chi2_min = {chi2_min:.4f}, "
          f"chi2/dof = {chi2_min/(N-1):.3f}")
    print(f"  Bimetric chi2 = {chi2_bim:.4f}")
    print(f"  dchi2 (standard) = {chi2_bim - chi2_min:.2f}")
    print(f"  dchi2 (pure)     = {chi2_bim - chi2_fixed:.2f}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--full', action='store_true',
                    help='use full STAT+SYS covariance (default: diagonal)')
    args = ap.parse_args()
    main(full=args.full)
