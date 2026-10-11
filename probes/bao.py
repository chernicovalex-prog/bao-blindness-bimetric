#!/usr/bin/env python3
"""
BAO probe: DESI DR2 (13 points) vs bimetric baseline
=====================================================

Port of archive v10.6 sections 1-2 onto the shared bimetric core.

Two covariance options:
  --cov full   : Cobaya 13x13 matrix (desi_gaussian_bao_ALL_GCcomb_cov.txt,
                 archive section 2)  -> canonical Delta chi2 = 0.178
  --cov table4 : block-diagonal from DESI DR2 Table IV sigmas/correlations
                 (archive section 1) -> Delta chi2 = 0.195 (rounded sigmas)

The free parameter f = c/(H0 r_d) is marginalised by direct minimisation
(1-D bounded), per probe.

Canonical results (v10.6, verified by DeepSeek independent code, rev 27):
  bimetric baseline theta=6.1, m_FP/H0=25, alpha=0, beta=3, OL=0.699:
  y(0) = 0.999512, Omega_m = 0.30435-0.30439, E^2(0) = 1 + 5e-15
  Delta chi2(BAO) = 0.178 (full cov) / 0.195 (Table IV)

Usage:  python3 -m probes.bao [--cov full|table4]
"""

import sys
import argparse
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, '.')
from bimetric import solve_y, omega_DE, E2_bimetric, compute_B_params, invert_B_params

# --- DESI DR2 BAO data vector (Table IV, arXiv:2503.14738) ---
# (z, type, value, sigma, correlation with the partner measurement)
DESI_DATA = [
    (0.295, 'DV', 7.942, 0.075, None),
    (0.510, 'DM', 13.588, 0.167, -0.459),
    (0.510, 'DH', 21.863, 0.425, -0.459),
    (0.706, 'DM', 17.351, 0.177, -0.404),
    (0.706, 'DH', 19.455, 0.330, -0.404),
    (0.934, 'DM', 21.576, 0.152, -0.416),
    (0.934, 'DH', 17.641, 0.193, -0.416),
    (1.321, 'DM', 27.601, 0.318, -0.434),
    (1.321, 'DH', 14.176, 0.221, -0.434),
    (1.484, 'DM', 30.512, 0.760, -0.500),
    (1.484, 'DH', 12.817, 0.516, -0.500),
    (2.330, 'DM', 38.988, 0.531, -0.431),
    (2.330, 'DH', 8.632, 0.101, -0.431),
]
DATA_VALS = np.array([d[2] for d in DESI_DATA])

# --- Full 13x13 covariance (Cobaya bao_data desi_bao_dr2, archive section 2) ---
# Order matches DESI_DATA except Ly-alpha rows: in the Cobaya file row 12 is
# DH/rd (8.6315) and row 13 is DM/rd (38.9890); DESI_DATA lists DM first.
# The matrix below is therefore used with a permuted data vector (see build).
COV_FULL_RAW = np.array([
    [5.78998687e-03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 2.83473742e-02, -3.26062007e-02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, -3.26062007e-02, 1.83928040e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 3.23752442e-02, -2.37445646e-02, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -2.37445646e-02, 1.11469198e-01, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 2.61732816e-02, -1.12938006e-02, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1.12938006e-02, 4.04183878e-02, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1.05336516e-01, -2.90308418e-02, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -2.90308418e-02, 5.04233092e-02, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 5.83020277e-01, -1.95215562e-01, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -1.95215562e-01, 2.68336193e-01, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.02136194e-02, -2.31395216e-02],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.31395216e-02, 2.82685779e-01],
], dtype=float)


def build_covariance(mode):
    """Return (C_inv, permutation) for the chosen covariance convention."""
    if mode == 'table4':
        n = len(DESI_DATA)
        C = np.zeros((n, n))
        C[0, 0] = DESI_DATA[0][3] ** 2
        for i in range(1, n, 2):
            r = DESI_DATA[i][4]
            C[i, i] = DESI_DATA[i][3] ** 2
            C[i + 1, i + 1] = DESI_DATA[i + 1][3] ** 2
            C[i, i + 1] = r * DESI_DATA[i][3] * DESI_DATA[i + 1][3]
            C[i + 1, i] = C[i, i + 1]
        return np.linalg.inv(C), np.arange(n)
    # full Cobaya matrix: swap the two Ly-alpha entries (rows/cols 11,12)
    perm = np.arange(13)
    perm[11], perm[12] = 12, 11
    C = COV_FULL_RAW[np.ix_(perm, perm)]
    return np.linalg.inv(C), np.arange(13)


def _integrate(y, x):
    return np.sum(0.5 * (y[:-1] + y[1:]) * (x[1:] - x[:-1]))


def comoving_distance(z, mFP2, sin2, cos2, alpha, beta, OL, Or, Om0, Nz=500):
    zs = np.linspace(0, z, Nz)
    E2s = np.array([E2_bimetric(zi, Om0, mFP2, sin2, cos2, alpha, beta,
                                OL, Or) for zi in zs])
    return _integrate(1.0 / np.sqrt(E2s), zs)


def bao_prediction(f, z, obs_type, mFP2, sin2, cos2, alpha, beta, OL, Or, Om0):
    """D_M/r_d, D_H/r_d or D_V/r_d (D_V per Eisenstein & Hu, used by DESI)."""
    if obs_type == 'DH':
        return f / np.sqrt(E2_bimetric(z, Om0, mFP2, sin2, cos2, alpha, beta,
                                       OL, Or))
    if obs_type == 'DM':
        return f * comoving_distance(z, mFP2, sin2, cos2, alpha, beta,
                                     OL, Or, Om0)
    if obs_type == 'DV':
        E2 = E2_bimetric(z, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or)
        DM = comoving_distance(z, mFP2, sin2, cos2, alpha, beta, OL, Or, Om0)
        return f * (z * DM ** 2 / np.sqrt(E2)) ** (1.0 / 3.0)
    raise ValueError(obs_type)


def chi2_bimetric(f, C_inv, mFP2, sin2, cos2, alpha, beta, OL, Or, Om0):
    model = np.array([bao_prediction(f, z, ot, mFP2, sin2, cos2, alpha, beta,
                                     OL, Or, Om0)
                      for z, ot, v, e, r in DESI_DATA])
    d = DATA_VALS - model
    return float(d @ C_inv @ d)


def chi2_lcdm(f, C_inv, Om0, Or=5.38e-5):
    model = np.zeros(len(DESI_DATA))
    for i, (z, ot, v, e, r) in enumerate(DESI_DATA):
        E2 = Om0 * (1 + z) ** 3 + Or * (1 + z) ** 4 + (1 - Om0 - Or)
        if ot == 'DH':
            model[i] = f / np.sqrt(E2)
        else:
            zs = np.linspace(0, z, 500)
            E2s = Om0 * (1 + zs) ** 3 + Or * (1 + zs) ** 4 + (1 - Om0 - Or)
            DM = _integrate(1.0 / np.sqrt(E2s), zs)
            if ot == 'DM':
                model[i] = f * DM
            else:  # DV
                model[i] = f * (z * DM ** 2 / np.sqrt(E2)) ** (1.0 / 3.0)
    d = DATA_VALS - model
    return float(d @ C_inv @ d)


def run(cov_mode='full', theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
        Omega_Lambda=0.699, Omega_r=5.38e-5):
    C_inv, _ = build_covariance(cov_mode)
    th = np.radians(theta_deg)
    sin2, cos2 = np.sin(th) ** 2, np.cos(th) ** 2
    mFP2 = mFP_H0 ** 2

    # Self-consistent Omega_m from the quartic at z=0 (archive section 1)
    y0 = solve_y(0.0, 0.3044, mFP2, sin2, cos2, alpha, beta, Omega_Lambda,
                 Omega_r)
    Om0 = 1.0 - Omega_r - omega_DE(y0, mFP2, sin2, alpha, beta, Omega_Lambda)
    E2_0 = E2_bimetric(0.0, Om0, mFP2, sin2, cos2, alpha, beta,
                       Omega_Lambda, Omega_r)

    f_b, chi2_b = minimize_scalar(
        lambda f: chi2_bimetric(f, C_inv, mFP2, sin2, cos2, alpha, beta,
                                Omega_Lambda, Omega_r, Om0),
        bounds=(20, 80), method='bounded').x, None
    res_b = minimize_scalar(
        lambda f: chi2_bimetric(f, C_inv, mFP2, sin2, cos2, alpha, beta,
                                Omega_Lambda, Omega_r, Om0),
        bounds=(20, 80), method='bounded')
    res_l = minimize_scalar(lambda f: chi2_lcdm(f, C_inv, 0.2972),
                            bounds=(20, 80), method='bounded')
    return dict(y0=y0, Om0=Om0, E2_0=E2_0,
                f_bim=res_b.x, chi2_bim=res_b.fun,
                f_lcdm=res_l.x, chi2_lcdm=res_l.fun,
                dchi2=res_b.fun - res_l.fun, cov=cov_mode)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--cov', choices=['full', 'table4'], default='full')
    args = ap.parse_args()
    r = run(args.cov)
    B = compute_B_params(6.1, 25.0, 0.0, 3.0, 0.699)
    print("BAO: DESI DR2 vs bimetric baseline (cov = %s)" % r['cov'])
    print("=" * 60)
    print(f"y(0) = {r['y0']:.6f}   Omega_m = {r['Om0']:.8f}   "
          f"E^2(0) = {r['E2_0']:.15f}")
    print(f"LCDM:     chi2 = {r['chi2_lcdm']:.4f} (f = {r['f_lcdm']:.2f})")
    print(f"Bimetric: chi2 = {r['chi2_bim']:.4f} (f = {r['f_bim']:.2f})")
    print(f"Delta chi2 = {r['dchi2']:.4f}")
    print("B-params: " + "  ".join(f"{k}={v:.3f}" for k, v in B.items())
          + f"  B1^2={B['B1']**2:.2f}")
    Bi = invert_B_params(B, 6.1)
    print(f"Round-trip: mFP2 {625.0:.6f}->{Bi['mFP2']:.6f}  "
          f"OL 0.699->{Bi['Omega_Lambda']:.6f}  "
          f"alpha 0->{Bi['alpha']:.2e}  beta 3->{Bi['beta']:.6f}")
