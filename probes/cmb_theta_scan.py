#!/usr/bin/env python3
"""
CMB theta-scan (compressed Planck 2018 + ACT DR6 likelihood)
============================================================

Port of archive v10.6 section 4a onto the shared bimetric core.
Data: Bansal & Huterer 2025 (arXiv:2502.07185), vector (R, l_a, omega_b)
with Gaussian covariance.

Canonical v10.6 results (rev 43, verified ✓✓ by independent DeepSeek code;
also reproduced verbatim by Qwen (rev 44) and GLM (rev 46) codes):
  LCDM: chi2 = 1.41 (h = 0.72085)
  theta = 6.1 deg:    0.04 (h = 0.72425, Om = 0.30312)
          10.0 deg:   9.34 (h = 0.72605, Om = 0.30904)
          12.2 deg:  31.62 (h = 0.72749, Om = 0.31375)
          15.0 deg:  93.65 (h = 0.72984, Om = 0.32129)
          18.3 deg: 244.71 (h = 0.73343, Om = 0.33260)
          20.0 deg: 370.33 (h = 0.73568, Om = 0.33955)
          25.0 deg:1029.42 (h = 0.74420, Om = 0.36514)
          30.0 deg:2368.93 (h = 0.75643, Om = 0.40006)

Usage:  python3 -m probes.cmb_theta_scan
"""

import sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, '.')
from bimetric import (Z_STAR, OMEGA_NU, OMEGA_R_H2,
                      solve_Om0, r_s, D_M_bimetric)

# --- Compressed CMB data (Bansal & Huterer 2025) ---
DATA_VEC = np.array([1.7504, 301.77, 0.022371])          # (R, l_a, omega_b)
C_CMB = 1e-8 * np.array([[1559.83,  -1325.41,   -36.45],
                         [-1325.41, 714691.80,  269.77],
                         [-36.45,    269.77,     2.10]])
C_INV = np.linalg.inv(C_CMB)


def chi2_cmb(vec):
    d = vec - DATA_VEC
    return float(d @ C_INV @ d)


def compute_observables(h, wb, Om0, OL, mFP2, sin2, cos2, alpha, beta, Or,
                        wnu_h):
    """(R, l_a, omega_b); R prefactor per Bansal & Huterer Eq. 4."""
    DM = D_M_bimetric(Z_STAR, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or,
                      wnu_h)
    R = np.sqrt(Om0 + wnu_h) * DM
    rs = r_s(Z_STAR, wb, h, Om0, OL, wnu_h)
    return np.array([R, np.pi * DM / rs, wb])


def minimize_chi2(Om0, OL, mFP2, sin2, cos2, alpha, beta, Or, wnu_h):
    """Nelder-Mead over (h, omega_b); coarse grid + polish.

    FIX-19/FIX-23: h in [0.55, 0.88], wb in [0.015, 0.025] (theta >= 25 deg
    hits h = 0.75 and wb = 0.02 boundaries of the old grid).
    """
    def objective(p):
        h, wb = p
        if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
            return 1e10
        return chi2_cmb(compute_observables(h, wb, Om0, OL, mFP2, sin2, cos2,
                                            alpha, beta, Or, wnu_h))

    best, best_init = 1e10, None
    for h in np.linspace(0.60, 0.85, 16):
        for wb in np.linspace(0.015, 0.025, 12):
            c2 = objective([h, wb])
            if c2 < best:
                best, best_init = c2, [h, wb]
    res = minimize(objective, best_init, method='Nelder-Mead',
                   options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 5000})
    h_best, wb_best = res.x
    vec = compute_observables(h_best, wb_best, Om0, OL, mFP2, sin2, cos2,
                              alpha, beta, Or, wnu_h)
    return chi2_cmb(vec), h_best, wb_best, vec


def run_point(theta_deg, alpha=0.0, beta=3.0, mFP_H0=25.0, OL=0.699):
    """One scan point with the 3-iteration h-consistency loop (archive v10.6)."""
    th = np.radians(theta_deg)
    sin2, cos2 = np.sin(th) ** 2, np.cos(th) ** 2
    mFP2 = mFP_H0 ** 2
    Or = OMEGA_R_H2 / 0.678 ** 2
    wnu_h = OMEGA_NU / 0.678 ** 2
    Om0 = solve_Om0(mFP2, sin2, cos2, alpha, beta, OL, Or, wnu_h)
    for _ in range(3):
        chi2_val, h_best, wb_best, vec = minimize_chi2(
            Om0, OL, mFP2, sin2, cos2, alpha, beta, Or, wnu_h)
        Or = OMEGA_R_H2 / h_best ** 2
        wnu_h = OMEGA_NU / h_best ** 2
        Om0 = solve_Om0(mFP2, sin2, cos2, alpha, beta, OL, Or, wnu_h)
    return dict(theta=theta_deg, chi2=chi2_val, h=h_best, wb=wb_best,
                Om=Om0, R=vec[0], la=vec[1])


def lcdm_sanity():
    """LCDM control (archive v10.6 convention, target chi2 = 1.41)."""
    Or = OMEGA_R_H2 / 0.678 ** 2
    wnu_h = OMEGA_NU / 0.678 ** 2
    Om0 = 0.309 - wnu_h
    chi2_val, h_best, wb_best, vec = minimize_chi2(
        Om0, 0.691, 0.0, 0.0, 1.0, 0.0, 0.0, Or, wnu_h)
    return dict(theta=0.0, chi2=chi2_val, h=h_best, wb=wb_best, Om=Om0,
                R=vec[0], la=vec[1])


if __name__ == '__main__':
    print("CMB THETA-SCAN (bimetric E^2, core of archive v10.6)")
    print("=" * 64)
    lc = lcdm_sanity()
    print(f"LCDM: chi2={lc['chi2']:.2f}, h={lc['h']:.5f}, wb={lc['wb']:.7f}\n")
    print(f"{'theta':>6} {'chi2':>10} {'h':>8} {'Omega_m':>8} "
          f"{'R':>8} {'l_a':>8}")
    print("-" * 64)
    for t in [6.1, 10.0, 12.2, 15.0, 18.3, 20.0, 25.0, 30.0]:
        r = run_point(t)
        print(f"{t:6.1f} {r['chi2']:10.2f} {r['h']:8.5f} {r['Om']:8.5f} "
              f"{r['R']:8.4f} {r['la']:8.3f}")
