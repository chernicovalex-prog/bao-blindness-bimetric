"""Independent high-precision implementation of the bimetric background.

Self-contained: derives the quartic from the B-parameters (H&M 2021a, Eqs. 3.4
and 5.5) and solves it with mpmath.polyroots at 40-digit precision -- no
import from core/scalarsec/selfsim. Cross-checked against the canonical
pipeline (Model.E2_bim) over a redshift grid; doubles as the large-m_FP
mpmath path (D-07 lesson: double precision loses the root at large masses).
"""
import numpy as np
from mpmath import mp, mpf, matrix, polyroots, sqrt as msqrt, findroot

mp.dps = 40


def B_params(theta_deg, mFP_H0, alpha, beta, Omega_L):
    t = np.tan(np.radians(theta_deg)) ** 2
    s2 = np.sin(np.radians(theta_deg)) ** 2
    m2 = mFP_H0 ** 2
    s2m2 = s2 * m2
    B0 = 3 * Omega_L - s2m2 * (3 + 3 * alpha + beta)
    B1 = s2m2 * (1 + 2 * alpha + beta)
    B2 = -s2m2 * (alpha + beta)
    B3 = s2m2 * beta
    B4 = 3 * t * Omega_L + s2m2 * (-1 + alpha - beta)
    return t, s2, m2, (B0, B1, B2, B3, B4)


def coeffs(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r):
    t, s2, m2, (B0, B1, B2, B3, B4) = B_params(theta_deg, mFP_H0, alpha, beta, Omega_L)
    Omz = Om * (1 + z) ** 3 + Omega_r * (1 + z) ** 4
    return [t * s2 * m2 * beta,
            3 * t * B2 - B4,
            3 * (t * B1 - B3),
            t * (B0 + 3 * Omz) - 3 * B2,
            -B1]


def y_root(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r):
    c = coeffs(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r)
    cm = matrix([mpf(repr(float(x))) for x in c])
    roots = polyroots(cm, maxsteps=200, error=False)
    reals = [r.real for r in roots if abs(float(r.imag)) < 1e-20 and float(r.real) > 0]
    if not reals:
        raise ValueError("no positive real root (theta=%g, m=%g, z=%g)" % (theta_deg, mFP_H0, z))
    return min(reals, key=float)


def E2_of_y(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r, y):
    s2 = np.sin(np.radians(theta_deg)) ** 2
    m2 = mFP_H0 ** 2
    d = 1 - float(y)
    return Om * (1 + z) ** 3 + Omega_r * (1 + z) ** 4 + Omega_L \
        - s2 * m2 * d * (1 + alpha * d + beta / 3.0 * d * d)


def E2_indep(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r=9.08e-5):
    y = y_root(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r)
    return float(E2_of_y(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, z, Omega_r, y))


def closure_Om(theta_deg, mFP_H0, alpha, beta, Omega_L, Omega_r=9.08e-5):
    """Self-consistent Om from E2(0) = 1 by bisection (mpmath findroot)."""
    f = lambda Om: E2_indep(theta_deg, mFP_H0, alpha, beta, Omega_L, Om, 0.0, Omega_r) - 1.0
    lo, hi = 0.05, 0.95
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0, "closure root not bracketed"
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-15:
            break
    return 0.5 * (lo + hi)


def crosscheck(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0, Omega_L=0.699,
               Omega_r=9.08e-5, n_z=25):
    """Independent vs canonical E^2(z); returns (max rel dev, Om_indep, Om_canon)."""
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from core import Model
    Om_i = closure_Om(theta_deg, mFP_H0, alpha, beta, Omega_L, Omega_r)
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om_c = m.solve_Om0()
    zg = np.geomspace(1e-3, 1090.0, n_z)
    dev = 0.0
    for z in zg:
        e_i = E2_indep(theta_deg, mFP_H0, alpha, beta, Omega_L, Om_c, float(z), Omega_r)
        e_c = m.E2_bim(float(z), Om_c)
        dev = max(dev, abs(e_i / e_c - 1.0))
    return dict(max_rel_dev=float(dev), Om_indep=Om_i, Om_canon=float(Om_c),
                theta_deg=theta_deg, mFP_H0=mFP_H0)


if __name__ == "__main__":
    import json
    for kw in (dict(theta_deg=6.1, mFP_H0=25.0),
               dict(theta_deg=30.0, mFP_H0=200.0),
               dict(theta_deg=6.1, mFP_H0=500.0)):
        r = crosscheck(**kw)
        print(json.dumps(r))
