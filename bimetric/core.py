"""
Bimetric cosmology core
=======================

Single source of truth for the background physics of the project
"BAO Blindness to Bimetric Gravity" (Chernikov 2026).

Ported verbatim from Bimetric Code Archive v10.6 (sections 1 and 4a),
which contains all 24 fixes (FIX-1 ... FIX-24) found during the
multi-AI cross-review (revisions 27-43 of the project journal).
FIX references are kept in the docstrings.

Physics references:
  - Hogas & Mortsell 2025 (arXiv:2507.03743): Friedmann eq. (A11),
    Omega_DE (A12), quartic for y = a_f/a_g (A13).
  - Hogas & Mortsell 2021 (arXiv:2101.08794): B-parameters (3.4a-e).
  - Bansal & Huterer 2025 (arXiv:2502.07185): compressed CMB likelihood,
    R prefactor (Eq. 4), massive-nu convention (App. C).

License: MIT
"""

import sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

# ==========================================================================
# Physical constants (archive v10.6 section 4a)
# ==========================================================================
Z_STAR = 1090.0              # recombination redshift (compressed CMB)
T_CMB = 2.725                # CMB temperature [K]
OMEGA_GAMMA_H2 = 2.469e-5    # photon density omega_gamma
N_EFF = 3.046
M_NU_SUM = 0.06              # sum of neutrino masses [eV]
OMEGA_NU = M_NU_SUM / 93.14  # ~6.44e-4
OMEGA_R_H2 = OMEGA_GAMMA_H2 * (1 + 7.0 / 8.0 * N_EFF * (4.0 / 11.0) ** (4.0 / 3.0))

Z_NR = 131.0                 # 1+z at which m_nu = 0.06 eV becomes non-relativistic
                             # (FIX-23, Bansal & Huterer App. C convention)

C_OVER_H0 = 2997.92458       # c/H0 in Mpc for h = 1


# ==========================================================================
# Massive neutrino density (FIX-23)
# ==========================================================================
def nu_density(z, wnu_h):
    """Massive-nu density in units of rho_crit,0.

    Matter-like (w^3) below 1+z = Z_NR, radiation-like (w^4/Z_NR) above;
    continuous at the switch. (Bansal & Huterer 2025, App. C.)
    """
    zp1 = 1.0 + z
    if zp1 <= Z_NR:
        return wnu_h * zp1 ** 3
    return wnu_h * zp1 ** 4 / Z_NR


# ==========================================================================
# Quartic for y = a_f / a_g  (Hogas & Mortsell 2025, Eq. A13)
# ==========================================================================
def quartic_coeffs(z, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or):
    """Coefficients of a0 + a1 y + a2 y^2 + a3 y^3 + a4 y^4 = 0."""
    Ot = Om0 * (1 + z) ** 3 + Or * (1 + z) ** 4
    a0 = -1.0 / 3.0 * cos2 * mFP2 * (1 + 2 * alpha + beta)
    a1 = (Ot + OL
          + mFP2 * (cos2 * (alpha + beta) - sin2 * (1 + alpha + beta / 3.0)))
    a2 = mFP2 * (-cos2 * beta + sin2 * (1 + 2 * alpha + beta))
    a3 = -(OL
           + 1.0 / 3.0 * mFP2 * (cos2 * (-1 + alpha - beta)
                                 + 3 * sin2 * (alpha + beta)))
    a4 = 1.0 / 3.0 * sin2 * mFP2 * beta
    return a0, a1, a2, a3, a4


def solve_y(z, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or):
    """Physical branch y(z): Brent root on [1e-6, 1.0].

    FIX-4/FIX-10: brentq range [1e-6, 1.0] (a narrow [0.5, 1.0] silently
    reverted to GR for large theta).
    FIX-21: mFP2 == 0 -> y = 1 immediately (GR limit, no warning).
    FIX-13/FIX-17: fallback logs a warning to stderr, never silent.
    """
    if mFP2 == 0:
        return 1.0
    a0, a1, a2, a3, a4 = quartic_coeffs(z, Om0, mFP2, sin2, cos2,
                                        alpha, beta, OL, Or)

    def f(y):
        return a0 + a1 * y + a2 * y ** 2 + a3 * y ** 3 + a4 * y ** 4

    try:
        return brentq(f, 1e-6, 1.0, xtol=1e-14)
    except ValueError:
        print(f"WARNING [core]: brentq fallback (y=1.0) at z={z}",
              file=sys.stderr)
        return 1.0


# ==========================================================================
# Bimetric E^2(z)  (A11 + A12, with massive-nu term, FIX-23)
# ==========================================================================
def omega_DE(y, mFP2, sin2, alpha, beta, OL):
    """Effective dark-energy density Omega_DE(z) from the quartic root."""
    dy = 1.0 - y
    return OL - sin2 * mFP2 * dy * (1 + alpha * dy + beta / 3.0 * dy ** 2)


def E2_bimetric(z, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or, wnu_h=0.0):
    """E^2(z) = H^2(z)/H_0^2 for the bimetric model (full, not LCDM-approx)."""
    y = solve_y(z, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or)
    return (Om0 * (1 + z) ** 3 + Or * (1 + z) ** 4
            + omega_DE(y, mFP2, sin2, alpha, beta, OL)
            + nu_density(z, wnu_h))


def E2_lcdm(z, Om0, Or=5.38e-5):
    return Om0 * (1 + z) ** 3 + Or * (1 + z) ** 4 + (1 - Om0 - Or)


# ==========================================================================
# Self-consistent Omega_m from E^2(0) = 1 (closure)
# ==========================================================================
def solve_Om0(mFP2, sin2, cos2, alpha, beta, OL, Or, wnu_h=0.0):
    """Omega_m from E^2(0) = 1; Brent on [0.01, 0.99]."""
    def eq(Om):
        y = solve_y(0.0, Om, mFP2, sin2, cos2, alpha, beta, OL, Or)
        return 1.0 - Or - omega_DE(y, mFP2, sin2, alpha, beta, OL) - Om - wnu_h
    return brentq(eq, 0.01, 0.99, xtol=1e-14)


# ==========================================================================
# Sound speed and sound horizon (FIX-1, FIX-20, FIX-24)
# ==========================================================================
def sound_speed(z, wb):
    """c_s in units of c.

    FIX-24: R_b = 31500 * wb * (T_CMB/2.7)**(-4) / (1+z).
    The exponent sign is negative (omega_gamma ~ T^4). The wrong sign
    biased r_s by 0.7% and h by ~0.011 (found by blind cross-check, rev 43).
    """
    R_b = 31500 * wb * (T_CMB / 2.7) ** (-4) / (1 + z)
    return 1.0 / np.sqrt(3 * (1 + R_b))


def r_s(z_rec, wb, h, Om0, OL_de, wnu_h=0.0):
    """Sound horizon at z_rec, in units of c/H0.

    Integrated in the scale factor a with the 1/a^2 Jacobian (FIX-1).
    E2_0 = 1.0 by construction (FIX-20): solve_Om0 guarantees E^2(0)=1
    through ODE(0), which differs from the input OL_de for theta > 0.
    Above z_rec the early universe is Vainshtein-screened, so the
    LCDM-like form with OL_de is used (archive v10.6, section 4a).
    """
    def integrand(a):
        z = 1.0 / a - 1.0
        cs = sound_speed(z, wb)
        E2 = (Om0 * (1 + z) ** 3 + (OMEGA_R_H2 / h ** 2) * (1 + z) ** 4
              + OL_de + nu_density(z, wnu_h))
        E2_0 = 1.0  # FIX-20
        return cs / (a * a * np.sqrt(E2 / E2_0))

    a_rec = 1.0 / (1.0 + z_rec)
    result, _ = quad(integrand, 1e-8, a_rec, limit=200)
    return result


# ==========================================================================
# Comoving distance with full bimetric E^2 (FIX-16: adaptive quad)
# ==========================================================================
def D_M_bimetric(z_rec, Om0, mFP2, sin2, cos2, alpha, beta, OL, Or, wnu_h=0.0):
    """D_M = integral_0^z_rec dz / E(z), units of c/H0.

    FIX-16: scipy.integrate.quad (adaptive). A 3000-point trapezoid missed
    the sharp bimetric transition at z < 10 (0.16% in D_M, up to 15% in chi2).
    """
    def integrand(z):
        return 1.0 / np.sqrt(E2_bimetric(z, Om0, mFP2, sin2, cos2,
                                         alpha, beta, OL, Or, wnu_h))
    result, _ = quad(integrand, 0, z_rec, limit=200)
    return result


# ==========================================================================
# B-parameters (Hogas & Mortsell 2021, eqs. 3.4a-e / 3.3)
# ==========================================================================
def compute_B_params(theta_deg, mFP_H0, alpha, beta, Omega_Lambda):
    """Forward conversion (physical -> B-parameters), eq. 3.4a-e."""
    th = np.radians(theta_deg)
    s2 = np.sin(th) ** 2
    m2 = mFP_H0 ** 2
    return {
        'B0': 3 * Omega_Lambda - s2 * m2 * (3 + 3 * alpha + beta),
        'B1': s2 * m2 * (1 + 2 * alpha + beta),
        'B2': -s2 * m2 * (alpha + beta),
        'B3': s2 * m2 * beta,
        'B4': 3 * np.tan(th) ** 2 * Omega_Lambda + s2 * m2 * (-1 + alpha - beta),
    }


def invert_B_params(B, theta_deg):
    """Inverse conversion (B-parameters -> physical), eq. 3.3."""
    th = np.radians(theta_deg)
    s2 = np.sin(th) ** 2
    den = B['B1'] + 2 * B['B2'] + B['B3']
    return {
        'mFP2': den / s2,
        'Omega_Lambda': B['B0'] / 3 + B['B1'] + B['B2'] + B['B3'] / 3,
        'alpha': -(B['B2'] + B['B3']) / den,
        'beta': B['B3'] / den,
    }
