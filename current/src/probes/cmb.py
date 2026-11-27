"""Compressed-CMB implementation with one consistent non-DE density.

The quartic, expansion history, closure and direct/spline paths share the
same density convention; the recombination redshift follows Hu--Sugiyama
under standard recombination physics.
           SHA-256 fe874f62fa322919aedc08785aba67c1313b8dd41510baed6327ffc40ecca5d8

Scientific contract implemented (SPEC_SCI_01_v1.md, section 1):

    Omega_ng(z; Omega_m, h) = Omega_m (1+z)^3
                               + Omega_rad(z; h)      [radiation_density, verbatim]
                               + Omega_nu(z; h)        [nu_density,       verbatim]

    c1 = t (B0 + 3 Omega_ng) - 3 B2                     [quartic]
    E^2(z) = Omega_ng(z) + Omega_DE[y(z)]               [expansion]

The SAME Omega_ng is used in the quartic, in E^2, in the z = 0 closure and in
both the direct and the spline execution paths.

Changes with respect to the source file (everything else is byte-identical
in content):

  (1) NEW  omega_ng(z, Om, h)          -- single source of the non-DE density.
  (2) NEW  solve_y_cmb(m, z, Om, h)    -- quartic root with c1 = t(B0+3*omega_ng).
                                           The default Model.solve_y is NOT touched.
  (3) bim_ode / E2_bim_of / solve_Om0_nu / _y_spline / vec_bimetric_m now
      use omega_ng and solve_y_cmb.
  (4) The CMB spline cache key contains h and every parameter entering omega_ng.

Unchanged: DATA_VEC, COV, COV_INV, chi2, z_star_hs, all physical constants
(N_EFF, OMEGA_NU, OMEGA_GAMMA_H2, OMEGA_R_H2, Z_NR), sound_speed, dm_rs, the
LambdaCDM reference (vec_lcdm / fit_lcdm_reference), the optimizer settings,
bounds, starts and the r_map conventions.  No GR fallback is introduced: if
no finite-branch root can be bracketed, an explicit ValueError is raised.
"""
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import interp1d, CubicSpline
import pathlib, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from core import Model, Y_LO


def z_star_hs(wb, wm, h):
    """Recombination redshift, Hu & Sugiyama (1996) fitting formula."""
    wbh2 = wb
    wmh2 = wm
    g1 = 0.0783 * wbh2 ** (-0.238) / (1.0 + 39.5 * wbh2 ** 0.763)
    g2 = 0.560 / (1.0 + 21.1 * wbh2 ** 1.81)
    return 1048.0 * (1.0 + 0.00124 * wbh2 ** (-0.738)) \
        * (1.0 + g1 * wmh2 ** g2)


Z_STAR = 1090.0   # fallback only; all fits use z_star_hs
T_CMB = 2.725
OMEGA_GAMMA_H2 = 2.469e-5
N_EFF = 3.046
OMEGA_NU = 0.06 / 93.14                      # ~ 6.44e-4
OMEGA_R_H2 = OMEGA_GAMMA_H2 * (1 + 7/8 * N_EFF * (4/11) ** (4/3))
# Massive-neutrino scheme (review round 7). The pre-3 scheme added the
# matter form below 1+z = 131 while the N_eff radiation term - which itself
# contains the neutrinos - was kept at all redshifts: below the join the
# neutrinos were counted twice (the spurious remnant reaches ~0.5% of E^2
# near z ~ 70 and shifts D_M(z_rec) by ~7e-4). The corrected scheme uses
# the photons-only radiation density below the join and the full N_eff
# density above it, with the two asymptotic forms matched continuously at
# their crossing, 1 + z = OMEGA_NU/(OMEGA_R_H2 - OMEGA_GAMMA_H2)
# (h-independent; T ~ m_nu/10). The residual transition-region
# approximation affects E^2 at the 1e-4 level.
Z_NR = OMEGA_NU / (OMEGA_R_H2 - OMEGA_GAMMA_H2)  # 1+z = 37.7

DATA_VEC = np.array([1.7504, 301.77, 0.022371])
COV = np.array([[1559.83, -1325.41, -36.45],
                [-1325.41, 714691.80, 269.77],
                [-36.45, 269.77, 2.10]]) * 1e-8
COV_INV = np.linalg.inv(COV)


def chi2(vec):
    d = vec - DATA_VEC
    return float(d @ COV_INV @ d)


def nu_density(z, wnu_h):
    """Massive-neutrino density (units of h^-2 relative, i.e. Omega_nu(z)).
    Below 1+z = Z_NR the neutrinos are non-relativistic matter; above Z_NR
    their relativistic density is already contained in N_eff = 3.046 of the
    radiation term, so no extra contribution is added (no double counting)."""
    zp1 = 1.0 + z
    if zp1 <= Z_NR:
        return wnu_h * zp1 ** 3
    return 0.0


def radiation_density(z, h):
    """Radiation density: photons only below the join, photons + N_eff above."""
    zp1 = 1.0 + z
    if zp1 <= Z_NR:
        return OMEGA_GAMMA_H2 / h ** 2 * zp1 ** 4
    return OMEGA_R_H2 / h ** 2 * zp1 ** 4


def sound_speed(z, wb):
    R_b = 31500.0 * wb * (T_CMB / 2.7) ** (-4) / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R_b))


# --------------------------------------------------------------------------
# Single source of the non-DE density for the CMB path.
# --------------------------------------------------------------------------

def omega_ng(z, Om, h):
    """Non-DE density Omega_ng(z; Omega_m, h) of SPEC section 1.

    Om is the bimetric matter density (Omega_m), the radiation and the
    massive-neutrino terms are the verbatim source functions
    radiation_density() and nu_density().  Used in the quartic, in E^2, in
    the closure and in both execution paths.
    """
    return (Om * (1 + z) ** 3
            + radiation_density(z, h)
            + nu_density(z, OMEGA_NU / h ** 2))


def solve_y_cmb(m, z, Om, h):
    """Smallest positive real root of the quartic, CMB (consistent) density.

    The coefficient c1 carries the FULL non-DE density, Omega_ng, i.e.
    c1 = t (B0 + 3 Omega_ng) - 3 B2.  The default Model.solve_y keeps the
    source behaviour and is not used on the CMB path.
    """
    if m.mFP2 == 0.0:
        return 1.0                      # GR limit
    Om_ng = omega_ng(z, Om, h)
    c4, c3, c2, c1, c0 = m.quartic_coeffs(Om_ng)
    f = lambda y: ((((c4 * y + c3) * y + c2) * y + c1) * y + c0)
    try:
        return brentq(f, Y_LO, 1.0, xtol=1e-14, rtol=1e-14)
    except ValueError as exc:
        # At sufficiently high z the positive root can be below Y_LO.
        # Extend only the numerical bracket; never replace the bimetric
        # background by GR and never accept a point without a sign change.
        f0, f1 = f(0.0), f(1.0)
        if not np.isfinite(f0) or not np.isfinite(f1) or f0 * f1 >= 0.0:
            raise ValueError("no finite-branch root bracket on [0, 1]") from exc
        return brentq(f, 0.0, 1.0, xtol=5e-324, rtol=1e-14)


def bim_ode_y(m, y):
    """Dark-energy branch evaluated at a given y (dy = 1 - y)."""
    dy = 1.0 - y
    return (m.Omega_L - m.sin2 * m.mFP2 * dy
            * (1 + m.alpha * dy + m.beta / 3.0 * dy ** 2))


def bim_ode(m, z, Om, h):
    """Dark-energy branch of the CMB path (consistent quartic)."""
    return bim_ode_y(m, solve_y_cmb(m, z, Om, h))


def E2_bim_of(m, Om, h, z):
    """E^2(z) on the CMB path: Omega_ng(z) + Omega_DE[y(z)]."""
    return omega_ng(z, Om, h) + bim_ode(m, z, Om, h)


def solve_Om0_nu(m, h):
    """z = 0 closure with the massive-neutrino term (CMB convention).

    The closure residual is 1 - Omega_ng(0) - Omega_DE[y(0)], with y(0) from
    the SAME consistent quartic solver that is used on the direct and spline
    execution paths.
    """
    def eq(Om):
        ODE = bim_ode(m, 0.0, Om, h)
        return 1.0 - omega_ng(0.0, Om, h) - ODE
    return brentq(eq, 0.01, 0.99, xtol=1e-13, rtol=1e-13)


def dm_rs(E2f, wb, E2f_rs=None, z_star=Z_STAR):
    """(D_M, r_s) in units of c/H0. D_M integrates the supplied E^2(z) out to
    recombination; r_s integrates from recombination to the early Universe,
    where the background is Vainshtein-screened (GR-like), so the screened
    form E2f_rs is used (defaults to E2f)."""
    DM, _ = quad(lambda z: 1.0 / np.sqrt(E2f(z)), 0.0, z_star,
                 epsabs=1e-11, epsrel=1e-10, limit=300)
    if E2f_rs is None:
        E2f_rs = E2f

    def rs_integrand(a):
        z = 1.0 / a - 1.0
        return sound_speed(z, wb) / (a * a * np.sqrt(E2f_rs(z)))
    r_s, _ = quad(rs_integrand, 1e-8, 1.0 / (1.0 + z_star),
                  epsabs=1e-11, epsrel=1e-10, limit=300)
    return DM, r_s


_ZG = np.linspace(0.0, 1200.0, 4001)
_AG = np.logspace(np.log10(1e-8), np.log10(1.0), 4000)
_Y_CACHE = {}


def _y_spline(m, Om, h):
    """Cubic spline of the quartic root y(z) on the fine grid.

    On the corrected CMB path the quartic depends on h (through Omega_ng), so
    the cache key contains h together with every parameter that enters
    Omega_ng: (theta, m_FP, alpha, beta, Omega_Lambda, Omega_m, h).  Two
    different h therefore can never share one spline.  Content-addressed:
    id() reuse after garbage collection must never alias two different
    models.  The spline is used ONLY to make quad evaluation fast.
    """
    key = (m.theta_deg, m.mFP_H0, m.alpha, m.beta, m.Omega_L,
           round(Om, 9), float(h))
    if key not in _Y_CACHE:
        yg = np.array([solve_y_cmb(m, z, Om, h) for z in _ZG])
        _Y_CACHE[key] = CubicSpline(_ZG, yg)
    return _Y_CACHE[key]


def vec_bimetric(theta_deg, h, wb, Omega_L=0.699, r_map="quartic"):
    m = Model(theta_deg=theta_deg, Omega_L=Omega_L)
    return vec_bimetric_m(m, h, wb, r_map)


def vec_bimetric_m(m, h, wb, r_map="quartic"):
    """Data vector with D_M and r_s from adaptive quadrature (quad) -- the
    same scheme as the LambdaCDM reference (dm_rs). The pre-8 uniform-grid
    trapezoid + cubic-interpolation D_M overestimated D_M(z_rec) by ~0.1%
    even for pure LambdaCDM and mixed schemes between the bimetric block
    and the reference; both are removed."""
    Om = solve_Om0_nu(m, h)
    yspl = _y_spline(m, Om, h)

    def E2_of(z):
        dy = 1.0 - float(yspl(z))
        ODE = (m.Omega_L - m.sin2 * m.mFP2 * dy
               * (1 + m.alpha * dy + m.beta / 3.0 * dy ** 2))
        return omega_ng(z, Om, h) + ODE

    zs = z_star_hs(wb, Om * h ** 2 + OMEGA_NU, h)
    DM, _ = quad(lambda z: 1.0 / np.sqrt(E2_of(z)), 0.0, zs,
                 epsabs=1e-11, epsrel=1e-10, limit=300)

    def E2scr_a(a):
        z = 1.0 / a - 1.0
        return omega_ng(z, Om, h) + m.Omega_L

    a_rec = 1.0 / (1.0 + zs)
    r_s, _ = quad(lambda a: sound_speed(1.0 / a - 1.0, wb)
                  / (a * a * np.sqrt(E2scr_a(a))), 1e-8, a_rec,
                  epsabs=1e-11, epsrel=1e-10, limit=300)

    if isinstance(r_map, float):
        Om_for_R = (1.0 - r_map * m.sin2) * Om
    elif r_map == "quartic":
        Om_for_R = Om
    else:
        Om_for_R = Om * m.cos2
    R = np.sqrt(Om_for_R + OMEGA_NU / h ** 2) * DM
    return np.array([R, np.pi * DM / r_s, wb]), Om


def fit_bimetric(theta_deg, Omega_L=0.699, r_map="quartic"):
    m = Model(theta_deg=theta_deg, Omega_L=Omega_L)
    def obj(p):
        h, wb = p
        if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
            return 1e10
        return chi2(vec_bimetric_m(m, h, wb, r_map)[0])
    best, x0 = 1e10, None
    for h in np.linspace(0.60, 0.85, 18):
        for wb in np.linspace(0.016, 0.024, 13):
            c = obj([h, wb])
            if c < best:
                best, x0 = c, [h, wb]
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 4000})
    h, wb = res.x
    vec, Om = vec_bimetric_m(m, h, wb, r_map)
    return chi2(vec), h, wb, Om, vec


def vec_lcdm(wm, h, wb):
    Om_t = wm / h ** 2
    OL = 1.0 - Om_t - radiation_density(0.0, h) - OMEGA_NU / h ** 2

    def E2f(z):
        return (Om_t * (1 + z) ** 3 + radiation_density(z, h) + OL
                + nu_density(z, OMEGA_NU / h ** 2))
    # omega_m for z_star_hs INCLUDES the massive neutrinos (as in
    # vec_bimetric_m and in Bansal & Huterer's compression definition);
    # the R prefactor is sqrt((Omega_cb + omega_nu/h^2)) with the same
    # closure convention as the bimetric path.
    zs = z_star_hs(wb, wm + OMEGA_NU, h)
    DM, r_s = dm_rs(E2f, wb, E2f_rs=E2f, z_star=zs)
    return np.array([np.sqrt(wm / h ** 2 + OMEGA_NU / h ** 2) * DM,
                     np.pi * DM / r_s, wb])


def fit_lcdm_reference():
    def obj(p):
        wm, h, wb = p
        if not (0.10 < wm < 0.25 and 0.55 < h < 0.85 and 0.015 < wb < 0.025):
            return 1e10
        return chi2(vec_lcdm(wm, h, wb))
    best, x0 = 1e10, None
    for h in np.linspace(0.62, 0.80, 10):
        for wb in np.linspace(0.020, 0.0235, 8):
            for wm in (0.14, 0.16):
                c = obj([wm, h, wb])
                if c < best:
                    best, x0 = c, [wm, h, wb]
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 6000})
    return res.fun, res.x


def split_test(theta_deg, z_screen=12.0):
    """Isolate the source of the theta penalty in the D_M integral: the
    bimetric residual relative to the screened background is zeroed above
    z_screen (zero_high) or below it (zero_low); the acoustic scale l_a is
    held at its full-mode value in every mode (the screened sound horizon).
    The legacy mixed-slot chi2 fields combine the modified R, full-mode ell_a
    and observed omega_b. They are not refit goodness-of-fit values, are not
    additive across modes and are not used for the manuscript conclusion;
    only the R shifts are interpreted."""
    c_full, h, wb, Om, vec = fit_bimetric(theta_deg)
    m = Model(theta_deg=theta_deg)

    def E2_full(z):
        return E2_bim_of(m, Om, h, z)

    def E2_screened(z):
        return omega_ng(z, Om, h) + m.Omega_L

    out = {"full": c_full}
    for name, keep_low in (("zero_high", True), ("zero_low", False)):
        def E2f(z, keep_low=keep_low):
            res = E2_full(z) - E2_screened(z)
            if (z <= z_screen) == keep_low:
                return E2_screened(z) + res
            return E2_screened(z)
        DM, _ = quad(lambda z: 1.0 / np.sqrt(E2f(z)), 0, z_star_hs(wb, Om*h**2 + OMEGA_NU, h), limit=200)
        R = np.sqrt(Om + OMEGA_NU / h ** 2) * DM
        out[name] = chi2(np.array([R, vec[1], DATA_VEC[2]]))
        out[name + "_R"] = R
    return out


def theta_scan(thetas=(6.1, 10.0, 12.2, 15.0, 18.3, 20.0, 25.0, 30.0),
               r_map="quartic"):
    chi2_ref, p_ref = fit_lcdm_reference()
    out = {}
    for th in thetas:
        c2, h, wb, Om, vec = fit_bimetric(th, r_map=r_map)
        out[th] = dict(chi2=c2, dchi2=c2 - chi2_ref, h=h, wb=wb, Om=Om,
                       R=vec[0], la=vec[1])
    return chi2_ref, p_ref, out


if __name__ == "__main__":
    c_ref, p_ref = fit_lcdm_reference()
    print("LambdaCDM reference: chi2 = %.3e  (wm, h, wb) = (%.6f, %.6f, %.6f)"
          % (c_ref, p_ref[0], p_ref[1], p_ref[2]))
    _, _, scan = theta_scan()
    for th, d in scan.items():
        print("  theta=%5.1f  chi2=%10.5f  dchi2=%10.5f  h=%.5f  Om=%.5f  R=%.5f  la=%.3f"
              % (th, d["chi2"], d["dchi2"], d["h"], d["Om"], d["R"], d["la"]))
    for th in (10.0, 30.0):
        print("  split %4.1f: %s" % (th, {k: round(v, 5) for k, v in split_test(th).items()}))
