"""Linear/WKB scalar-sector diagnostic (Akrami et al. 2015, Eq. 38).

Convention line: perturbations Phi ~ exp(i w N), N = ln a; instability when
X < 0; growth integral in the canonical 1/E convention,
    ln Gamma(k) = (k/H0) * I1,  I1 = int_{z*}^{1090} sqrt(-X)/E dz.
The dimensionally inconsistent 1/E^2 variant is not used.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

from core import Model

Z_MAX = 1090.0
I1_CANON = 0.945210109914
ZSTAR_CANON = 7.880536
X_CANON = {0.0: 0.9999955766, 5.0: 0.8399925, 10.0: -0.465722, 1090.0: -1.163845}


def akrami_mapping(m):
    """Mapping to Akrami conventions: alpha_A^2 = tan^2(theta),
    m^2 = B1 + 2 B2 + B3, beta_n^(A) = B_n / m^2 (gauge H0 = 1)."""
    alpha_A2 = m.tan2
    B = m.B_params()
    m2 = B[1] + 2 * B[2] + B[3]
    beta_A = B / m2
    return alpha_A2, m2, beta_A


def y_of_z(m, z, Om0):
    return m.solve_y(z, Om0)


def dy_dN(m, z, Om0):
    """dy/dN with N = ln a (exact implicit differentiation of the quartic;
    the y-factor in dP/dz is mandatory)."""
    y = m.solve_y(z, Om0)
    Om = Om0 * (1 + z) ** 3 + m.Omega_r * (1 + z) ** 4
    dOm_dz = 3 * Om0 * (1 + z) ** 2 + 4 * m.Omega_r * (1 + z) ** 3
    c4, c3, c2, c1, _ = m.quartic_coeffs(Om)
    dP_dy = c1 + 2 * c2 * y + 3 * c3 * y ** 2 + 4 * c4 * y ** 3
    dP_dz = 3 * m.tan2 * y * dOm_dz
    return -(1 + z) * (-dP_dz / dP_dy)


def X_of_z(m, z, Om0):
    """X(z) = (aH/k)^2 w^2, Akrami et al. 2015 Eq. (38)."""
    B = m.B_params()
    y = m.solve_y(z, Om0)
    yp = dy_dN(m, z, Om0)
    T1 = ((B[1] + 4 * B[2] * y + 3 * B[3] * y ** 2) * yp
          / (3 * y * (B[1] + 2 * B[2] * y + B[3] * y ** 2)))
    m2_rho_1pw = 3 * Om0 * (1 + z) ** 3 + 4 * m.Omega_r * (1 + z) ** 4
    T2 = ((1 + m.tan2 * y ** 2) * (B[1] - B[3] * y ** 2) * yp ** 2
          / (3 * m.tan2 * y ** 3 * m2_rho_1pw))
    return 1.0 + T1 - T2


def find_zstar(m, Om0, bracket=(5.0, 10.0)):
    """First positive-to-negative crossing of X(z) (brentq)."""
    return brentq(lambda z: X_of_z(m, z, Om0), bracket[0], bracket[1],
                  xtol=1e-12, rtol=1e-14)


def growth_integral(m, Om0, zstar):
    I1, err = quad(lambda z: np.sqrt(-X_of_z(m, z, Om0)) / np.sqrt(m.E2_bim(z, Om0)),
                   zstar, Z_MAX, epsabs=1e-12, epsrel=1e-12, limit=200)
    return I1, err


def log10_gamma(k, h, I1):
    H0 = 100.0 * h / 299792.458        # Mpc^-1
    return (k / H0) * I1 / np.log(10.0)


def run_scalar(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
               Omega_L=0.699, Omega_m0=None, Omega_r=9.08e-5):
    """Full scalar-sector evaluation at a model point. By default Omega_m0 is
    the self-consistent z = 0 closure value from the same code (primary);
    pass Omega_m0=0.304 for the legacy scalar-canon convention, which shifts
    z* to 7.880536 and I1 to 0.945210110 (documented difference, see paper)."""
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    if Omega_m0 is None:
        Omega_m0 = m.solve_Om0()
    alpha_A2, m2, beta_A = akrami_mapping(m)
    Xs = {z: X_of_z(m, z, Omega_m0) for z in X_CANON}
    zstar = find_zstar(m, Omega_m0)
    Hzstar = np.sqrt(m.E2_bim(zstar, Omega_m0))
    I1, err = growth_integral(m, Omega_m0, zstar)
    gamma = {(k, h): log10_gamma(k, h, I1)
             for k in (0.001, 0.002, 0.003, 0.01, 0.05, 0.10)
             for h in (0.700, 0.678)}
    return dict(alpha_A2=alpha_A2, m2=m2, beta_A=beta_A, X=Xs, zstar=zstar,
                Hzstar=Hzstar, I1=I1, I1_err=err, log10_gamma=gamma)


def E2_horizon(m, Om0, z):
    """Full bimetric E^2 for the horizon-entry condition k = aH.

    A quartic/bracket failure is a numerical failure and is propagated to
    the caller; it must never be silently converted into a GR background.
    """
    return m.E2_bim(z, Om0)


def k_over_aH(m, Om0, z, k):
    """Sub-horizon ratio k/(aH) for a mode k [Mpc^-1] at redshift z, with
    a = 1/(1+z): aH = H0 E(z)/(1+z), H0 = 100h/299792.458 Mpc^-1, h = 0.7
    for the canonical table. k/(aH) > 1 is the sub-horizon (WKB) side."""
    return (k * (1.0 + z)
            / ((100.0 * 0.7 / 299792.458) * np.sqrt(E2_horizon(m, Om0, z))))


def z_horizon_entry(m, Om0, k):
    """Redshift at which the mode k enters the horizon, k = aH."""
    f = lambda z: k_over_aH(m, Om0, z, k) - 1.0
    zg = np.geomspace(0.1, 1e7, 4000)
    fg = np.array([f(z) for z in zg])
    for i in range(len(zg) - 1):
        if fg[i] * fg[i + 1] < 0.0:
            return brentq(f, zg[i], zg[i + 1], xtol=1e-8)
    return None


def k_map(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0, Omega_L=0.699,
          Omega_r=9.08e-5, I1=0.9448705501542771, z_rec=1090.0,
          k_grid=(0.001, 0.05, 0.1, 0.5, 1.0, 5.0)):
    """The k-dilemma map (Sec. 6.2): horizon entry redshift z_h(k) of a mode
    against its formal growth exponent log10 Gamma(k) = (k/H0) I1/ln(10).
    z* is k-independent; Gamma is linear in k. No scale is simultaneously
    WKB-valid across the unstable window (z_h >>~ z*) and benign (Gamma << 1):
    validity to z_rec needs k >~ aH(z_rec), Gamma = 1 needs k ~ 6e-4 Mpc^-1."""
    m = Model(theta_deg=theta_deg, mFP_H0=float(mFP_H0), alpha=alpha, beta=beta,
              Omega_L=Omega_L, Omega_r=Omega_r)
    Om = m.solve_Om0()
    H0 = 100.0 * 0.7 / 299792.458
    rows = []
    for k in k_grid:
        rows.append(dict(k=float(k),
                         z_h=z_horizon_entry(m, Om, float(k)),
                         k_over_aH_rec=k_over_aH(m, Om, z_rec, float(k)),
                         log10Gamma=(float(k) / H0) * I1 / np.log(10.0)))
    k_rec = ((100.0 * 0.7 / 299792.458)
             * np.sqrt(E2_horizon(m, Om, z_rec)) / (1.0 + z_rec))
    k_g1 = H0 * np.log(10.0) / I1
    return dict(rows=rows, k_at_zrec=k_rec, k_gamma1=k_g1,
                z_h_gamma1=z_horizon_entry(m, Om, k_g1))


if __name__ == "__main__":
    r = run_scalar()
    print("scalar sector, baseline:")
    for z, x in r["X"].items():
        print("  X(%6.1f) = %+.8f" % (z, x))
    print("  z* = %.8f   H(z*)/H0 = %.4f" % (r["zstar"], r["Hzstar"]))
    print("  I1 = %.12f (+/- %.1e)" % (r["I1"], r["I1_err"]))
    for (k, h), lg in sorted(r["log10_gamma"].items()):
        print("  log10 Gamma(k=%.3f, h=%.3f) = %.2f" % (k, h, lg))
