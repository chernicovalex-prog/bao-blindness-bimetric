"""Bimetric background core (Hassan-Rosen, H&M 2021/2025 conventions).

Single canonical implementation. Units H0 = 1 throughout; m_FP is given in
units of H0. The scale-factor ratio y = a_f/a_g solves the quartic
constraint (H&M 2025, Eq. A13); the physical finite branch is the smallest
positive real root (H&M 2021, Sec. 5.2).
"""
import numpy as np
from scipy.optimize import brentq

C_KM_S_MPC = 299792.458          # c in (km/s)/Mpc
OMEGA_R_NU = 9.08e-5             # photons + neutrinos (N_eff = 3.046)
OMEGA_R_PHOTONS = 5.38e-5        # photons only (DESI DR2 fiducial)

Y_LO = 1e-12                     # bracket lower end (contains the full branch)


class Model:
    """Bimetric benchmark point."""

    def __init__(self, theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
                 Omega_L=0.699, Omega_r=OMEGA_R_NU):
        self.theta_deg = float(theta_deg)
        self.mFP_H0 = float(mFP_H0)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.Omega_L = float(Omega_L)
        self.Omega_r = float(Omega_r)
        th = np.radians(theta_deg)
        self.sin2 = np.sin(th) ** 2
        self.cos2 = np.cos(th) ** 2
        self.tan2 = np.tan(th) ** 2
        self.mFP2 = self.mFP_H0 ** 2

    def B_params(self):
        """B-parameters, H&M 2021 Eq. (3.4a-e)."""
        s2m2 = self.sin2 * self.mFP2
        B0 = 3 * self.Omega_L - s2m2 * (3 + 3 * self.alpha + self.beta)
        B1 = s2m2 * (1 + 2 * self.alpha + self.beta)
        B2 = -s2m2 * (self.alpha + self.beta)
        B3 = s2m2 * self.beta
        B4 = 3 * self.tan2 * self.Omega_L + s2m2 * (-1 + self.alpha - self.beta)
        return np.array([B0, B1, B2, B3, B4])

    def quartic_coeffs(self, Om):
        """Quartic coefficients (c4, c3, c2, c1, c0), H&M 2021 Eq. (5.5),
        with the full non-DE density Om = Omega_m (1+z)^3 + Omega_r (1+z)^4."""
        B0, B1, B2, B3, B4 = self.B_params()
        t = self.tan2
        c4 = t * B3
        c3 = 3 * t * B2 - B4
        c2 = 3 * (t * B1 - B3)
        c1 = t * (B0 + 3.0 * Om) - 3 * B2
        c0 = -B1
        return np.array([c4, c3, c2, c1, c0], dtype=float)

    def solve_y(self, z, Om0):
        """Smallest positive real root of the quartic at redshift z."""
        if self.mFP2 == 0.0:
            return 1.0                      # GR limit
        Om = Om0 * (1 + z) ** 3 + self.Omega_r * (1 + z) ** 4
        c4, c3, c2, c1, c0 = self.quartic_coeffs(Om)
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

    def E2_bim(self, z, Om0):
        """E^2(z) = H^2/H0^2, full bimetric (H&M 2025 Eqs. A11-A12)."""
        y = self.solve_y(z, Om0)
        dy = 1.0 - y
        ODE = (self.Omega_L - self.sin2 * self.mFP2 * dy
               * (1 + self.alpha * dy + self.beta / 3.0 * dy ** 2))
        return Om0 * (1 + z) ** 3 + self.Omega_r * (1 + z) ** 4 + ODE

    def solve_Om0(self):
        """Self-consistent matter density from the z = 0 closure E^2(0) = 1."""
        def eq(Om):
            y = self.solve_y(0.0, Om)
            dy = 1.0 - y
            ODE = (self.Omega_L - self.sin2 * self.mFP2 * dy
                   * (1 + self.alpha * dy + self.beta / 3.0 * dy ** 2))
            return 1.0 - self.Omega_r - ODE - Om
        return brentq(eq, 0.01, 0.99, xtol=1e-14, rtol=1e-14)

    def w_DE(self, z, Om0, dz=1e-4):
        """Effective dark-energy equation of state (H&M 2021 Eq. 5.9)."""
        def ODE(y):
            d = 1.0 - y
            return (self.Omega_L - self.sin2 * self.mFP2 * d
                    * (1 + self.alpha * d + self.beta / 3.0 * d ** 2))
        y0 = self.solve_y(z, Om0)
        yp = self.solve_y(z + dz, Om0)
        if z < dz:
            dODE = (ODE(yp) - ODE(y0)) / dz
        else:
            ym = self.solve_y(z - dz, Om0)
            dODE = (ODE(yp) - ODE(ym)) / (2 * dz)
        return -1.0 + (1 + z) / (3.0 * ODE(y0)) * dODE


def E2_lcdm(z, Om0, Omega_r=OMEGA_R_NU):
    return Om0 * (1 + z) ** 3 + Omega_r * (1 + z) ** 4 + (1 - Om0 - Omega_r)
