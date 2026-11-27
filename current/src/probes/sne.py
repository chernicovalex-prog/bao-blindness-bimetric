"""SNe Ia probe: Pantheon+ (1701 light curves; 1590 after the zHD > 0.01 cut).

Two statistically consistent covariance treatments: the diagonal approximation
(m_b_corr_err_DIAG; upper-bound diagnostic) and the full STAT+SYS covariance
(primary). The absolute magnitude M is profiled analytically in both. The
bimetric Omega_m is fixed by the quartic closure; only M is profiled for the
bimetric model, while LambdaCDM is also fitted in Omega_m.
"""
import numpy as np
import pathlib, sys
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from core import Model, OMEGA_R_NU

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
C_H0 = 2997.92458          # c/H0 in Mpc (h = 1)
Z_FINE = np.linspace(0.0, 5.0, 10000)

DF = None


def load_pantheon():
    """Pantheon+SH0ES public data release (see README for provenance)."""
    global DF
    cols = ["CID", "IDSURVEY",
            "zHD", "zHDERR", "zCMB", "zCMBERR", "zHEL", "zHELERR",
            "m_b_corr", "m_b_corr_err_DIAG", "MU_SH0ES", "MU_SH0ES_ERR_DIAG",
            "CEPH_DIST", "IS_CALIBRATOR", "USED_IN_SH0ES",
            "x0", "x1", "c", "HOST_MASS"]
    DF = {}
    raw = [l.split() for l in open(DATA / "Pantheon+SH0ES.dat") if l.strip()
           and not l.startswith("#")]
    # numeric block starts at token 2; keep the columns this probe needs
    need = {"zHD": 2, "zCMB": 4, "zHEL": 6, "m_b_corr": 8, "m_b_corr_err_DIAG": 9}
    DF["CID"] = [t[0] for t in raw]
    for name, i in need.items():
        DF[name] = np.array([float(t[i]) for t in raw])
    cov = np.loadtxt(DATA / "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1)
    n = int(round(np.sqrt(cov.size)))
    DF["COV_FULL"] = cov.reshape(n, n)
    return DF


def bimetric_mu(zcut):
    """Bimetric distance modulus at the baseline point (Omega_m from the
    quartic closure, Omega_r = 9.08e-5; zHD in the integral (Pantheon+ standard), zHEL in 1+z)."""
    m = Model(Omega_r=OMEGA_R_NU)
    Om = m.solve_Om0()
    E2g = np.array([m.E2_bim(z, Om) for z in Z_FINE])
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2g), Z_FINE, initial=0)
    DMi = interp1d(Z_FINE, ic, kind="cubic", fill_value="extrapolate")
    zc, zh = DF["zHD"][zcut], DF["zHEL"][zcut]
    DL = (1.0 + zh) * C_H0 * DMi(zc)
    return 5.0 * np.log10(DL) + 25.0, Om


def lcdm_mu(Om, zcut):
    E2g = Om * (1 + Z_FINE) ** 3 + OMEGA_R_NU * (1 + Z_FINE) ** 4         + (1 - Om - OMEGA_R_NU)
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2g), Z_FINE, initial=0)
    zc, zh = DF["zHD"][zcut], DF["zHEL"][zcut]
    DL = (1.0 + zh) * C_H0 * np.interp(zc, Z_FINE, ic)
    return 5.0 * np.log10(DL) + 25.0


def chi2_diag(mu, zcut, sigma_int=0.0):
    r = DF["m_b_corr"][zcut] - mu
    s2 = DF["m_b_corr_err_DIAG"][zcut] ** 2 + sigma_int ** 2
    w = 1.0 / s2
    return float(np.sum(w * r ** 2) - np.sum(w * r) ** 2 / np.sum(w))


def chi2_full(mu, zcut):
    cov = DF["COV_FULL"][np.ix_(zcut, zcut)]
    Cinv = np.linalg.inv(cov)
    r = DF["m_b_corr"][zcut] - mu
    Cr = Cinv @ r
    return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv))


def fit_lcdm(chi2, zcut):
    res = minimize_scalar(lambda Om: chi2(lcdm_mu(Om, zcut), zcut),
                          bounds=(0.2, 0.5), method="bounded",
                          options={"xatol": 1e-12})
    return float(res.x), float(res.fun)


def run(zcut_value=0.01, sigma_int=0.0):
    zcut = DF["zHD"] > zcut_value
    mu_b, Om_b = bimetric_mu(zcut)
    c2d = lambda mu, zc: chi2_diag(mu, zc, sigma_int)
    Om_d, chi2_d = fit_lcdm(c2d, zcut)
    chi2_bim_d = c2d(mu_b, zcut)
    chi2_fixed_d = c2d(lcdm_mu(Om_b, zcut), zcut)

    cov = DF["COV_FULL"][np.ix_(zcut, zcut)]
    Cinv = np.linalg.inv(cov)
    def c2f(mu, zc):
        r = DF["m_b_corr"][zc] - mu
        Cr = Cinv @ r
        return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv))
    Om_f, chi2_f = fit_lcdm(c2f, zcut)
    chi2_bim_f = c2f(mu_b, zcut)
    chi2_fixed_f = c2f(lcdm_mu(Om_b, zcut), zcut)

    return dict(N=int(zcut.sum()), Om_bim=Om_b,
                diag=dict(Om_lcdm=Om_d, chi2_min=chi2_d, chi2_bim=chi2_bim_d,
                          dchi2=chi2_bim_d - chi2_d,
                          pure=chi2_bim_d - chi2_fixed_d,
                          chi2_dof=chi2_d / (int(zcut.sum()) - 2)),
                full=dict(Om_lcdm=Om_f, chi2_min=chi2_f, chi2_bim=chi2_bim_f,
                          dchi2=chi2_bim_f - chi2_f,
                          pure=chi2_bim_f - chi2_fixed_f,
                          chi2_dof=chi2_f / (int(zcut.sum()) - 2)))


def run_zcut_on(zcol="zHD", zcut_value=0.01):
    """Robustness of the event selection: the z > 0.01 low-z cut applied to
    the zCMB column instead of zHD while the distance integral stays on the
    zHD primary convention (only the event selection changes)."""
    return run_on_mask(DF[zcol] > zcut_value)


def run_on_mask(mask):
    """The standard z>0.01 comparison restricted to an explicit mask."""
    mu_b, Om_b = bimetric_mu(mask)
    c2d = lambda mu, zc: chi2_diag(mu, zc, 0.0)
    Om_d, chi2_d = fit_lcdm(c2d, mask)
    chi2_bim_d = c2d(mu_b, mask)
    chi2_fixed_d = c2d(lcdm_mu(Om_b, mask), mask)
    cov = DF["COV_FULL"][np.ix_(mask, mask)]
    Cinv = np.linalg.inv(cov)
    def c2f(mu, zc):
        r = DF["m_b_corr"][zc] - mu
        Cr = Cinv @ r
        return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv))
    Om_f, chi2_f = fit_lcdm(c2f, mask)
    chi2_bim_f = c2f(mu_b, mask)
    chi2_fixed_f = c2f(lcdm_mu(Om_b, mask), mask)
    return dict(N=int(mask.sum()), Om_bim=Om_b,
                diag=dict(Om_lcdm=Om_d, chi2_min=chi2_d, chi2_bim=chi2_bim_d,
                          dchi2=chi2_bim_d - chi2_d,
                          pure=chi2_bim_d - chi2_fixed_d,
                          chi2_dof=chi2_d / (int(mask.sum()) - 2)),
                full=dict(Om_lcdm=Om_f, chi2_min=chi2_f, chi2_bim=chi2_bim_f,
                          dchi2=chi2_bim_f - chi2_f,
                          pure=chi2_bim_f - chi2_fixed_f,
                          chi2_dof=chi2_f / (int(mask.sum()) - 2)))


def bimetric_mu_at(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0):
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_r=OMEGA_R_NU)
    Om = m.solve_Om0()
    E2g = np.array([m.E2_bim(z, Om) for z in Z_FINE])
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2g), Z_FINE, initial=0)
    DMi = interp1d(Z_FINE, ic, kind="cubic", fill_value="extrapolate")
    zc, zh = DF["zHD"], DF["zHEL"]
    DL = (1.0 + zh) * C_H0 * DMi(zc)
    return 5.0 * np.log10(DL) + 25.0


def scan_grid(full=True):
    """(alpha,beta) and m_FP scans at the baseline point; full or diagonal."""
    zcut = DF["zHD"] > 0.01
    if full:
        cov = DF["COV_FULL"][np.ix_(zcut, zcut)]
        Cinv = np.linalg.inv(cov)
        def f(mu):
            r = DF["m_b_corr"][zcut] - mu
            Cr = Cinv @ r
            return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv))
    else:
        def f(mu):
            return chi2_diag(mu, zcut, 0.0)
    out = {}
    for a in (-1.0, -0.5, 0.0, 0.5, 1.0):
        for b in (1, 2, 3, 4):
            try:
                mu = bimetric_mu_at(alpha=a, beta=b)[zcut]
            except ValueError:
                out[("ab", a, b)] = None
                continue
            res = minimize_scalar(lambda Om: f(lcdm_mu(Om, zcut)),
                                  bounds=(0.2, 0.5), method="bounded",
                                  options={"xatol": 1e-12})
            out[("ab", a, b)] = (f(mu) - res.fun, res.x)
    for mv in (15.0, 25.0, 45.0, 100.0, 1000.0, 10000.0, 100000.0):
        mu = bimetric_mu_at(mFP_H0=mv)[zcut]
        res = minimize_scalar(lambda Om: f(lcdm_mu(Om, zcut)),
                              bounds=(0.2, 0.5), method="bounded",
                              options={"xatol": 1e-12})
        out[("mFP", mv)] = (f(mu) - res.fun, res.x)
    return out


def sigma_int_scan(values=(0.0, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20)):
    out = {}
    for s in values:
        out[s] = run(0.01, s)["diag"]
    return out


if __name__ == "__main__":
    load_pantheon()
    print("N total:", len(DF["CID"]))
    r = run(0.01)
    print("z>0.01 (N=%d)  Om_bim=%.5f" % (r["N"], r["Om_bim"]))
    for scheme in ("diag", "full"):
        d = r[scheme]
        print("  %-4s Om_lcdm=%.4f chi2_min=%.3f chi2_bim=%.3f dchi2=%.3f "
              "pure=%.3f chi2/dof=%.3f"
              % (scheme, d["Om_lcdm"], d["chi2_min"], d["chi2_bim"],
                 d["dchi2"], d["pure"], d["chi2_dof"]))
    r2 = run(0.023)
    d = r2["full"]
    print("z>0.023 (N=%d) full: Om=%.4f dchi2=%.3f pure=%.3f"
          % (r2["N"], d["Om_lcdm"], d["dchi2"], d["pure"]))
    print("sigma_int scan (diag):")
    for s, d in sigma_int_scan().items():
        print("  si=%.2f dchi2=%.3f pure=%.3f Om=%.4f chi2/dof=%.3f"
              % (s, d["dchi2"], d["pure"], d["Om_lcdm"], d["chi2_dof"]))
