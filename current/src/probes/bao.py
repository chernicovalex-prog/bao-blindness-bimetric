"""BAO probe: DESI DR2 (13 points, full official covariance, Cobaya release).

The scale f = c/(H0 r_d) is profiled by chi^2 minimisation (Brent,
bounds [20, 80]). All BAO numbers use the photon-only radiation density
Omega_r = 5.38e-5 (DESI DR2 fiducial convention); the effect of the
neutrino contribution at z < 2.33 is < 0.01 in Delta chi^2.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
import pathlib, json, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from core import Model, E2_lcdm, OMEGA_R_PHOTONS

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
ZG = np.linspace(0.0, 2.5, 2001)

OBS = None          # list of (z, type, value)
CINV = None


def load_desi():
    """Official DESI DR2 mean vector and full 13x13 covariance."""
    global OBS, CINV
    obs = []
    for line in open(DATA / "desi_gaussian_bao_ALL_GCcomb_mean.txt"):
        if line.startswith("#") or not line.strip():
            continue
        z, val, qty = line.split()
        obs.append((float(z), qty.split("_")[0], float(val)))
    assert len(obs) == 13
    cov = np.loadtxt(DATA / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    assert cov.shape == (13, 13)
    OBS = obs
    CINV = np.linalg.inv(cov)
    return obs, np.linalg.inv(cov)


def model_from_E2(E2g, obs_list=None):
    """Model (DH, DM, DV) vector as a function of the scale f."""
    obs_list = OBS if obs_list is None else obs_list
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2g), ZG, initial=0)
    Ei = interp1d(ZG, E2g, kind="cubic")
    DMi = interp1d(ZG, ic, kind="cubic")

    def model(f):
        out = []
        for z, typ, _ in obs_list:
            if typ == "DH":
                out.append(f / np.sqrt(Ei(z)))
            elif typ == "DM":
                out.append(f * DMi(z))
            else:
                out.append(f * (z * DMi(z) ** 2 / np.sqrt(Ei(z))) ** (1.0 / 3.0))
        return np.array(out)
    return model


def chi2(model, f):
    d = np.array([o[2] for o in OBS]) - model(f)
    return float(d @ CINV @ d)


def fit_f(model):
    res = minimize_scalar(lambda f: chi2(model, f), bounds=(20, 80), method="bounded")
    return float(res.x), float(res.fun)


def E2_bim_grid(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0, Om0=None):
    m = Model(theta_deg=theta_deg, mFP_H0=mFP_H0, alpha=alpha, beta=beta,
              Omega_r=OMEGA_R_PHOTONS)
    if Om0 is None:
        Om0 = m.solve_Om0()
    E2g = np.array([m.E2_bim(z, Om0) for z in ZG])
    return E2g, Om0


def fit_lcdm():
    """LCDM best fit: Omega_m and f profiled."""
    def obj(o):
        return fit_f(model_from_E2(E2_lcdm(ZG, o, OMEGA_R_PHOTONS)))[1]
    res = minimize_scalar(obj, bounds=(0.25, 0.40), method="bounded",
                          options={"xatol": 1e-12})
    Om = float(res.x)
    f, c2 = fit_f(model_from_E2(E2_lcdm(ZG, Om, OMEGA_R_PHOTONS)))
    return Om, f, c2


def baseline():
    """Baseline Delta chi^2 with the Omega_m decomposition."""
    E2g, Om_bim = E2_bim_grid()
    f_b, c_b = fit_f(model_from_E2(E2g))
    Om_l, f_l, c_l = fit_lcdm()
    c_l_at_bim = fit_f(model_from_E2(E2_lcdm(ZG, Om_bim, OMEGA_R_PHOTONS)))[1]
    E2g_lbest, _ = E2_bim_grid(Om0=Om_l)
    c_b_at_l = fit_f(model_from_E2(E2g_lbest))[1]
    return dict(Om_bim=Om_bim, f_bim=f_b, chi2_bim=c_b,
                Om_lcdm=Om_l, f_lcdm=f_l, chi2_lcdm=c_l,
                Om_shift_penalty=c_l_at_bim - c_l,
                bim_correction=c_b - c_l_at_bim,
                dchi2=c_b - c_l,
                pure=c_b_at_l - c_l)


def scan_alpha_beta():
    out = {}
    c_l = fit_lcdm()[2]
    for a in (-1.0, -0.5, 0.0, 0.5, 1.0):
        for b in (1, 2, 3, 4):
            try:
                E2g, Om = E2_bim_grid(alpha=a, beta=b)
            except ValueError:
                out[(a, b)] = ("NO_FINITE_BRANCH", None)   # B1 = 0 locus
                continue
            c_b = fit_f(model_from_E2(E2g))[1]
            out[(a, b)] = (c_b - c_l, Om)
    return out


def scan_mFP(mvals=(15, 25, 45, 100, 1000, 10000, 100000)):
    out = {}
    c_l = fit_lcdm()[2]
    for mv in mvals:
        E2g, Om = E2_bim_grid(mFP_H0=mv)
        out[mv] = (fit_f(model_from_E2(E2g))[1] - c_l, Om)
    return out


def scan_theta(thetas=(6.1, 10.0, 15.0, 20.0, 25.0, 31.5)):
    """Self-consistent scan (Omega_m refit per theta) and the fixed-Omega_m
    variant (bimetric Omega_m frozen at its baseline value)."""
    c_l_ref = fit_lcdm()[2]
    _, Om_bim0 = E2_bim_grid()
    out = {}
    for th in thetas:
        E2g, Om = E2_bim_grid(theta_deg=th)
        d_self = fit_f(model_from_E2(E2g))[1] - c_l_ref
        E2g_f, _ = E2_bim_grid(theta_deg=th, Om0=Om_bim0)
        d_fixed = fit_f(model_from_E2(E2g_f))[1] - c_l_ref
        out[th] = (d_self, d_fixed, Om)
    return out


def synthetic_null_refit(n=2000, seed=20261014):
    """LambdaCDM-null realisations with BOTH nuisance profiles refit per
    realisation: f for both models, and Omega_m for LambdaCDM (matching the
    observed-data comparison, where Omega_m is refit for LambdaCDM)."""
    E2g, Om_b = E2_bim_grid()
    model_b = model_from_E2(E2g)
    Om_l, f_l, c_l = fit_lcdm()
    model_l = model_from_E2(E2_lcdm(ZG, Om_l, OMEGA_R_PHOTONS))
    mu_l = model_l(f_l)
    cov = np.linalg.inv(CINV)
    rng = np.random.default_rng(seed)
    obs_d = baseline()["dchi2"]
    count = 0
    original_obs = list(OBS)
    for i in range(n):
        d = rng.multivariate_normal(mu_l, cov)
        OBS[:] = [(z, t, v) for (z, t, _), v in zip(original_obs, d)]
        dl = fit_lcdm()[2]
        db = fit_f(model_b)[1]
        if db - dl >= obs_d:
            count += 1
    OBS[:] = original_obs
    return count / n


def synthetic_null(n=10000, seed=20261014):
    """Fraction of LambdaCDM-null realisations with Delta chi^2 >= observed."""
    Om_l, f_l, c_l = fit_lcdm()
    E2g, Om_b = E2_bim_grid()
    model_l = model_from_E2(E2_lcdm(ZG, Om_l, OMEGA_R_PHOTONS))
    model_b = model_from_E2(E2g)
    dvec = np.array([o[2] for o in OBS])
    cov = np.linalg.inv(CINV)
    rng = np.random.default_rng(seed)
    obs_d = baseline()["dchi2"]
    count = 0
    original_obs = list(OBS)
    for i in range(n):
        d = rng.multivariate_normal(model_l(f_l), cov)
        OBS[:] = [(z, t, v) for (z, t, _), v in zip(original_obs, d)]
        dl = fit_f(model_l)[1]
        db = fit_f(model_b)[1]
        if db - dl >= obs_d:
            count += 1
    return count / n


if __name__ == "__main__":
    load_desi()
    bl = baseline()
    print("BAO baseline (official full covariance):")
    for k, v in bl.items():
        print("  %-16s = %.6f" % (k, v))
    ab = scan_alpha_beta()
    vals = [d for d, _ in ab.values() if isinstance(d, float)]
    print("(alpha,beta) scan: Delta chi^2 in [%.4f, %.4f] over %d points;"
          " degenerate locus B1=0 at (-1,1): %s"
          % (min(vals), max(vals), len(vals), ab.get((-1.0, 1))))
    ms = scan_mFP()
    print("m_FP scan:", {k: round(v[0], 4) for k, v in ms.items()})
    th = scan_theta()
    for t, (ds_, df_, om_) in th.items():
        print("  theta=%5.1f  self-consistent=%.4f  fixed-Om=%.4f  Om=%.5f" % (t, ds_, df_, om_))
    print("synthetic null refit p = %.4f" % synthetic_null_refit(n=2000, seed=20261014))
