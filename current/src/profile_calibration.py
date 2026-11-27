"""Monte Carlo calibration of the descriptive anchor profile.

Synthetic data are drawn from the bimetric benchmark itself (theta=6.1 deg,
Omega_Lambda = 0.699) with the real covariances; the full per-probe profiling
chain is rerun per realisation; the coverage of the dchi=1 / 3.84 bands around
the profile minimum at the true anchor is measured. Model predictions are
data-independent and precomputed once; per realisation only quadratic forms
remain (BAO: closed-form f per grid point and per LCDM Om; SNe: closed-form M
and a 1D LCDM minimisation; CMB: warm-started Nelder-Mead on (h, wb)).

Covariance convention: in BOTH the generator and the
likelihood the roles are canonical, exactly as in probes/*.py - synthetic
noise ~ N(0, C) with the published covariance C, chi^2 = r^T C^{-1} r. (The
first version held C under a Ci_* name and inverted both roles; that pair was
self-consistent - the profile LR at the truth is chi^2_1 for any matched
(noise, precision) pair - but the generator produced unrealistic data.)

GRID_MODE="offset" (argv[4]) shifts the anchor grid by +0.002 so the true
anchor is NOT a grid point; it measures the robustness of the coverage to the
grid resolution and to the interpolation of the minimum.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, ".")
sys.path.insert(0, "probes")
import probes.bao as bao
import probes.sne as sne
import probes.cmb as cmb
from core import Model
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d, CubicSpline

THETA, OL_TRUE = 6.1, 0.699
N_MC = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 20261118
TAG = sys.argv[3] if len(sys.argv) > 3 else "main"
# GRID_MODE="fine"   -> the quoted grid (spacing 0.004) with the true anchor
#                       as a grid point (canonical calibration, as in the paper);
# GRID_MODE="offset" -> the same grid shifted by +0.002 so the true anchor is
#                       NOT a grid point (offset-grid robustness run).
GRID_MODE = sys.argv[4] if len(sys.argv) > 4 else "fine"
_g = sorted(set(np.round(np.arange(0.662, 0.740, 0.004), 4)) | {round(OL_TRUE, 4)})
if GRID_MODE == "offset":
    _g = sorted(set(np.round(np.asarray(_g) + 0.002, 4)))
    assert round(OL_TRUE, 4) not in _g
OL_GRID = np.array(_g)

bao.load_desi()
sne.load_pantheon()
zg = bao.ZG
OR_BAO = 5.38e-5

dat_b = [l.split() for l in open("../data/desi_gaussian_bao_ALL_GCcomb_mean.txt") if not l.startswith("#")]
zb = np.array([float(x[0]) for x in dat_b])
kindb = np.array([x[2] for x in dat_b])
Cbi = np.linalg.inv(np.loadtxt("../data/desi_gaussian_bao_ALL_GCcomb_cov.txt"))
Cb = np.linalg.inv(Cbi)
Lb = np.linalg.cholesky(Cb)
kidx = np.array([{"DM_over_rs": 0, "DH_over_rs": 1, "DV_over_rs": 2}[k] for k in kindb])

zcut = sne.DF["zHD"] > 0.01
C_sne = sne.DF["COV_FULL"][np.ix_(zcut, zcut)]  # covariance (mag^2), as in sne.py
Cinv_sne = np.linalg.inv(C_sne)                # precision, as in sne.py
Lsne = np.linalg.cholesky(C_sne)               # synthetic noise ~ N(0, C_sne)

def bao_pred(OL):
    """BAO prediction using exactly the anchor_scan.bao_chi2 evaluation path."""
    m = Model(theta_deg=THETA, Omega_L=float(OL), Omega_r=OR_BAO)
    Om = m.solve_Om0()
    E2 = np.array([m.E2_bim(z, Om) for z in zg])
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2), zg, initial=0.0)
    Ei = interp1d(zg, E2, kind="cubic")
    DMi = interp1d(zg, ic, kind="cubic")
    DH_i = 1.0 / np.sqrt(Ei(zb))
    DM_i = DMi(zb)
    DV_i = (zb * DM_i ** 2 * DH_i) ** (1.0 / 3.0)
    return np.choose(kidx, [DM_i, DH_i, DV_i])

def sne_mu(OL):
    """SNe prediction using exactly the anchor_scan.sne_chi2 evaluation path."""
    m = Model(theta_deg=THETA, Omega_L=float(OL), Omega_r=sne.OMEGA_R_NU)
    Om = m.solve_Om0()
    E2 = np.array([m.E2_bim(z, Om) for z in sne.Z_FINE])
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2), sne.Z_FINE, initial=0.0)
    DMi = interp1d(sne.Z_FINE, ic, kind="cubic", fill_value="extrapolate")
    DL = (1.0 + sne.DF["zHEL"]) * sne.C_H0 * DMi(sne.DF["zHD"])
    return 5.0 * np.log10(DL) + 25.0

print("precomputing model grid...", flush=True)
cmb_fit = [cmb.fit_bimetric(THETA, float(OL), r_map="quartic") for OL in OL_GRID]
cmb_h = np.array([c[1] for c in cmb_fit]); cmb_wb = np.array([c[2] for c in cmb_fit])
def cmb_vec(j, h=None, wb=None):
    OL = float(OL_GRID[j])
    m = Model(theta_deg=THETA, Omega_L=OL, Omega_r=9.08e-5)
    return cmb.vec_bimetric_m(m, cmb_h[j] if h is None else h,
                              cmb_wb[j] if wb is None else wb, "quartic")[0]
bao_grid = np.array([bao_pred(OL) for OL in OL_GRID])
sne_grid = np.array([sne_mu(OL) for OL in OL_GRID])
cmb_grid = np.array([cmb_vec(j) for j in range(len(OL_GRID))])
i_true = int(np.argmin(np.abs(OL_GRID - OL_TRUE)))
bao_true = bao_grid[i_true]; sne_true = sne_grid[i_true]; cmb_true = cmb_grid[i_true]
print("precompute done", flush=True)

def quadmin(pred, obs, Cinv):
    A = float(pred @ Cinv @ pred); b = float(obs @ Cinv @ pred); c = float(obs @ Cinv @ obs)
    return c - b * b / A

def bao_lcdm_best(obs):
    def g(Om):
        E2 = Om * (1 + zg) ** 3 + OR_BAO * (1 + zg) ** 4 + 1.0 - Om - OR_BAO
        ic = cumulative_trapezoid(1.0 / np.sqrt(E2), zg, initial=0.0)
        Ei = interp1d(zg, E2, kind="cubic")
        DMi = interp1d(zg, ic, kind="cubic")
        DH_i = 1.0 / np.sqrt(Ei(zb))
        DM_i = DMi(zb)
        DV_i = (zb * DM_i ** 2 * DH_i) ** (1.0 / 3.0)
        pred = np.choose(kidx, [DM_i, DH_i, DV_i])
        return quadmin(pred, obs, Cbi)
    return minimize_scalar(g, bounds=(0.25, 0.45), method="bounded", options={"xatol": 1e-9}).fun

def sne_lcdm_best(obs):
    def g(Om):
        E2 = Om * (1 + sne.Z_FINE) ** 3 + sne.OMEGA_R_NU * (1 + sne.Z_FINE) ** 4 + 1.0 - Om - sne.OMEGA_R_NU
        ic = cumulative_trapezoid(1.0 / np.sqrt(E2), sne.Z_FINE, initial=0.0)
        DMi = interp1d(sne.Z_FINE, ic, kind="cubic", fill_value="extrapolate")
        DL = (1.0 + sne.DF["zHEL"]) * sne.C_H0 * DMi(sne.DF["zHD"])
        mu = 5.0 * np.log10(DL) + 25.0
        r = obs[zcut] - mu[zcut]
        Cr = Cinv_sne @ r
        return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv_sne))
    return minimize_scalar(g, bounds=(0.25, 0.45), method="bounded", options={"xatol": 1e-9}).fun

def _cmb_vec(j, h, wb):
    m = Model(theta_deg=THETA, Omega_L=float(OL_GRID[j]), Omega_r=9.08e-5)
    return cmb.vec_bimetric_m(m, h, wb, "quartic")[0]

# one-time: linearise v(h, wb) about the real-data fit; per realisation the
# chi2 is then an exactly-solvable 2x2 quadratic. bias vs the exact NM fit is
# checked once on the real data vector below.
CMB_V0, CMB_J = [], []
for j in range(len(OL_GRID)):
    v0 = _cmb_vec(j, cmb_h[j], cmb_wb[j])
    hh = cmb_h[j] + 1e-3
    ww = cmb_wb[j] + 1e-4
    Jh = (_cmb_vec(j, hh, cmb_wb[j]) - v0) / 1e-3
    Jw = (_cmb_vec(j, cmb_h[j], ww) - v0) / 1e-4
    CMB_V0.append(v0); CMB_J.append(np.column_stack([Jh, Jw]))
CMB_V0 = np.array(CMB_V0); CMB_J = np.array(CMB_J)
Ci_cmb = cmb.COV_INV
CMB_JtC = np.einsum("jka,kb->jab", CMB_J, Ci_cmb)          # (nOL, 2, 3)
CMB_A2 = np.einsum("jak,jkb->jab", CMB_JtC, CMB_J)          # (nOL, 2, 2)
CMB_BIAS = 0.0
for j in range(len(OL_GRID)):
    def g(p):
        h, wb = p
        if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
            return 1e10
        r = cmb_grid[j] - _cmb_vec(j, h, wb)
        return float(r @ Ci_cmb @ r)
    nm = minimize(g, [cmb_h[j], cmb_wb[j]], method="Nelder-Mead",
                  options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 400}).fun
    lin = float((cmb_grid[j] - CMB_V0[j]) @ Ci_cmb @ (cmb_grid[j] - CMB_V0[j]))
    CMB_BIAS = max(CMB_BIAS, abs(nm - lin))
print("CMB linearisation max |NM - linear| on real data: %.2e" % CMB_BIAS, flush=True)

def cmb_prof(j, obs):
    r0 = obs - CMB_V0[j]
    b = CMB_JtC[j] @ r0                                       # (2,)
    dp = np.linalg.solve(CMB_A2[j], b)
    r = r0 - CMB_J[j] @ dp
    return float(r @ Ci_cmb @ r)

def quad_vertex(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    i = int(np.argmin(y))
    i = min(max(i, 1), len(x) - 2)
    a, b, c = np.polyfit(x[i - 1:i + 2], y[i - 1:i + 2], 2)
    xv = -b / (2 * a)
    return float(xv), float(c - b * b / (4 * a))


def vertex_statistic_at_true(x, y, x_true):
    cs = CubicSpline(x, y)
    _, yv = quad_vertex(x, y)
    return float(cs(x_true) - yv)


CK = "../results/profile_calibration_checkpoint_%s.json" % TAG
rng = np.random.default_rng(SEED)
prof_argmin, cov1, cov384, dchi_true = [], [], [], []
_start = 0
if os.path.exists(CK):
    _c = json.load(open(CK))
    _start = _c["done"]
    prof_argmin = list(_c["prof_argmin"]); cov1 = [bool(x) for x in _c["cov1"]]
    cov384 = [bool(x) for x in _c["cov384"]]; dchi_true = list(_c["dchi_true"])
    rng.bit_generator.state = _c["rng_state"]
    print("resumed from %d" % _start, flush=True)
Lcmb = np.linalg.cholesky(cmb.COV)
for i in range(_start, N_MC):
    d_bao = bao_true + Lb @ rng.standard_normal(len(zb))
    d_sne = sne_true.copy()
    d_sne[zcut] = sne_true[zcut] + Lsne @ rng.standard_normal(int(zcut.sum()))
    d_cmb = cmb_true + Lcmb @ rng.standard_normal(3)
    bb = bao_lcdm_best(d_bao)
    sb = sne_lcdm_best(d_sne)
    prof = np.array([quadmin(bao_grid[j], d_bao, Cbi) - bb
                     + (lambda mu: float((d_sne[zcut] - mu[zcut]) @ Cinv_sne @ (d_sne[zcut] - mu[zcut])
                                         - np.sum(Cinv_sne @ (d_sne[zcut] - mu[zcut])) ** 2 / np.sum(Cinv_sne)) - sb)(sne_grid[j])
                     + cmb_prof(j, d_cmb) for j in range(len(OL_GRID))])
    xv, yv = quad_vertex(OL_GRID, prof)
    dchi_true_value = vertex_statistic_at_true(OL_GRID, prof, OL_TRUE)
    prof_argmin.append(float(xv))
    cov1.append(bool(dchi_true_value <= 1.0))
    cov384.append(bool(dchi_true_value <= 3.84))
    dchi_true.append(float(dchi_true_value))
    if (i + 1) % 20 == 0:
        print("%d/%d" % (i + 1, N_MC), flush=True)
        json.dump(dict(done=i + 1,
                       coverage_dchi1=round(float(np.mean(cov1)), 4),
                       coverage_dchi384=round(float(np.mean(cov384)), 4),
                       dchi_true=list(map(float, dchi_true)),
                       prof_argmin=list(map(float, prof_argmin)),
                       cov1=list(map(bool, cov1)), cov384=list(map(bool, cov384)),
                       rng_state=rng.bit_generator.state),
                  open(CK + ".tmp", "w"))
        os.replace(CK + ".tmp", CK)

# Always persist the final state, including segments whose length is not a
# multiple of 20. Without this, a short exact-stream segment has no resumable
# per-realisation arrays even when its computation completed successfully.
json.dump(dict(done=N_MC,
               coverage_dchi1=round(float(np.mean(cov1)), 4),
               coverage_dchi384=round(float(np.mean(cov384)), 4),
               dchi_true=list(map(float, dchi_true)),
               prof_argmin=list(map(float, prof_argmin)),
               cov1=list(map(bool, cov1)), cov384=list(map(bool, cov384)),
               rng_state=rng.bit_generator.state),
          open(CK + ".tmp", "w"))
os.replace(CK + ".tmp", CK)

out = dict(N=N_MC, seed=SEED, seeds=[SEED], theta=THETA, OL_true=OL_TRUE,
           grid_mode=GRID_MODE,
           statistic_mode="quadratic_vertex_plus_cubic_spline_at_true",
           covariance_scheme="canonical: noise ~ N(0, C), chi2 uses C^{-1} (as in probes/*.py)",
           cmb_linbias_real=float(CMB_BIAS),
           OL_grid=list(map(float, OL_GRID)), grid_spacing=0.004,
           coverage_dchi1=round(float(np.mean(cov1)), 4),
           coverage_dchi384=round(float(np.mean(cov384)), 4),
           dchi_true_quantiles={str(q): round(float(np.quantile(dchi_true, q)), 4) for q in (0.5, 0.68, 0.9, 0.95)},
           argmin_mean=round(float(np.mean(prof_argmin)), 5),
           argmin_quantiles={str(q): round(float(np.quantile(prof_argmin, q)), 4) for q in (0.16, 0.5, 0.84)})
_final = "../results/profile_calibration_%s.json" % TAG
with open(_final + ".tmp", "w") as _f:
    json.dump(out, _f, indent=1)
os.replace(_final + ".tmp", _final)
print(json.dumps(out, indent=1))
