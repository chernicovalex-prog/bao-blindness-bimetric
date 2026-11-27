"""Paper numbers not produced elsewhere in the pipeline (regenerated here).

Covers: Table of mass bounds (Sec. 6.4), BIC (Sec. 7.6), growth-integral tail and
scaling statements (Secs. 6.2, 6.4), unscreened-window growth exponent (Sec. 6.5),
the m_FP-independent theta-dependence of u*, SNe flat scans (Sec. 4.2), Table A1
(B parameters), Table A4 (D_M accumulation), Tables A2/A3 rows.
Writes results/extra_numbers.json and results/extra_numbers_log.txt; exit code 1
if any check fails.
"""
import sys, json, pathlib
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

SRC = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SRC)); sys.path.insert(0, str(SRC / "probes"))
from core import Model, E2_lcdm
import scalarsec, selfsim
import sne

LOG, OUT = [], {}
FAILED = [0]


def check(label, got, want, tol, rel=False):
    err = abs(got - want) / (abs(want) if rel else 1.0)
    ok = err <= tol
    LOG.append("%s  %s: got %.7g, paper %.7g (%s err %.2e <= %.1e)"
               % ("PASS" if ok else "FAIL", label, got, want, "rel" if rel else "abs", err, tol))
    FAILED[0] += (not ok)


def scheme_b(mf, Omega_r=9.08e-5):
    """Scheme B of the mass-bound table: Omega_m = 0.315, Omega_Lambda = 0.685 PRESCRIBED
    (the z = 0 quartic closure is not imposed here), benchmark theta, (alpha, beta)."""
    return Model(theta_deg=6.1, mFP_H0=float(mf), Omega_L=0.685, Omega_r=Omega_r), 0.315


def zstar_B(mf, Omega_r=9.08e-5):
    m, Om = scheme_b(mf, Omega_r)
    return selfsim.zstar_any(m, Om)


def mass_bounds():
    rows = {}
    H0_eV = 1.438e-33
    # X(0) = 0 edges and the minimum of X(0) in the low-mass window
    X0 = lambda mf: scalarsec.X_of_z(*scheme_b(mf)[:1], 0.0, 0.315)
    lo = brentq(X0, 1.05, 1.25, xtol=1e-9); hi = brentq(X0, 1.33, 1.6, xtol=1e-9)
    from scipy.optimize import minimize_scalar
    mn = minimize_scalar(X0, bounds=(1.15, 1.45), method="bounded", options={"xatol": 1e-9})
    check("Table 5 X(0)=0 lower edge m_FP/H0", lo, 1.1196, 5e-4)
    check("Table 5 X(0)=0 upper edge m_FP/H0", hi, 1.4363, 5e-4)
    check("Table 5 minimum of X(0)", mn.fun, -0.0760, 5e-4)
    check("Table 5 location of the minimum", mn.x, 1.291, 2e-3)
    rows["X0_lo"], rows["X0_hi"], rows["X0_min"] = lo, hi, [mn.x, mn.fun]
    for z_t, want, eV in ((2.33, 5.937, 8.5e-33), (4.2, 11.437, 1.6e-32)):
        mf = brentq(lambda x: zstar_B(x) - z_t, 3.0, 40.0, xtol=1e-8)
        check("Table 5 m_FP/H0 for z* = %g" % z_t, mf, want, 5e-4)
        check("Table 5 m_FP[eV] for z* = %g" % z_t, mf * H0_eV, eV, 0.04, rel=True)
        rows["z*=%g" % z_t] = mf
    mrec = brentq(lambda x: zstar_B(x) - 1090.0, 2e4, 8e4, xtol=1e-3)
    mrec0 = brentq(lambda x: zstar_B(x, 0.0) - 1090.0, 2e4, 8e4, xtol=1e-3)
    check("Table 5 m_FP/H0 for z* = 1090", mrec, 41153, 5e-4, rel=True)
    check("Table 5 m_FP/H0 for z* = 1090 without radiation", mrec0, 34556, 5e-4, rel=True)
    check("Table 5 m_FP[eV] for z* = 1090", mrec * H0_eV, 5.9e-29, 0.02, rel=True)
    check("Sec 6.4 radiation shift of the recombination bound (19%)", 100 * (mrec / mrec0 - 1), 19.1, 0.5)
    rows["z*=1090"], rows["z*=1090_noRad"] = mrec, mrec0
    OUT["mass_bounds"] = rows


def tails_and_scaling():
    m = Model(); Om = m.solve_Om0()
    zs = scalarsec.find_zstar(m, Om)
    I1, _ = scalarsec.growth_integral(m, Om, zs)
    # tail beyond z = 1090 with the N_eff radiation
    f = lambda z: np.sqrt(-scalarsec.X_of_z(m, z, Om)) / np.sqrt(m.E2_bim(z, Om))
    zmax = 5e4        # solve_y is bracketed below ~1e5; the rest is analytic (E ~ z^2, X -> const)
    tail, _ = quad(f, 1090.0, zmax, epsabs=1e-12, epsrel=1e-10, limit=400)
    tail += np.sqrt(-scalarsec.X_of_z(m, zmax, Om)) / np.sqrt(m.E2_bim(zmax, Om)) * zmax
    check("Sec 6.2 tail of I1 beyond z = 1090 (to infinity; the draft quoted 0.067)", tail, 0.0734, 5e-4)
    check("Sec 6.2 tail fraction of I1 (draft: 7.1%)", 100 * tail / I1, 7.77, 0.05)
    lg = lambda I: I * 0.05 / (100 * 0.7 / 299792.458) / np.log(10)
    check("Sec 6.2 log10 Gamma(k=0.05) with the tail (draft: 94.1)", lg(I1 + tail), 94.70, 0.05)
    check("Sec 6.4 radiation-free tail 2/sqrt(Om(1+z_rec))", 2 / np.sqrt(Om * 1091.0), 0.1098, 5e-4)
    # unscreened window growth exponent
    zV = brentq(lambda z: np.sqrt(m.E2_bim(z, Om)) - m.mFP_H0, 0.05, 3000)
    Iu, _ = quad(f, zs, zV, epsabs=1e-12, epsrel=1e-12, limit=300)
    for h, want in ((0.700, 10.8), (0.678, 11.2)):
        H0 = 100 * h / 299792.458
        check("Sec 6.5 unscreened-window log10 Gamma (h=%.3f)" % h, 0.05 / H0 * Iu / np.log(10), want, 0.06)
    # I1 * m^{1/3} variation along the m_FP scan
    S = json.load(open(SRC.parent / "results" / "scalar_mfp_scan.json"))
    mm = np.array([r["mFP_H0"] for r in S]); ii = np.array([r["I1"] for r in S])
    q = ii * mm ** (1.0 / 3.0)
    check("Sec 6.4 I1 m^(1/3) at m=1.5", q[0], 2.87, 6e-3)
    check("Sec 6.4 I1 m^(1/3) at m=200", q[-1], 2.44, 6e-3)
    check("Sec 6.4 I1 m^(1/3) variation over the full scan (15%)", 100 * (q[0] - q[-1]) / q[0], 15.0, 0.6)
    sel = mm >= 10
    check("Sec 6.4 I1 m^(1/3) variation above m=10 (14%)", 100 * (q[sel].max() - q[sel].min()) / q[sel].max(), 14.0, 1.0)
    check("Sec 6.2 I1 at m=1.5", ii[0], 2.50, 6e-3)
    check("Sec 6.2 I1 at m=200", ii[-1], 0.42, 6e-3)
    un = np.array([r["unscreened_pct"] for r in S])
    check("Sec 6.2 unscreened share min", un.min(), 9.7, 0.05)
    check("Sec 6.2 unscreened share max", un.max(), 13.9, 0.05)
    check("Sec 6.2 z* at m=1.5", S[0]["zstar"], 0.06, 6e-3)
    check("Sec 6.2 z* at m=200", S[-1]["zstar"], 34.0, 0.5)
    # theta dependence of u* and w* = t u*
    us = {}
    for th in (4.0, 6.1, 8.0):
        t = np.tan(np.radians(th)) ** 2
        u = selfsim.find_ustar(t, 0.0, 3.0, 0.0)
        us[th] = (u, t * u)
    OUT["ustar_theta"] = us
    check("Sec 6.4 u*(4 deg) ~ 70", us[4.0][0], 70.0, 4.0)
    check("Sec 6.4 u*(8 deg) ~ 17.5", us[8.0][0], 17.5, 1.0)
    check("Sec 6.4 w* spread over 4-8 deg (1.7%)", 100 * (max(v[1] for v in us.values()) / min(v[1] for v in us.values()) - 1), 1.7, 0.2)
    OUT["I1_tail"] = tail


def bic():
    check("Sec 7.6 Delta BIC (BAO alone, dk=-1)", 0.179 - np.log(13), -2.39, 6e-3)
    check("Sec 7.6 Delta BIC (four parameters, dk=+3)", 0.179 + 3 * np.log(13), 7.87, 6e-3)


def sne_scans():
    sne.load_pantheon()
    for full, name, want in ((False, "diag", 7.315), (True, "full", 2.936)):
        sc = sne.scan_grid(full=full)
        vals = [v[0] for k, v in sc.items() if v is not None]
        check("Sec 4.2 SNe %s: (alpha,beta) + m_FP scans flat (n=%d, min)" % (name, len(vals)), min(vals), want, 2e-3)
        check("Sec 4.2 SNe %s: (alpha,beta) + m_FP scans flat (max)" % name, max(vals), want, 2e-3)
        OUT["sne_scan_" + name] = [min(vals), max(vals), len(vals)]
    # theta scan at fixed Omega_Lambda = 0.699 with the self-consistent quartic Omega_m
    # (Sec. 4.2, Table 2): flat because the late-time background depends on cos^2(theta) Omega_m
    zcut = sne.DF["zHD"] > 0.01
    cov = sne.DF["COV_FULL"][np.ix_(zcut, zcut)]
    Ci = np.linalg.inv(cov)
    def c2f(mu):
        r = sne.DF["m_b_corr"][zcut] - mu
        Cr = Ci @ r
        return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Ci))
    from scipy.optimize import minimize_scalar
    kw = dict(bounds=(0.2, 0.5), method="bounded", options={"xatol": 1e-12})
    ref_d = minimize_scalar(lambda Om: sne.chi2_diag(sne.lcdm_mu(Om, zcut), zcut, 0.0), **kw).fun
    ref_f = minimize_scalar(lambda Om: c2f(sne.lcdm_mu(Om, zcut)), **kw).fun
    dd, df = [], []
    for th in (6.1, 10.0, 15.0, 20.0, 25.0, 31.5):
        mu = sne.bimetric_mu_at(theta_deg=th)[zcut]
        dd.append(sne.chi2_diag(mu, zcut, 0.0) - ref_d)
        df.append(c2f(mu) - ref_f)
    check("Table 2 SNe diag theta scan 6.1-31.5 (min)", min(dd), 7.315, 2e-3)
    check("Table 2 SNe diag theta scan 6.1-31.5 (max)", max(dd), 7.320, 2e-3)
    check("Table 2 SNe full theta scan 6.1-31.5 (min)", min(df), 2.936, 2e-3)
    check("Table 2 SNe full theta scan 6.1-31.5 (max)", max(df), 2.939, 2e-3)
    OUT["sne_theta_scan"] = dict(theta=[6.1, 10.0, 15.0, 20.0, 25.0, 31.5], diag=dd, full=df)


def tables_A():
    m = Model(); Om = m.solve_Om0(); B = m.B_params()
    for b, w in zip(B, (-40.248, 28.230, -21.173, 21.173, -28.206)):
        check("Table A1 B-parameter", b, w, 6e-4)
    check("Table A1 B1^2", B[1] ** 2, 796.94, 5e-4)
    # Table A2 (selected rows) and A3
    A2 = {0.1: (0.999350, 1.099643, -0.1036), 1.0: (0.996093, 3.107718, -0.7704), 2.33: (0.981968, 11.821478, -1.0368),
          10.0: (0.556527, 403.379176, -0.9191)}
    for z, (y, e2, d) in A2.items():
        check("Table A2 y(z=%g)" % z, m.solve_y(z, Om), y, 2e-6)
        check("Table A2 E2_bim(z=%g)" % z, m.E2_bim(z, Om), e2, 2e-6, rel=True)
        check("Table A2 dE2[%%](z=%g)" % z, 100 * (m.E2_bim(z, Om) / E2_lcdm(z, Om) - 1), d, 2e-4)
    for z, w in ((0.5, -1.0169), (1.0, -1.0411), (2.0, -1.1536), (3.0, -1.4610)):
        check("Table A3 w_DE(z=%g)" % z, m.w_DE(z, Om), w, 3e-4)
    # Table A4: D_M(z_rec) accumulation
    zr = 1090.0
    dm = lambda a, b: quad(lambda z: 1 / np.sqrt(m.E2_bim(z, Om)), a, b, epsabs=1e-11, epsrel=1e-11, limit=300)[0]
    tot = dm(0, zr)
    for zc, w in ((1, 24.32), (5, 57.19), (10, 69.42), (12, 72.19), (50, 87.85), (500, 98.51)):
        check("Table A4 share of D_M below z=%g [%%]" % zc, 100 * dm(0, zc) / tot, w, 0.02)


def main():
    for fn in (mass_bounds, tails_and_scaling, bic, sne_scans, tables_A):
        print("running", fn.__name__, flush=True)
        fn()
    (SRC.parent / "results" / "extra_numbers.json").write_text(json.dumps(OUT, indent=1, default=float))
    txt = "\n".join(LOG) + "\n" + ("ALL EXTRA CHECKS PASSED (%d)\n" % len(LOG) if not FAILED[0]
                                 else "%d EXTRA CHECKS FAILED\n" % FAILED[0])
    (SRC.parent / "results" / "extra_numbers_log.txt").write_text(txt)
    print(txt)
    return FAILED[0]


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
