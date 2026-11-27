"""Self-checks for bimetric_r28. Run after the probe scripts; every check
raises on failure. These are the canonical internal consistency gates:
E^2(0) closure, B-parameter round-trip, scalar-sector canon, and the
probe-level canon values used in the paper.
"""
import numpy as np
import json
import pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))

from core import Model
import scalarsec
import selfsim
import bao, sne, cmb


def check(name, cond):
    print(("PASS  " if cond else "FAIL  ") + name)
    if not cond:
        raise SystemExit("check failed: " + name)


def main():
    m = Model()
    # 1. z = 0 closure: E^2(0) = 1 to solver precision
    Om = m.solve_Om0()
    resid = abs(m.E2_bim(0.0, Om) - 1.0)
    check("E2(0)-1 residual < 1e-12 (got %.2e)" % resid, resid < 1e-12)

    # 2. B-parameter round-trip (H&M 2021 Eq. 3.3 vs 3.4)
    B = m.B_params()
    mfp2 = (B[1] + 2 * B[2] + B[3]) / m.sin2
    OL = B[0] / 3 + B[1] + B[2] + B[3] / 3
    alpha = -(B[2] + B[3]) / (B[1] + 2 * B[2] + B[3])
    beta = B[3] / (B[1] + 2 * B[2] + B[3])
    rt = max(abs(mfp2 - m.mFP2), abs(OL - m.Omega_L),
             abs(alpha - m.alpha), abs(beta - m.beta))
    check("B round-trip < 1e-10 (got %.2e)" % rt, rt < 1e-10)

    # 3. scalar-sector canon (self-consistent closure; legacy Omega_m0 = 0.304
    #    convention documented in the paper: z* = 7.880536, I1 = 0.945210110)
    s = scalarsec.run_scalar()
    check("z* == 7.877101 (self-consistent closure)",
          abs(s["zstar"] - 7.877101) < 1e-5)
    check("I1 == 0.944870550", abs(s["I1"] - 0.944870550) < 1e-8)
    s_legacy = scalarsec.run_scalar(Omega_m0=0.304)
    check("legacy z* == 7.880536", abs(s_legacy["zstar"] - 7.880536) < 1e-6)
    check("log10 Gamma(0.05,0.7) == 87.87",
          abs(s["log10_gamma"][(0.05, 0.7)] - 87.87) < 0.01)

    # 4. probe-level canon
    bao.load_desi()
    bl = bao.baseline()
    check("BAO |dchi2 - 0.1792| < 5e-4 (got %.4f)" % bl["dchi2"],
          abs(bl["dchi2"] - 0.1792) < 5e-4)

    sne.load_pantheon()
    rs = sne.run(0.01)
    check("SNe full dchi2 == 2.94 (zHD primary) (got %.4f)" % rs["full"]["dchi2"],
          abs(rs["full"]["dchi2"] - 2.9357) < 0.01)
    check("SNe diag dchi2 == 7.32 (zHD primary) (got %.4f)" % rs["diag"]["dchi2"],
          abs(rs["diag"]["dchi2"] - 7.3152) < 0.02)
    check("SNe full Om == 0.331422 (got %.6f)" % rs["full"]["Om_lcdm"],
          abs(rs["full"]["Om_lcdm"] - 0.331422) < 5e-5)
    check("SNe diag Om == 0.350241 (got %.6f)" % rs["diag"]["Om_lcdm"],
          abs(rs["diag"]["Om_lcdm"] - 0.350241) < 5e-5)
    check("SNe full pure correction == 0.64 (got %.4f)" % rs["full"]["pure"],
          abs(rs["full"]["pure"] - 0.6370) < 0.02)
    # Cross-convention values quoted in Sec. 3: each covariance uses its own
    # LCDM best-fit Omega_m, not the other covariance's best fit.
    _zcut = sne.DF["zHD"] > 0.01
    from scipy.integrate import cumulative_trapezoid as _ctrap
    from scipy.interpolate import interp1d as _interp1d
    def _mu_at_om(om0):
        _e2 = np.array([m.E2_bim(z, om0) for z in sne.Z_FINE])
        _ic = _ctrap(1.0 / np.sqrt(_e2), sne.Z_FINE, initial=0.0)
        _dmi = _interp1d(sne.Z_FINE, _ic, kind="cubic", fill_value="extrapolate")
        _dl = ((1.0 + sne.DF["zHEL"][_zcut]) * sne.C_H0
               * _dmi(sne.DF["zHD"][_zcut]))
        return 5.0 * np.log10(_dl) + 25.0
    _cross_full = (sne.chi2_full(_mu_at_om(rs["full"]["Om_lcdm"]), _zcut)
                   - rs["full"]["chi2_min"])
    _cross_diag = (sne.chi2_diag(_mu_at_om(rs["diag"]["Om_lcdm"]), _zcut)
                   - rs["diag"]["chi2_min"])
    check("SNe cross-convention full/diag == +0.4645/+1.0433 (got %.4f/%.4f)"
          % (_cross_full, _cross_diag),
          abs(_cross_full - 0.4645375) < 1e-3 and abs(_cross_diag - 1.0433292) < 1e-3)

    c2, pref, scan = cmb.theta_scan(thetas=(6.1,))
    base = scan[6.1]
    check("CMB baseline chi2 == 0.0003 +/- 0.005 (got %.5f)" % base["chi2"],
          abs(base["chi2"] - 0.00035) < 5e-3)
    from scipy.stats import chi2 as _chi2_dist
    _p_cmb = float(_chi2_dist.sf(base["chi2"], 1))
    check("CMB baseline one-DOF tail p == 0.990 (got %.6f)" % _p_cmb,
          abs(_p_cmb - 0.9895891) < 1e-5)

    # 4b. D_M integrator gate (review round 6): the spline+quad path inside
    #     vec_bimetric_m must match a direct quad with the exact solve_y to
    #     1e-5 relative. This guards the pre-8 uniform-grid trapezoid +
    #     cubic-interpolation D_M, which overestimated D_M(z_rec) by ~0.1%
    #     even for pure LambdaCDM and mixed schemes with the reference.
    from scipy.integrate import quad as _q
    mb = Model(theta_deg=6.1, Omega_L=0.699)
    hb, wbb = 0.6864, 0.0224
    Omb = cmb.solve_Om0_nu(mb, hb)
    zsb = cmb.z_star_hs(wbb, Omb * hb ** 2 + cmb.OMEGA_NU, hb)
    yspl = cmb._y_spline(mb, Omb, hb)

    def _E2_gate(z, use_spline):
        dy = 1.0 - (float(yspl(z)) if use_spline else cmb.solve_y_cmb(mb, z, Omb, hb))
        return (cmb.omega_ng(z, Omb, hb)
                + mb.Omega_L - mb.sin2 * mb.mFP2 * dy
                * (1 + mb.alpha * dy + mb.beta / 3.0 * dy ** 2))

    dm_spl, _ = _q(lambda z: 1.0 / np.sqrt(_E2_gate(z, True)), 0, zsb,
                   epsabs=1e-12, epsrel=1e-12, limit=200)
    dm_ex, _ = _q(lambda z: 1.0 / np.sqrt(_E2_gate(z, False)), 0, zsb,
                  epsabs=1e-12, epsrel=1e-12, limit=200)
    check("CMB D_M: spline+quad vs exact quad < 1e-5 rel (got %.1e)"
          % abs(dm_spl / dm_ex - 1), abs(dm_spl / dm_ex - 1) < 1e-5)

    # 5. Planck regression: the fitted LCDM reference must reproduce the
    #    compressed data vector (guards the neutrino treatment and the
    #    integration scheme; convention-proof: evaluated at the pipeline's
    #    own best fit, not at frozen historical parameters)
    c2r, pref = cmb.fit_lcdm_reference()
    vp = cmb.vec_lcdm(*pref)
    check("Planck regression chi2_ref ~ 0 (got %.1e)" % c2r, c2r < 1e-10)
    check("Planck regression |R - 1.7504| < 1e-4 (got %.4f)" % vp[0],
          abs(vp[0] - 1.7504) < 1e-4)
    check("Planck regression |l_a - 301.77| < 1e-2 (got %.3f)" % vp[1],
          abs(vp[1] - 301.77) < 1e-2)

    # review-11 gates: convention self-consistency of the reference
    # (a) bimetric theta -> 0 limit must coincide with vec_lcdm (both now
    #     carry omega_nu in R and in z_rec): theta = 0.1 deg, matched params
    m_gr = Model(theta_deg=0.1, Omega_L=0.699, Omega_r=9.08e-5)
    om_gr = cmb.solve_Om0_nu(m_gr, 0.675822)
    vec_gr = cmb.vec_bimetric_m(m_gr, 0.675822, 0.022371, "quartic")[0]
    vec_lc = cmb.vec_lcdm(om_gr * 0.675822 ** 2, 0.675822, 0.022371)
    rel_r = abs(vec_gr[0] / vec_lc[0] - 1.0)
    rel_la = abs(vec_gr[1] / vec_lc[1] - 1.0)
    check("theta->0 limit == vec_lcdm: |dR/R| = %.2e, |dla/la| = %.2e (< 1e-4)"
          % (rel_r, rel_la), rel_r < 1e-4 and rel_la < 1e-4)
    # (b) late-time bimetric E^2 coincides with LambdaCDM(Omega_eff) to 3e-5
    #     (the dark-degeneracy content: Omega_eff = cos^2(theta) Omega_m)
    zs_g = np.linspace(0.01, 2.33, 240)
    om_b = m.solve_Om0()
    Om_eff = (1.0 - m.sin2) * om_b
    dev = max(abs(m.E2_bim(z, om_b) / (m.Omega_L + Om_eff * (1 + z) ** 3
                                       + m.Omega_r * (1 + z) ** 4) - 1.0)
              for z in zs_g)
    check("late-time E^2 == LambdaCDM(Omega_eff = cos^2 theta Om_m): max dev %.2e (< 4e-5)"
          % dev, dev < 4e-5)

    # 6. printed quartic constants (Appendix B) reproduce the closure
    t = m.tan2
    c4 = t * B[3]
    c3 = 3 * t * B[2] - B[4]
    c2 = 3 * (t * B[1] - B[3])
    Om0_full = Om + m.Omega_r
    c1 = t * (B[0] + 3 * Om0_full) - 3 * B[2]
    PRINTED = dict(c4=0.241813, c3=27.480744, c2=-62.550547, c1=63.068553, t=0.0114210,
                   slope=0.034263, Om0=0.3044455)
    for name, val in PRINTED.items():
        got = dict(c4=c4, c3=c3, c2=c2, c1=c1, t=t, slope=3.0 * t, Om0=Om0_full)[name]
        check("App-B constant %s (printed %.6g == code %.6g)" % (name, val, got),
              abs(val - got) < 5e-7)
    from scipy.optimize import brentq as _bq
    def E2_from_printed(z, Om):
        Omz = Om * (1 + z) ** 3 + m.Omega_r * (1 + z) ** 4
        c1z = t * (B[0] + 3 * Omz) - 3 * B[2]
        pol = lambda y: ((((PRINTED["c4"] * y + PRINTED["c3"]) * y + PRINTED["c2"]) * y + c1z) * y - B[1])
        y = _bq(pol, 1e-12, 1.0)
        dy = 1 - y
        ODE = (m.Omega_L - m.sin2 * m.mFP2 * dy
               * (1 + m.alpha * dy + m.beta / 3.0 * dy ** 2))
        return Om * (1 + z) ** 3 + m.Omega_r * (1 + z) ** 4 + ODE
    resid = abs(E2_from_printed(0.0, Om) - 1.0)
    check("closure from printed App-B constants < 5e-6 (got %.2e)" % resid,
          resid < 5e-6)

    # 7. Gamma table reconstruction (single I1)
    I1v = 0.944870550
    for (k, h), expect in {(0.001, 0.700): 1.76, (0.05, 0.700): 87.87,
                           (0.10, 0.700): 175.74, (0.001, 0.678): 1.81,
                           (0.05, 0.678): 90.72, (0.10, 0.678): 181.45}.items():
        H0v = 100.0 * h / 299792.458
        got = (k / H0v) * I1v / np.log(10.0)
        check("Gamma table k=%g h=%.3f (paper %.2f)" % (k, h, expect),
              abs(got - expect) < 0.005 + 1e-9)


    # 8. Self-similar reduction of the instability boundary (Sec. 6.4) and the
    #    ghost-free branch. The reduced quartic is m_FP-free except for the
    #    exact linear eps-terms; the boundary is universal in u and the branch
    #    constants (u*, y*, c, kappa) follow from the reduced problem alone.
    ss = selfsim.branch_constants()
    check("selfsim u* == 30.337301 (reduced problem; got %.6f)" % ss["ustar"],
          abs(ss["ustar"] - 30.337301) < 1e-5)
    check("selfsim y* == 0.705465 (got %.6f)" % ss["ystar"],
          abs(ss["ystar"] - 0.705465) < 1e-5)
    check("selfsim c == 0.582199 (got %.6f)" % ss["c"],
          abs(ss["c"] - 0.582199) < 1e-5)
    check("selfsim kappa == 0.374329 (got %.6f)" % ss["kappa"],
          abs(ss["kappa"] - 0.374329) < 1e-5)
    slope = selfsim.eps_linearity_slope()
    check("selfsim d u*/deps == -0.626809 (got %.6f)" % slope,
          abs(slope - (-0.626809)) < 1e-5)

    # exact eps-shifts of the reduced-quartic coefficients (Lemma 1):
    # da3 = -3 t eps, da1 = +3 t eps, all others zero, no O(eps^2)
    t_ = ss["tan2"]
    for eps_t in (0.1, 1.7, 27.5):
        d = selfsim.eps_shift_coeffs(t_, eps_t)
        check("selfsim eps-shift (eps=%.2f): da3 = -3t eps, da1 = +3t eps, rest 0"
              % eps_t,
              abs(d[1] + 3 * t_ * eps_t) < 1e-14 and abs(d[3] - 3 * t_ * eps_t) < 1e-14
              and abs(d[0]) < 1e-14 and abs(d[2]) < 1e-14 and abs(d[4]) < 1e-14)

    # ghost-free branch (Könnig 2015, Sec. III A): the Higuchi bound on the
    # finite branch is dy/dN >= 0 everywhere. dy/dN(0) = 1.46e-3, the
    # minimum over 0 <= z <= 1090 is 5.1e-6 at z = 1090 (dense grid).
    gf = selfsim.ghost_free_scan()
    check("ghost-free: dy/dN(0) == 1.4648e-3 (got %.4e)" % gf["dydN"][0],
          abs(gf["dydN"][0] - 1.46475715e-3) < 1e-8)
    check("ghost-free: min dy/dN over [0,1090] > 0 (got %.3e at z = %.0f)"
          % (gf["min"], gf["z_at_min"]), gf["min"] > 0.0)

    # deep self-similar crosscheck: the full Model at m_FP = 2000, Omega_r = 0
    # must reproduce the reduced-problem u* and the plateau r = H(z*)/m_FP = c
    xc = selfsim.full_model_crosscheck()
    check("selfsim full-vs-reduced u* at m=2000 (%.6f vs %.6f)"
          % (xc["ustar_full"], xc["ustar_reduced"]),
          abs(xc["ustar_full"] - xc["ustar_reduced"]) < 2e-4)
    check("selfsim plateau r == c at m=2000 (%.8f vs %.8f)" % (xc["r_full"], xc["c"]),
          abs(xc["r_full"] - xc["c"]) < 1e-6)
    pl = selfsim.plateau_scan()
    rvals = np.array([r for _, r, _ in pl])
    check("selfsim plateau flat on m=500-8000 (max|r-c| = %.2e)"
          % np.max(np.abs(rvals - ss["c"])), np.max(np.abs(rvals - ss["c"])) < 2e-6)

    # z* scaling: ln(1+z*) vs ln m_FP has slope 2/3 in the matter-dominated
    # self-similar regime
    sl = selfsim.zstar_scaling()
    check("selfsim z* scaling slope ~ 2/3 (got %.4f)" % sl, abs(sl - 2.0 / 3.0) < 0.0025)

    # H(z*) law at the benchmark: sqrt(c^2 m_FP^2 + kappa Omega_Lambda) vs full
    hl = selfsim.Hstar_law_check()
    check("selfsim H(z*) law vs full model (rel %+.5f, < 1e-3)" % hl["rel"],
          abs(hl["rel"]) < 1e-3)

    # analytic asymptote: X(u) -> -1 exactly as u -> infinity
    # (y ~ (1+2a+b)/(3 t u) => dy/dN = 3 y, T1 -> 1, T2 -> 3)
    xa = selfsim.X_asymptote()
    check("selfsim X(u=1e8) == -1 asymptote (got %.8f)" % xa["X"],
          abs(xa["X"] + 1.0) < 1e-4)
    deep = selfsim.X_deep_negative()
    check("selfsim X(u) < 0 on the deep grid 4e3-1e8",
          all(X < 0.0 for _, X in deep))

    # master equation: Omega_m(1+z*)^3 = sin^2 m_FP^2 u*(eps) [1 + O(rho)]
    # reproduces the full model across the self-similar regime
    mc = selfsim.zstar_master_check()
    devmax = max(r["dev"] for r in mc)
    check("selfsim master equation max dev over m=3-200 = %.4f (< 6e-3)" % devmax,
          devmax < 6e-3)

    # k-dilemma map (Sec. 6.2, Table 7): horizon entry z_h(k) (k = aH =
    # H0 E/(1+z), definition-checked) vs log10 Gamma(k)
    km = scalarsec.k_map()
    row05 = [r for r in km["rows"] if abs(r["k"] - 0.05) < 1e-12][0]
    check("k-map z_h(0.05) == 20857.86 (got %.4f)" % row05["z_h"],
          abs(row05["z_h"] - 20857.8615) < 0.5)
    check("k-map k(z_h=1090) == 4.903e-3 (got %.4e), k(Gamma=1) == 5.690e-4 (got %.4e),"
          " z_h(Gamma=1) == 18.5 (got %.4f)"
          % (km["k_at_zrec"], km["k_gamma1"], km["z_h_gamma1"]),
          abs(km["k_at_zrec"] - 4.9029e-3) < 1e-5 and
          abs(km["k_gamma1"] - 5.6901107e-4) < 1e-8 and
          abs(km["z_h_gamma1"] - 18.507) < 0.1)
    check("k-map log10Gamma(5.0) == 8787.18 (got %.2f)" %
          [r for r in km["rows"] if r["k"] == 5.0][0]["log10Gamma"],
          abs([r for r in km["rows"] if r["k"] == 5.0][0]["log10Gamma"] - 8787.1753) < 1e-2)
    # definition gate: k = aH must hold at the reported z_h (relative 1e-6)
    mdef = Model(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0,
                 Omega_L=0.699, Omega_r=9.08e-5)
    omdef = mdef.solve_Om0()
    for kk in (0.001, 0.05, 1.0):
        zh = scalarsec.z_horizon_entry(mdef, omdef, kk)
        lhs = kk
        rhs = (100.0 * 0.7 / 299792.458) * np.sqrt(
            scalarsec.E2_horizon(mdef, omdef, zh)) / (1.0 + zh)
        check("k-map definition k = aH at z_h(k=%.3g): rel %.2e < 1e-6" % (kk, abs(lhs/rhs-1)),
              abs(lhs / rhs - 1.0) < 1e-6)
    # rational boundary function X(y) = -Q(y)/[(3y^2-6y+4)A(y)^2]:
    # both denominator factors positive on (0,1); Sturm count gives a single
    # root of the degree-10 numerator Q at y* (numeric Sturm, benchmark t)
    ys = (0.3, 0.5, 0.9)
    t6 = selfsim.branch_constants()["tan2"]
    xdiff = max(abs(selfsim.X_rational(6.1, y)
                    - selfsim.X_of_u(t6, 0.0, 3.0, selfsim.u_of_y(6.1, y), 0.0))
                for y in ys)
    check("selfsim rational X(y) == X_of_u(u(y)) (max diff %.2e)" % xdiff, xdiff < 1e-10)
    yg2 = np.linspace(1e-9, 1.0 - 1e-9, 20001)
    A6 = (9.0 * t6 * yg2 ** 4 - 18.0 * t6 * yg2 ** 3 + 12.0 * t6 * yg2 ** 2
          + 8.0 * yg2 ** 3 - 9.0 * yg2 ** 2 + 4.0)
    check("selfsim X(y) denominator positive on (0,1): min A = %.4f, min 3y^2-6y+4 = %.4f"
          % (A6.min(), (3 * yg2 ** 2 - 6 * yg2 + 4).min()),
          A6.min() > 0.0)
    nsturm = selfsim.sturm_root_count(selfsim.q_coeffs(6.1), 1e-12, 1.0 - 1e-12)
    y_sturm = [r.real for r in np.roots(selfsim.q_coeffs(6.1))
               if abs(r.imag) < 1e-8 and 0.0 < r.real < 1.0]
    check("selfsim Sturm: Q has a single root in (0,1) (count %d), at y = %.9f == y*"
          % (nsturm, y_sturm[0] if y_sturm else -1),
          nsturm == 1 and len(y_sturm) == 1 and abs(y_sturm[0] - ss["ystar"]) < 1e-8)

    # zCMB purge guard: the primary SNe canon is the zHD-selection value; a
    # primary path silently reverting to zCMB selection would shift dchi2 by
    # ~3.6 and must fail here, not propagate (bug class of the shared-fit fix)
    _r_zc = sne.run_zcut_on("zCMB", 0.01)
    check("zCMB purge guard: primary SNe dchi2 == zHD selection (%.4f == 2.9357);"
          " selection-only zCMB cross-check == %.4f (== 3.10, Sec. 4.2; the 6.50"
          " variant additionally changes the integral and is quoted, not recomputed)"
          % (rs["full"]["dchi2"], _r_zc["full"]["dchi2"]),
          abs(rs["full"]["dchi2"] - 2.9357) < 0.01 and abs(_r_zc["full"]["dchi2"] - 3.10) < 0.02)

    # independent high-precision background implementation:
    # self-contained mpmath/polyroots path from the B-parameters agrees with
    # the canonical Brent solver to machine precision (and works at large m_FP)
    import independent_bg as _ibg
    _xi = _ibg.crosscheck(theta_deg=6.1, mFP_H0=25.0)
    check("independent bg (mpmath polyroots): benchmark max rel dev %.2e (<1e-10),"
          " closure |dOm| = %.2e" % (_xi["max_rel_dev"], abs(_xi["Om_indep"] - _xi["Om_canon"])),
          _xi["max_rel_dev"] < 1e-10 and abs(_xi["Om_indep"] - _xi["Om_canon"]) < 1e-8)
    _x2 = _ibg.crosscheck(theta_deg=30.0, mFP_H0=200.0)
    check("independent bg at (30 deg, m=200): max rel dev %.2e (<1e-8)" % _x2["max_rel_dev"],
          _x2["max_rel_dev"] < 1e-8)

    import os as _os2
    # MC calibration with canonical covariance and exact-stream regeneration:
    # fine coverage 0.683/0.9535, offset 0.680/0.9530.
    _pc = json.load(open(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                        "profile_calibration.json")))
    check("profile calibration (fine grid): N=%d, coverage(dchi=1)=%.3f (==0.683), coverage(dchi=3.84)=%.3f (==0.9535)"
          % (_pc["N"], _pc["coverage_dchi1"], _pc["coverage_dchi384"]),
          _pc["N"] == 2000 and abs(_pc["coverage_dchi1"] - 0.683) < 0.03
          and abs(_pc["coverage_dchi384"] - 0.9535) < 0.02
          and _pc.get("exact_stream") is True and _pc["seeds"] == [20261118]
          and "canonical" in _pc["covariance_scheme"]
          and float(_pc["cmb_linbias_real"]) == 0.0)
    # Synthetic CMB spot-check of the linearised profiling
    _lb = json.load(open(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                        "cmb_linbias_check.json")))
    check("CMB linearisation spot-check on %d synthetic realisations: max bias %.2e (<1e-3)"
          % (_lb["n_rep"], _lb["max_bias"]), _lb["max_bias"] < 1e-3)
    # Offset-grid robustness (true anchor NOT a grid point)
    _og = json.load(open(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                        "profile_calibration_offset.json")))
    check("profile calibration (offset grid +0.002): N=%d, coverage(dchi=1)=%.3f, (dchi=3.84)=%.3f (consistent with fine grid)"
          % (_og["N"], _og["coverage_dchi1"], _og["coverage_dchi384"]),
          _og["grid_mode"] == "offset" and round(0.699, 4) not in [round(x, 4) for x in _og["OL_grid"]]
          and abs(_og["coverage_dchi1"] - _pc["coverage_dchi1"]) < 0.04
          and abs(_og["coverage_dchi384"] - _pc["coverage_dchi384"]) < 0.03)

    # shared-Omega_m reference consistency (bug-fixed 2026-11-15): a coupled
    # fit can never lie below the sum of uncoupled minima (cost >= 0)
    import os as _os2
    _ls = json.load(open(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                        "anchor_scan.json")))["lcdm_shared"]
    check("lcdm_shared consistency: cost = %.4f >= 0, Om = %.4f, joint = %.2f > sum of minima %.2f"
          % (_ls["cost"], _ls["Om"], _ls["joint"], _ls["chi2_sum_min"]),
          _ls["cost"] >= 0.0 and _ls["joint"] > _ls["chi2_sum_min"] - 1e-9)

    # postponed batch: numbers that were quoted but ungated (all verified)
    # (a) synthetic BAO null, fixed seed (Sec. 3.2, p = 0.682)
    bao.load_desi()
    p_null = bao.synthetic_null_refit(n=2000, seed=20261014)
    check("synthetic BAO null p == 0.682 (seed 20261014, got %.4f)" % p_null,
          abs(p_null - 0.6820) < 5e-4)
    # (b) split-test shift parameter R (Sec. 5.2)
    r30 = cmb.split_test(30.0)
    _, _, _, _, vec30 = cmb.fit_bimetric(30.0)
    check("split test 30 deg: R = %.5f/%.5f/%.5f (== 1.94441/1.93500/1.78743)"
          % (vec30[0], r30["zero_high_R"], r30["zero_low_R"]),
          abs(vec30[0] - 1.9444060163205277) < 1e-4
          and abs(r30["zero_high_R"] - 1.9349999274549476) < 1e-4
          and abs(r30["zero_low_R"] - 1.7874254728192567) < 1e-4)
    r10 = cmb.split_test(10.0)
    _, _, _, _, vec10 = cmb.fit_bimetric(10.0)
    check("split test 10 deg: R = %.5f/%.5f/%.5f (== 1.76330/1.76210/1.74789)"
          % (vec10[0], r10["zero_high_R"], r10["zero_low_R"]),
          abs(vec10[0] - 1.7632984607673479) < 1e-4
          and abs(r10["zero_high_R"] - 1.7621042790222547) < 1e-4
          and abs(r10["zero_low_R"] - 1.7478903790407319) < 1e-4)
    # (c) Table A4 D_M accumulation (6 rows)
    from scipy.integrate import quad as _quad
    _m4 = Model(theta_deg=6.1, mFP_H0=25.0, alpha=0.0, beta=3.0, Omega_L=0.699,
                Omega_r=9.08e-5)
    _om4 = _m4.solve_Om0()
    _DM = _quad(lambda z: 1.0/np.sqrt(_m4.E2_bim(z, _om4)), 0, 1090.0, epsrel=1e-12,
                limit=300)[0]
    _a4ok = True
    for _zc, _e in ((1, 24.32), (5, 57.19), (10, 69.42), (12, 72.19), (50, 87.85),
                    (500, 98.51)):
        _d = _quad(lambda z: 1.0/np.sqrt(_m4.E2_bim(z, _om4)), 0, float(_zc),
                   epsrel=1e-12, limit=300)[0]
        if abs(100.0*_d/_DM - _e) > 0.01:
            _a4ok = False
    check("Table A4 D_M accumulation rows (6/6)", _a4ok)
    # (d) lambda = 1 anchor profile minimum (Sec. 7.4.2(i))
    import os as _os
    _anc_p1 = json.load(open(_os.path.join(_os.path.dirname(__file__), "..", "results",
                                           "anchor_scan.json")))
    _p1 = np.array(_anc_p1["profile_1D_lambda1"])
    _i1 = _p1[:, 4].argmin()
    check("lambda=1 anchor profile: min %.4f at %.3f (== 6.80 at 0.680)"
          % (_p1[_i1, 4], _p1[_i1, 0]),
          abs(_p1[_i1, 4] - 6.80) < 0.02 and abs(_p1[_i1, 0] - 0.680) < 0.004)
    # (e) du*/deps by implicit differentiation (independent of finite differences)
    _t6 = np.tan(np.radians(6.1)) ** 2
    _u0 = ss["ustar"]
    _hh = 1e-6
    _dXdu = (selfsim.X_of_u(_t6, 0.0, 3.0, _u0+_hh, 0.0)
             - selfsim.X_of_u(_t6, 0.0, 3.0, _u0-_hh, 0.0)) / (2*_hh)
    _dXde = (selfsim.X_of_u(_t6, 0.0, 3.0, _u0, _hh)
             - selfsim.X_of_u(_t6, 0.0, 3.0, _u0, -_hh)) / (2*_hh)
    _du_impl = -_dXde / _dXdu
    check("du*/deps via implicit differentiation = %.6f (== -0.626809)" % _du_impl,
          abs(_du_impl - (-0.626809)) < 1e-5)
    # (f) exact identities P(1) = 3 t u, P'(1) = 3[1 + t(1 - 2 eps + u)]
    _y_t = 0.37
    _pp = selfsim.P_poly(selfsim.reduced_coeffs(_t6, 0.0, 3.0, 7.7, 0.31), _y_t)
    _id_ok = abs(selfsim.P_poly(selfsim.reduced_coeffs(_t6, 0.0, 3.0, 7.7, 0.31), 1.0)
                 - 3*_t6*7.7) < 1e-12
    _c1f = selfsim.reduced_coeffs(_t6, 0.0, 3.0, 7.7, 0.31)
    _dP1 = 4*_c1f[0] + 3*_c1f[1] + 2*_c1f[2] + _c1f[3]
    _id_ok = _id_ok and abs(_dP1 - 3*(1 + _t6*(1 - 2*0.31 + 7.7))) < 1e-12
    check("exact identities P(1) = 3 t u and P'(1) = 3[1+t(1-2eps+u)]", _id_ok)

    # Mistral-review gates: convention documentation made quantitative
    # (a) Pantheon+ event counts behind the N quoted in Sec. 4
    N_all = len(sne.DF["zHD"])
    N_zHD = int((sne.DF["zHD"] > 0.01).sum())
    N_zCMB = int((sne.DF["zCMB"] > 0.01).sum())
    N_zHD23 = int((sne.DF["zHD"] > 0.023).sum())
    N_zCMB23 = int((sne.DF["zCMB"] > 0.023).sum())
    check("Pantheon+ counts: 1701 total, N(zHD>0.01)=%d (==1590), N(zCMB>0.01)=%d (==1588),"
          " N(zHD>0.023)=%d (==1371), N(zCMB>0.023)=%d (==1367)"
          % (N_zHD, N_zCMB, N_zHD23, N_zCMB23),
          N_all == 1701 and N_zHD == 1590 and N_zCMB == 1588 and N_zHD23 == 1371 and N_zCMB23 == 1367)
    # (b) CMB closure convention: neutrinos excluded from Omega_m, added as
    # omega_nu/h^2; radiation Omega_r(h) = omega_r/h^2
    m_cmb = Model(theta_deg=6.1, Omega_L=0.699, Omega_r=9.08e-5)
    om_nu_687 = cmb.solve_Om0_nu(m_cmb, 0.687)
    om_nu_68701 = cmb.solve_Om0_nu(m_cmb, 0.68701)
    check("CMB closure (neutrinos excluded, +omega_nu/h^2): Om(0.687)=%.5f (==0.30303),"
          " Om(0.68701)=%.5f (Table A5)" % (om_nu_687, om_nu_68701),
          abs(om_nu_687 - 0.30302830469870806) < 1e-5
          and abs(om_nu_68701 - 0.30302834595588546) < 1e-5)
    wr_h2 = cmb.OMEGA_GAMMA_H2 * (1.0 + 0.2271 * 3.046)   # omega_r h^2, N_eff included
    r_above = wr_h2 / 0.678 ** 2                          # the SNe/scalar fixed fraction
    r_below = cmb.OMEGA_GAMMA_H2 / 0.678 ** 2             # photons only (below the join)
    check("Omega_r convention: omega_r/h^2 at h=0.678 = %.5f (==9.08e-5, the fixed SNe/scalar"
          " fraction); photons-only = %.5f (CMB below the join)" % (r_above, r_below),
          abs(r_above - 9.08e-5) < 1e-6 and abs(r_below - 5.37e-5) < 1e-6)
    # (c) A5 theta-scan row at 6.1 deg, live
    chi2_ref6, p_ref6, scan6 = cmb.theta_scan()
    a5 = scan6[6.1]
    check("A5 6.1 deg row, live: chi2=%.5f (==0.00017), h=%.5f (==0.68688), Om=%.5f (==0.30303)"
          % (a5["chi2"], a5["h"], a5["Om"]),
          abs(a5["chi2"] - 0.0001702627819509077) < 2e-5
          and abs(a5["h"] - 0.6868767439372414) < 1e-5
          and abs(a5["Om"] - 0.303027796030993) < 1e-5)

    # D-08: tied-rd comparators are not live in the current article. Their old
    # artifacts must not be copied into this candidate or required by the gates.
    check("D-08 joint_tie deferred: no stale artifact required",
          not _os2.path.exists(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                               "joint_tie.json")))
    check("D-08 joint_tie_l1 deferred: no stale artifact required",
          not _os2.path.exists(_os2.path.join(_os2.path.dirname(__file__), "..", "results",
                                               "joint_tie_l1.json")))

    # lambda*(theta): at 10 deg the R(lambda, theta) sweep has an interior
    # minimum (compressed CMB likelihood saturates for theta >= 10 deg)
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "probes"))
    import cmb as _cmb
    from scipy.optimize import minimize as _min
    def _fit10(lam):
        mm = _cmb.Model(theta_deg=10.0, Omega_L=0.699)
        def _obj(p):
            h, wb = p
            if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
                return 1e10
            return _cmb.chi2(_cmb.vec_bimetric_m(mm, h, wb, lam)[0])
        rr = _min(_obj, [0.689, 0.022371], method="Nelder-Mead",
                  options={"xatol": 1e-9, "fatol": 1e-11, "maxiter": 3000})
        return float(rr.fun)
    c10_0 = _fit10(0.0)
    c10_5 = _fit10(0.5)
    check("lambda sweep: chi2(10 deg, lam=0) = %.4f (> 1) vs chi2(lam=0.5) = %.4f (< 0.5):"
          " interior minimum exists" % (c10_0, c10_5), c10_0 > 1.0 and c10_5 < 0.5)

    # equilibrium anchor: k_eq = aH(z_eq), z_eq = Om/Om_r - 1
    z_eq = omdef / 9.08e-5 - 1.0
    k_eq = (100.0 * 0.7 / 299792.458) * np.sqrt(
        scalarsec.E2_horizon(mdef, omdef, z_eq)) / (1.0 + z_eq)
    check("k-map equilibrium anchor z_eq == 3350 (got %.1f), k_eq == 0.0106 (got %.5f)"
          % (z_eq, k_eq), abs(z_eq - 3350.9) < 1.0 and abs(k_eq - 0.010592) < 1e-4)

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
