"""Anchor profiling: Delta chi^2(Omega_Lambda) per probe and jointly,
separating properties of the benchmark point from properties of the model.

For each (theta, Omega_Lambda): the bimetric nuisance parameters of each probe
are profiled (BAO: f; SNe: M with the LCDM Omega_m refit; CMB: h, omega_b),
and the quartic closure supplies Omega_m(theta, Omega_Lambda). Delta chi^2 is
taken against each probe's best-fit LambdaCDM (constant). Joint = BAO(full) +
SNe(full) + CMB. Profiling, not marginalisation; descriptive arithmetic.
"""
import numpy as np
import json, pathlib, sys
from scipy.optimize import brentq, minimize, minimize_scalar

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))
from core import Model, OMEGA_R_PHOTONS
import bao, sne, cmb

OLs_1D = [0.60, 0.62, 0.64, 0.66, 0.68, 0.699, 0.72, 0.74, 0.76, 0.78]
thetas_2D = [4.0, 6.1, 8.0, 12.2, 18.3, 25.0, 30.0]
OLs_2D = [0.60, 0.64, 0.68, 0.699, 0.72, 0.76, 0.78]


def bim_point(theta_deg, OL):
    m = Model(theta_deg=theta_deg, Omega_L=OL)
    return m, m.solve_Om0()


def bao_chi2(theta_deg, OL):
    # BAO convention: Omega_r = 5.38e-5 (photons only, DESI DR2 fiducial),
    # matching src/probes/bao.py and Sec. 3 of the paper. (The pre-8 anchor
    # scan used the SNe/CMB radiation density here, which shifted the BAO
    # column by ~0.014 in Delta chi^2 at the benchmark.)
    try:
        m = Model(theta_deg=theta_deg, Omega_L=OL, Omega_r=OMEGA_R_PHOTONS)
        Om = m.solve_Om0()
    except ValueError:
        return None
    E2g = np.array([m.E2_bim(z, Om) for z in bao.ZG])
    return bao.fit_f(bao.model_from_E2(E2g))[1]


def sne_chi2(theta_deg, OL):
    try:
        m, Om = bim_point(theta_deg, OL)
    except ValueError:
        return None
    E2g = np.array([m.E2_bim(z, Om) for z in sne.Z_FINE])
    from scipy.integrate import cumulative_trapezoid
    from scipy.interpolate import interp1d
    ic = cumulative_trapezoid(1.0 / np.sqrt(E2g), sne.Z_FINE, initial=0)
    DMi = interp1d(sne.Z_FINE, ic, kind="cubic", fill_value="extrapolate")
    zcut = sne.DF["zHD"] > 0.01
    zc, zh = sne.DF["zHD"][zcut], sne.DF["zHEL"][zcut]
    DL = (1.0 + zh) * sne.C_H0 * DMi(zc)
    mu = 5.0 * np.log10(DL) + 25.0
    cov = sne.DF["COV_FULL"][np.ix_(zcut, zcut)]
    Cinv = np.linalg.inv(cov)
    r = sne.DF["m_b_corr"][zcut] - mu
    Cr = Cinv @ r
    return float(r @ Cr - np.sum(Cr) ** 2 / np.sum(Cinv))


OLs_1D = [0.60, 0.62, 0.64, 0.66, 0.68, 0.699, 0.72, 0.74, 0.76, 0.78]
OLs_FINE = [round(0.66 + 0.004 * i, 3) for i in range(21)]        # 0.66-0.74
THETAS_FINE = [4.0, 4.5, 5.0, 5.5, 6.0, 6.1, 6.5, 7.0, 7.5, 8.0]
thetas_2D = [4.0, 6.1, 8.0, 12.2, 18.3, 25.0, 30.0]
OLs_2D = [0.60, 0.64, 0.68, 0.699, 0.72, 0.76, 0.78]
# extension for the interior-minimum check at large theta (Table 4)
THETAS_EXT = [18.3, 25.0, 30.0]
OLs_EXT = [0.80, 0.82, 0.86, 0.90]


def probe_row(theta_deg, OL, r_map="quartic"):
    """[OL, BAO, SNe, CMB, joint] ABSOLUTE chi^2 of each probe and their sum."""
    cb = bao_chi2(theta_deg, OL)
    cs = sne_chi2(theta_deg, OL)
    cc = cmb.fit_bimetric(theta_deg, OL, r_map=r_map)[0]
    return [OL, cb, cs, cc, None if None in (cb, cs, cc) else cb + cs + cc]


def _save(out, path):
    json.dump(out, open(path, "w"), indent=2)


def main():
    bao.load_desi()
    sne.load_pantheon()
    chi2_ref, p_ref = cmb.fit_lcdm_reference()
    bao_lcdm = bao.fit_lcdm()[2]
    sne_lcdm = sne.run(0.01)["full"]["chi2_min"]
    refs = dict(bao=bao_lcdm, sne=sne_lcdm, cmb=chi2_ref)
    print("refs:", refs, flush=True)

    out = {"refs": refs}
    path = ROOT / "results/anchor_scan.json"

    def row1d(OL, r_map="quartic"):
        cb = bao_chi2(6.1, OL)
        cs = sne_chi2(6.1, OL)
        cc = cmb.fit_bimetric(6.1, OL, r_map=r_map)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        joint = None if None in vals else sum(vals)
        return [OL] + vals + [joint]

    row = []
    for OL in OLs_1D:
        r = row1d(OL)
        row.append(r)
        print("1D OL=%.3f BAO=%.3f SNe=%.3f CMB=%.3f joint=%.3f" %
              (r[0], r[1], r[2], r[3], r[4]), flush=True)
    out["profile_1D"] = row

    # lambda = 1 (effective-density map) diagnostic at theta = 6.1 deg:
    # BAO and SNe do not depend on the R mapping, so only CMB is recomputed
    row_l1 = []
    for OL in OLs_1D:
        cb = bao_chi2(6.1, OL)
        cs = sne_chi2(6.1, OL)
        cc = cmb.fit_bimetric(6.1, OL, r_map=1.0)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        row_l1.append([OL] + vals + [sum(vals)])
        print("1D-l1 OL=%.3f CMB=%.3f joint=%.3f" % (OL, vals[2], sum(vals)),
              flush=True)
    out["profile_1D_lambda1"] = row_l1

    # fine 1D profile (0.004 spacing, 0.66-0.74) at theta = 6.1 deg
    fine = []
    for OL in OLs_FINE:
        fine.append(row1d(OL))
        print("fine OL=%.3f joint=%.3f" % (fine[-1][0], fine[-1][4]), flush=True)
    out["profile_1D_fine"] = fine

    # fine theta scan at OL = 0.699
    tfine = []
    for th in THETAS_FINE:
        cb = bao_chi2(th, 0.699)
        cs = sne_chi2(th, 0.699)
        cc = cmb.fit_bimetric(th, 0.699)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        tfine.append([th] + vals + [sum(vals)])
        print("theta_fine th=%.1f joint=%.3f" % (th, sum(vals)), flush=True)
    out["theta_fine"] = tfine

    # 2D grid
    grid = []
    for th in thetas_2D:
        for OL in OLs_2D:
            cb = bao_chi2(th, OL)
            cs = sne_chi2(th, OL)
            cc = cmb.fit_bimetric(th, OL)[0]
            vals = [x - r for x, r in zip((cb, cs, cc),
                                          (refs["bao"], refs["sne"], refs["cmb"]))]
            joint = None if None in vals else sum(vals)
            grid.append([th, OL] + vals + [joint])
            print("2D th=%.1f OL=%.3f joint=%s" % (th, OL,
                  "%.3f" % joint if joint is not None else "None"), flush=True)
    out["grid_2D"] = grid

    # extension to Omega_Lambda = 0.90 at theta = 18.3, 25, 30 deg
    # (interior-minimum verification for Table 4)
    ext = []
    for th in THETAS_EXT:
        for OL in OLs_EXT:
            cb = bao_chi2(th, OL)
            cs = sne_chi2(th, OL)
            cc = cmb.fit_bimetric(th, OL)[0]
            vals = [x - r for x, r in zip((cb, cs, cc),
                                          (refs["bao"], refs["sne"], refs["cmb"]))]
            joint = None if None in vals else sum(vals)
            ext.append([th, OL] + vals + [joint])
            print("ext th=%.1f OL=%.3f joint=%s" % (th, OL,
                  "%.3f" % joint if joint is not None else "None"), flush=True)
    out["grid_2D_ext"] = ext

    # joint LambdaCDM reference with one shared Omega_m
    Om_s, joint_s, cost = lcdm_shared_cost()
    out["lcdm_shared"] = dict(Om=Om_s, joint=joint_s, cost=cost,
                              chi2_sum_min=joint_s - cost)   # sum of the per-probe LCDM minima
    print("lcdm_shared: Om=%.4f joint=%.3f cost=%.4f" % (Om_s, joint_s, cost),
          flush=True)

    # continuous Omega_Lambda minimisation at each theta (Table 7, Table 8)
    rows = [[th, round(ol, 4), round(j, 2)] for th, ol, j in theta_profile(refs=refs)]
    for r in rows:
        print("theta_profile th=%.1f OL*=%.4f joint=%.2f" % tuple(r), flush=True)
    out["theta_profile_continuous"] = dict(
        note="Omega_Lambda minimised continuously at each theta", rows=rows)

    # continuous profile at theta = 6.1 with the lambda = 1 R-slot filling
    l1 = lambda1_profile()
    print("lambda1 profile:", l1, flush=True)
    out["profile_lambda1_continuous"] = dict(
        OL_star=round(l1["OL_star"], 4), joint_min=round(l1["joint_min"], 2),
        note="continuous Omega_Lambda profile at theta=6.1, lambda=1")

    _save(out, path)
    print("ANCHOR SCAN DONE", flush=True)



# ---------------------------------------------------------------------------
# Addendum (pre-release 2 review): joint reference, refinement, map bracket
# ---------------------------------------------------------------------------

def lcdm_shared_cost(oms=None):
    """Joint LambdaCDM with a SINGLE shared Omega_m across the three probes,
    nuisance-profiled per probe (f; M; h, omega_b with omega_m = Omega_m h^2).
    Returns (Omega_m at minimum, min joint chi^2, cost vs per-probe refs)."""
    from scipy.optimize import minimize as _min
    if oms is None:
        oms = np.linspace(0.28, 0.36, 17)
    sne.load_pantheon()
    bao.load_desi()
    zcut = sne.DF["zHD"] > 0.01   # primary zHD convention (Sec. 4.1)
    cov = sne.DF["COV_FULL"][np.ix_(zcut, zcut)]
    Cinv = np.linalg.inv(cov)
    S = float(np.sum(Cinv))

    def sne_c2(Om):
        mu = sne.lcdm_mu(Om, zcut)
        r = sne.DF["m_b_corr"][zcut] - mu
        Cr = Cinv @ r
        return float(r @ Cr - np.sum(Cr) ** 2 / S)

    def cmb_c2(Om):
        def obj(p):
            h, wb = p
            if not (0.5 < h < 0.9 and 0.015 < wb < 0.025):
                return 1e10
            # Omega_cb = Om - omega_nu/h^2: the shared Om is full matter,
            # the CMB vector takes neutrinos separately (Sec. 7.4.2(i))
            return cmb.chi2(cmb.vec_lcdm(Om * h ** 2 - cmb.OMEGA_NU, h, wb))
        best, x0 = 1e10, None
        for h in np.linspace(0.60, 0.80, 13):
            for wb in np.linspace(0.020, 0.023, 7):
                c = obj([h, wb])
                if c < best:
                    best, x0 = c, [h, wb]
        return float(_min(obj, x0, method="Nelder-Mead",
                          options={"xatol": 1e-10, "fatol": 1e-12,
                                   "maxiter": 1500}).fun)

    rows = []
    for Om in oms:
        cb = bao.fit_f(bao.model_from_E2(
            bao.E2_lcdm(bao.ZG, Om, bao.OMEGA_R_PHOTONS)))[1]
        rows.append([float(Om), cb, sne_c2(Om), cmb_c2(Om)])
    joint = [r[1] + r[2] + r[3] for r in rows]
    imin = int(np.argmin(joint))
    refs = (bao.fit_lcdm()[2], sne.run(0.01)["full"]["chi2_min"],
            cmb.fit_lcdm_reference()[0])
    cost = joint[imin] - sum(refs)
    return rows[imin][0], joint[imin], cost


def lambda1_profile(theta_deg=6.1, lo=0.66, hi=0.74):
    """Continuous Omega_Lambda profile of the joint chi^2 at fixed theta with
    the lambda = 1 (effective-density) R-slot filling (pre-release 6;
    stored as profile_lambda1_continuous)."""
    from scipy.optimize import minimize_scalar
    import probes.cmb as _cmb
    def joint_l1(OL):
        cb = bao_chi2(theta_deg, OL) - 10.279
        cs = sne_chi2(theta_deg, OL) - 1402.921
        cc = _cmb.fit_bimetric(theta_deg, OL, r_map=1.0)[0]
        return cb + cs + cc
    r = minimize_scalar(joint_l1, bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-7})
    return dict(theta_deg=float(theta_deg), OL_star=float(r.x),
                joint_min=float(r.fun))


def theta_profile(thetas=(0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.1, 8.0, 10.0,
                      12.2, 18.3, 25.0, 30.0), lo=0.66, hi=0.90, refs=None):
    """Continuous Omega_Lambda profile of the joint chi^2 at each theta
    (pre-release 6: replaces the fixed-0.699 slice and the 0.02-spaced grid
    for the per-theta minima; stored as theta_profile_continuous)."""
    from scipy.optimize import minimize_scalar
    # joint Delta chi^2 = sum of absolute chi^2 minus the three per-probe LCDM minima
    base = sum(refs.values()) if refs is not None else 0.0
    out = []
    for th in thetas:
        r = minimize_scalar(lambda OL: probe_row(th, OL)[4] - base, bounds=(lo, hi),
                            method="bounded", options={"xatol": 1e-7})
        out.append([th, float(r.x), float(r.fun)])
    return out


if __name__ == "__main__":
    main()
