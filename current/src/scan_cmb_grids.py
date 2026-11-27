"""CMB scans over the potential parameters (Table 2, Sec. 5.2 of the paper).

Regenerates results/cmb_grids.json:
  grid_ab  : (alpha, beta) grid at theta = 6.1 deg, m_FP/H0 = 25 (19 non-degenerate
             points; the degenerate locus 1 + 2 alpha + beta = 0 carries no finite
             branch and is excluded)
  scan_mfp : m_FP/H0 scan at the benchmark (alpha, beta)
Each entry is [chi2, h, wb, Omega_m(CMB closure)]; chi2 equals Delta chi2 against
the saturated free LambdaCDM reference (chi2_ref = 0).
"""
import sys, json, pathlib
import numpy as np
from scipy.optimize import minimize

SRC = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))
import cmb
from core import Model

ALPHAS = (-1.0, -0.5, 0.0, 0.5, 1.0)
BETAS = (1, 2, 3, 4)
MFPS = (5, 15, 25, 45, 100)


def fit_point(m, x0=(0.687, 0.022371)):
    def obj(p):
        h, wb = p
        if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
            return 1e10
        return cmb.chi2(cmb.vec_bimetric_m(m, h, wb, "quartic")[0])
    res = minimize(obj, list(x0), method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 4000})
    h, wb = res.x
    vec, Om = cmb.vec_bimetric_m(m, h, wb, "quartic")
    return [float(cmb.chi2(vec)), float(h), float(wb), float(Om)]


def main():
    out = {"grid_ab": {}, "scan_mfp": {}}
    for a in ALPHAS:
        for b in BETAS:
            if abs(1.0 + 2.0 * a + b) < 1e-12:
                continue          # degenerate locus B1 = 0: no finite branch
            m = Model(theta_deg=6.1, mFP_H0=25.0, alpha=a, beta=float(b), Omega_L=0.699)
            out["grid_ab"]["(%s, %d)" % (a, b)] = fit_point(m)
            print("alpha=%5.2f beta=%d  chi2=%.6f" % (a, b, out["grid_ab"]["(%s, %d)" % (a, b)][0]),
                  flush=True)
    for mf in MFPS:
        m = Model(theta_deg=6.1, mFP_H0=float(mf), alpha=0.0, beta=3.0, Omega_L=0.699)
        out["scan_mfp"][str(mf)] = fit_point(m)
        print("m_FP/H0=%4d  chi2=%.6f" % (mf, out["scan_mfp"][str(mf)][0]), flush=True)
    json.dump(out, open(SRC.parent / "results" / "cmb_grids.json", "w"), indent=1)
    print("stored results/cmb_grids.json")


if __name__ == "__main__":
    main()
