"""lambda-sweep of the compressed CMB likelihood at each theta (Sec. 5.2, Table A5).
Writes results/cmb_lambda_theta_sweep.json."""
import sys, json, pathlib
SRC = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))
from scipy.optimize import minimize
import cmb

LAMBDAS = [0.0, 0.25, 0.5, 0.75, 1.0]
THETAS = [6.1, 10.0, 12.2, 15.0, 18.3, 20.0, 25.0, 30.0]

def fit_warm(theta, lam, x0):
    m = cmb.Model(theta_deg=theta, Omega_L=0.699)
    def obj(p):
        h, wb = p
        if not (0.55 < h < 0.88 and 0.015 < wb < 0.025):
            return 1e10
        return cmb.chi2(cmb.vec_bimetric_m(m, h, wb, lam)[0])
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-9, "fatol": 1e-11, "maxiter": 3000})
    return float(res.fun), list(res.x)

out = {}
for th in THETAS:
    c0, x0 = fit_warm(th, 0.0, [0.687, 0.022371])
    row = [c0]; xs = [x0]
    for lam in LAMBDAS[1:]:
        c, x = fit_warm(th, lam, xs[-1])
        row.append(c); xs.append(x)
    out[str(th)] = dict(lambdas=LAMBDAS, chi2=row, h=[x[0] for x in xs])
    print("theta=%5.1f  " % th + "  ".join("lam=%.2f:%9.4f" % (l, c) for l, c in zip(LAMBDAS, row)), flush=True)

json.dump(out, open(SRC.parent / "results" / "cmb_lambda_theta_sweep.json", "w"), indent=1)
print("DONE")
