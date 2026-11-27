"""m_FP scan of the scalar-sector diagnostic (Sec. 6.2 of the paper).

Reproduces results/scalar_mfp_scan.json: for each m_FP/H0 in the grid,
the self-consistent closure Omega_m is used, the first X(z) > 0 -> X < 0
crossing z* is located on a wide geometric grid, I1 is integrated to
z = 1090, and the unscreened fraction is evaluated between z* and z_V
(defined by H(z_V) = m_FP). Deterministic (fixed grids, one core).
"""
import json
import pathlib
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from core import Model
import scalarsec

MFP_GRID = [1.5, 2, 3, 6, 10, 15, 20, 25, 40, 60, 100, 200]
Z_MAX = 1090.0
X_GRID = np.geomspace(0.01, 1e4, 200)


def find_zstar(m, Om):
    Xg = np.array([scalarsec.X_of_z(m, z, Om) for z in X_GRID])
    for i in range(len(X_GRID) - 1):
        if Xg[i] > 0 and Xg[i + 1] < 0:
            return brentq(lambda z: scalarsec.X_of_z(m, z, Om),
                          X_GRID[i], X_GRID[i + 1])
    return None


def main():
    rows = []
    for mf in MFP_GRID:
        m = Model(theta_deg=6.1, mFP_H0=float(mf), alpha=0.0, beta=3.0,
                  Omega_L=0.699, Omega_r=9.08e-5)
        Om = m.solve_Om0()
        zs = find_zstar(m, Om)
        if zs is None:
            rows.append(dict(mFP_H0=mf, zstar=None))
            continue
        I1, _ = scalarsec.growth_integral(m, Om, zs)
        zv = brentq(lambda z: np.sqrt(m.E2_bim(z, Om)) - mf, 0.05, 3000)
        fu = lambda z: np.sqrt(-scalarsec.X_of_z(m, z, Om)) / np.sqrt(m.E2_bim(z, Om))
        Iu, _ = quad(fu, zs, min(zv, Z_MAX), epsabs=1e-10, epsrel=1e-10, limit=300)
        rows.append(dict(mFP_H0=mf, zstar=zs, zV=zv, I1=I1,
                         unscreened_pct=100.0 * Iu / I1))
    out = pathlib.Path(__file__).resolve().parents[1] / "results" / "scalar_mfp_scan.json"
    out.write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(r)

    # the k-dilemma map (Sec. 6.2, Table 7)
    km = scalarsec.k_map()
    outk = pathlib.Path(__file__).resolve().parents[1] / "results" / "scalar_k_map.json"
    outk.write_text(json.dumps(km, indent=1))
    print("k-map:", json.dumps(km))


if __name__ == "__main__":
    main()
