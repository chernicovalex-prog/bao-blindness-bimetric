"""Independent density/recombination specification gates, not article assertions.

Run only as an authorized focused check. Outputs JSON with T1–T6.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "src" / "probes"))
import cmb as C

M = C.Model(theta_deg=6.1, Omega_L=0.699)
H_LIST = (0.60, 0.70, 0.80)
Z_LIST = (0, 1, 5, 10, C.Z_NR - 2, C.Z_NR - 1, C.Z_NR, 1090)
TOL = 1e-10


def residual(coeffs, y):
    c4, c3, c2, c1, c0 = coeffs
    p = ((((c4 * y + c3) * y + c2) * y + c1) * y + c0)
    den = max(1.0, abs(c0) + abs(c1 * y) + abs(c2 * y ** 2)
              + abs(c3 * y ** 3) + abs(c4 * y ** 4))
    return abs(p) / den


def t1():
    worst, points = 0.0, []
    for h in H_LIST:
        om = C.solve_Om0_nu(M, h)
        for z in Z_LIST:
            y = C.solve_y_cmb(M, z, om, h)
            r = residual(M.quartic_coeffs(C.omega_ng(z, om, h)), y)
            worst = max(worst, r)
            points.append({"h": h, "z": z, "r_rel": r})
    return {"status": "PASS" if worst <= TOL else "FAIL", "max_r_rel": worst,
            "num_points": len(points)}


def t2():
    worst = 0.0
    for h in H_LIST:
        om = C.solve_Om0_nu(M, h)
        y0 = C.solve_y_cmb(M, 0.0, om, h)
        e20 = C.omega_ng(0.0, om, h) + C.bim_ode_y(M, y0)
        worst = max(worst, abs(e20 - 1.0))
    return {"status": "PASS" if worst <= TOL else "FAIL", "max_residual": worst}


def dm_direct(om, h, zs):
    def E2(z):
        return C.omega_ng(z, om, h) + C.bim_ode_y(M, C.solve_y_cmb(M, z, om, h))
    return quad(lambda z: 1.0 / np.sqrt(E2(z)), 0.0, zs,
                epsabs=1e-11, epsrel=1e-10, limit=300)[0]


def dm_spline(om, h, zs):
    spl = C._y_spline(M, om, h)
    def E2(z):
        dy = 1.0 - float(spl(z))
        ode = M.Omega_L - M.sin2 * M.mFP2 * dy * (1 + M.alpha * dy + M.beta / 3 * dy ** 2)
        return C.omega_ng(z, om, h) + ode
    return quad(lambda z: 1.0 / np.sqrt(E2(z)), 0.0, zs,
                epsabs=1e-11, epsrel=1e-10, limit=300)[0]


def t3():
    h, wb, om = 0.687014059239027, 0.02237204963071906, 0.3030283627026459
    zs = C.z_star_hs(wb, om * h ** 2 + C.OMEGA_NU, h)
    direct, spline = dm_direct(om, h, zs), dm_spline(om, h, zs)
    rel = abs(direct - spline) / abs(direct)
    return {"status": "PASS" if rel <= 1e-5 else "FAIL", "rel_difference": rel,
            "h": h, "omega_b": wb, "Omega_m": om, "z_rec": zs,
            "DM_direct": direct, "DM_spline": spline}


def t4():
    C._Y_CACHE.clear()
    om = C.solve_Om0_nu(M, H_LIST[0])
    spl60 = C._y_spline(M, om, H_LIST[0])
    spl80 = C._y_spline(M, om, H_LIST[2])
    again = C._y_spline(M, om, H_LIST[0])
    different = spl60 is not spl80
    repeat_same = again is spl60
    values_differ = any(float(spl60(z)) != float(spl80(z)) for z in (0, 1, 10, 100, 1090))
    ok = different and repeat_same and values_differ
    return {"status": "PASS" if ok else "FAIL",
            "distinct_objects": different,
            "same_h_same_object": repeat_same,
            "different_h_values_differ": values_differ,
            "cache_size": len(C._Y_CACHE)}


def t5():
    core = ROOT / "src" / "core.py"
    digest = hashlib.sha256(core.read_bytes()).hexdigest()
    # Reviewed docstring-only correction: H&M finite branch is Sec. 5.2.
    expected = "6636ea2520b5eaa5f9f6718cd7594c4bd77f918e9959a7e79c8a0876ad7735ae"
    return {"status": "PASS" if digest == expected else "FAIL",
            "core_py_matches_reviewed": digest == expected,
            "sha256": digest}


def t6():
    rows, exceed = [], 0
    for h in H_LIST:
        om = C.solve_Om0_nu(M, h)
        for z in Z_LIST:
            y = M.solve_y(z, om)
            r = residual(M.quartic_coeffs(C.omega_ng(z, om, h)), y)
            rows.append({"h": h, "z": z, "r_rel": r})
            exceed += int(r > TOL)
    return {"status": "PASS" if exceed else "FAIL", "num_exceeding": exceed,
            "num_points": len(rows), "max_r_rel": max(r["r_rel"] for r in rows)}


def main():
    result = {
        "T1_quartic_residual": t1(),
        "T2_closure": t2(),
        "T3_direct_vs_spline": t3(),
        "T4_cache": t4(),
        "T5_non_cmb_regression": t5(),
        "T6_negative_control": t6(),
    }
    print(json.dumps(result, indent=2))
    return 0 if all(row["status"] == "PASS" for row in result.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
