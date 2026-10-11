#!/usr/bin/env python3
"""
Self-checks + regression test against canonical numbers.

Run:  python3 tests/run_checks.py           (fast: core + BAO + 2 CMB points)
      python3 tests/run_checks.py --full    (full 8-angle CMB scan, minutes)

Exit status 0 = all checks pass, 1 = failure. Used by CI.
"""

import sys
import json
import numpy as np

sys.path.insert(0, '.')

from bimetric import (quartic_coeffs, solve_y, E2_bimetric, solve_Om0,
                      compute_B_params, invert_B_params,
                      OMEGA_R_H2, OMEGA_NU)


def report(name, ok, detail=""):
    print(f"  [{'OK' if ok else 'FAIL'}] {name} {detail}")
    return ok


def check_quartic_identity():
    """P(1) = Omega_tot for random parameters (identity of the A13 coeffs)."""
    rng = np.random.default_rng(20250501)
    worst = 0.0
    for _ in range(300):
        th = rng.uniform(0, 0.5 * np.pi)
        s2, c2 = np.sin(th) ** 2, np.cos(th) ** 2
        al, be = rng.uniform(-1, 1), rng.uniform(-1, 4)
        m2 = 10.0 ** rng.uniform(0, 4)
        OL = rng.uniform(0, 1)
        Ot = 10.0 ** rng.uniform(-3, 6)
        a0, a1, a2, a3, a4 = quartic_coeffs(0.0, 0.0, m2, s2, c2, al, be,
                                            OL, 0.0)
        # quartic_coeffs evaluated at Ot=0: identity reads
        # a0 + (Ot + a1_off) + a2 + a3 + a4 = Ot  <=>  P(1) = Omega_tot
        S = a0 + a1 + a2 + a3 + a4 + Ot
        worst = max(worst, abs(S - Ot) / max(1.0, abs(Ot)))
    return worst < 1e-10, f"max rel dev = {worst:.2e}"


def check_brent_vs_roots():
    th = np.radians(15.0)
    s2, c2 = np.sin(th) ** 2, np.cos(th) ** 2
    Or, wnu = OMEGA_R_H2 / 0.68 ** 2, OMEGA_NU / 0.68 ** 2
    Om0 = solve_Om0(625.0, s2, c2, 0.0, 3.0, 0.699, Or, wnu)
    worst, bad = 0.0, 0
    for z in (0.0, 0.5, 3.0, 30.0, 300.0, 1000.0):
        yb = solve_y(z, Om0, 625.0, s2, c2, 0.0, 3.0, 0.699, Or)
        Ot = Om0 * (1 + z) ** 3 + Or * (1 + z) ** 4
        a0, a1, a2, a3, a4 = quartic_coeffs(z, Om0, 625.0, s2, c2, 0.0, 3.0,
                                            0.699, Or)
        roots = np.roots([a4, a3, a2, a1, a0])
        inside = [r.real for r in roots
                  if abs(r.imag) < 1e-8 * max(1, abs(r)) and 1e-12 < r.real < 1]
        if len(inside) != 1:
            bad += 1
        else:
            worst = max(worst, abs(inside[0] - yb))
    return (bad == 0 and worst < 1e-7), f"|brent-roots| = {worst:.2e}"


def check_baseline():
    # Canonical baseline convention (archive section 1): Or = 5.38e-5, no
    # massive-nu term in the closure -> Omega_m = 0.30439, y(0) = 0.999512.
    th = np.radians(6.1)
    s2, c2 = np.sin(th) ** 2, np.cos(th) ** 2
    Om0 = solve_Om0(625.0, s2, c2, 0.0, 3.0, 0.699, 5.38e-5, 0.0)
    y0 = solve_y(0.0, Om0, 625.0, s2, c2, 0.0, 3.0, 0.699, 5.38e-5)
    E20 = E2_bimetric(0.0, Om0, 625.0, s2, c2, 0.0, 3.0, 0.699, 5.38e-5, 0.0)
    B = compute_B_params(6.1, 25.0, 0.0, 3.0, 0.699)
    Bi = invert_B_params(B, 6.1)
    ok = (abs(y0 - 0.999512) < 1e-6 and abs(Om0 - 0.30439) < 1e-4
          and abs(E20 - 1.0) < 1e-12 and abs(B['B1'] - 28.230) < 1e-3
          and abs(Bi['mFP2'] - 625.0) < 1e-9)
    return ok, (f"y0={y0:.6f} Om={Om0:.6f} E2(0)-1={E20-1:.1e} "
                f"B1={B['B1']:.3f}")


def check_bao():
    from probes import bao
    r = bao.run('full')
    return abs(r['dchi2'] - 0.178) < 0.005, f"dchi2 = {r['dchi2']:.4f} (canon 0.178)"


def check_cmb_point(theta, canon, tol):
    from probes.cmb_theta_scan import run_point
    r = run_point(theta)
    return abs(r['chi2'] - canon) < tol, f"chi2({theta}) = {r['chi2']:.3f} (canon {canon})"


def main():
    full = '--full' in sys.argv
    canon = json.load(open('canonical_numbers.json'))['cmb_theta_scan']
    ok = True
    print("Self-checks (bimetric core, archive v10.6 physics)")
    ok &= report("quartic identity P(1)=Omega_tot", *check_quartic_identity())
    ok &= report("brent vs numpy.roots (15 deg)", *check_brent_vs_roots())
    ok &= report("baseline y0/Om0/E2(0)/B/round-trip", *check_baseline())
    ok &= report("BAO full-cov dchi2", *check_bao())
    ok &= report("CMB LCDM sanity", *_lcdm())
    ok &= report("CMB 6.1 deg", *check_cmb_point(6.1, canon['6.1']['value'],
                                                 canon['6.1']['tol']))
    ok &= report("CMB 30 deg", *check_cmb_point(30.0, canon['30.0']['value'],
                                                canon['30.0']['tol']))
    if full:
        for t in ['10.0', '12.2', '15.0', '18.3', '20.0', '25.0']:
            e = canon[t]
            ok &= report(f"CMB {t} deg",
                         *check_cmb_point(float(t), e['value'], e['tol']))
    print("=" * 60)
    print("ALL OK" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)


def _lcdm():
    from probes.cmb_theta_scan import lcdm_sanity
    r = lcdm_sanity()
    return abs(r['chi2'] - 1.41) < 0.05, f"chi2(LCDM) = {r['chi2']:.3f} (canon 1.41)"


if __name__ == '__main__':
    main()
