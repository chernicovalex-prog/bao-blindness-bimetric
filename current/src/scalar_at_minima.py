#!/usr/bin/env python3
"""D-09 scalar diagnostics at fixed anchors and regenerated profile minima."""
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import scalarsec

ANCHOR = ROOT / "results/anchor_scan.json"
OUT = ROOT / "results/scalar_at_minima.json"
FIXED_THETAS = (4.0, 6.0, 6.1, 8.0)
PROFILE_THETAS = (12.2, 18.3, 25.0, 30.0)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(theta, omega_l, family):
    result = scalarsec.run_scalar(theta_deg=float(theta), mFP_H0=25.0,
                                  alpha=0.0, beta=3.0, Omega_L=float(omega_l))
    return {"theta_deg": float(theta), "Omega_L": float(omega_l), "family": family,
            "mFP_H0": 25.0, "zstar": result["zstar"], "Hzstar": result["Hzstar"],
            "I1": result["I1"], "I1_err": result["I1_err"],
            "log10Gamma_k0.05_h0.7": result["log10_gamma"][(0.05, 0.700)]}


def main():
    anchor = json.loads(ANCHOR.read_text())
    profile = {float(row[0]): float(row[1])
               for row in anchor["theta_profile_continuous"]["rows"]}
    rows = [evaluate(theta, 0.699, "fixed_anchor_OL_0.699") for theta in FIXED_THETAS]
    rows += [evaluate(theta, profile[theta], "continuous_profile_minimum")
             for theta in PROFILE_THETAS]
    result = {
        "schema": "scalar_at_minima_v1",
        "source_anchor_json_sha256": sha(ANCHOR),
        "profile_OL_rounding": "four decimals, exactly as stored in theta_profile_continuous.rows",
        "scalar_protocol": "mFP_H0=25, alpha=0, beta=3, self-consistent closure, k=0.05 Mpc^-1, h=0.700",
        "rows": rows,
    }
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(result, indent=2) + "\n")
    os.replace(tmp, OUT)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
