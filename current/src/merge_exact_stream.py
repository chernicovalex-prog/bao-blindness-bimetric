#!/usr/bin/env python3
"""Concatenate exact-stream realisations and compute statistics once centrally."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EVD = ROOT / "evidence"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(value, f, indent=2)
    os.replace(tmp, path)


def statistics(outputs, arrays):
    first = outputs[0]
    cov1 = arrays["cov1"]
    cov384 = arrays["cov384"]
    dchi = arrays["dchi_true"]
    argmin = arrays["prof_argmin"]
    return {"N": len(dchi), "seed": first["seed"], "seeds": [first["seed"]],
            "theta": first["theta"], "OL_true": first["OL_true"],
            "grid_mode": first["grid_mode"], "statistic_mode": first["statistic_mode"],
            "covariance_scheme": first["covariance_scheme"],
            "cmb_linbias_real": first["cmb_linbias_real"],
            "OL_grid": first["OL_grid"], "grid_spacing": first["grid_spacing"],
            "coverage_dchi1": round(float(np.mean(cov1)), 4),
            "coverage_dchi384": round(float(np.mean(cov384)), 4),
            "dchi_true_quantiles": {str(q): round(float(np.quantile(dchi, q)), 4)
                                    for q in (0.5, 0.68, 0.9, 0.95)},
            "argmin_mean": round(float(np.mean(argmin)), 5),
            "argmin_quantiles": {str(q): round(float(np.quantile(argmin, q)), 4)
                                 for q in (0.16, 0.5, 0.84)}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["control", "full"], required=True)
    args = parser.parse_args()
    if args.mode == "control":
        tags = ["exactctrl0b", "exactctrl1b"]
        ends = [10, 20]
        merged_path = ROOT / "results/exact_stream_control_merged.json"
        report_path = EVD / "G07_EQUIVALENCE_CONTROL.json"
    else:
        tags = ["exactF0", "exactF1", "exactF2", "exactF3"]
        ends = [500, 1000, 1500, 2000]
        merged_path = ROOT / "results/profile_calibration.json"
        report_path = EVD / "G07_EXACT_STREAM_MERGE.json"

    checkpoints, outputs, provenance = [], [], []
    for tag, end in zip(tags, ends):
        ckpt = ROOT / f"results/profile_calibration_checkpoint_{tag}.json"
        final = ROOT / f"results/profile_calibration_{tag}.json"
        prov = ROOT / f"results/exact_stream_{tag}_provenance.json"
        for path in (ckpt, final, prov):
            if not path.is_file():
                print(json.dumps({"status": "SEGMENTS_NOT_COMPLETE", "missing": str(path)}))
                return 2
        saved = json.loads(ckpt.read_text())
        if saved["done"] != end:
            raise RuntimeError(f"{tag}: done={saved['done']} != {end}")
        checkpoints.append(saved)
        outputs.append(json.loads(final.read_text()))
        provenance.append(json.loads(prov.read_text()))

    signatures = [p["code_data_signature"] for p in provenance]
    if any(sig != signatures[0] for sig in signatures[1:]):
        raise RuntimeError("segment code/data signatures differ")
    arrays = {}
    for key in ("prof_argmin", "cov1", "cov384", "dchi_true"):
        arrays[key] = [value for saved in checkpoints for value in saved[key]]
    expected = ends[-1]
    if len(arrays["dchi_true"]) != expected:
        raise RuntimeError("merged realisation count mismatch")
    stats = statistics(outputs, arrays)
    stats.update({"exact_stream": True,
                  "segments": [{"tag": tag, "start": p["start"], "end": p["end"]}
                               for p, tag in zip(provenance, tags)],
                  "merge_method": "concatenate per-realisation arrays in stream order; statistics computed once centrally",
                  "code_data_signature": signatures[0]})
    save_json(merged_path, stats)

    report = {"status": "MERGE_COMPLETE", "mode": args.mode, "records": expected,
              "merged_file": str(merged_path), "merged_sha256": sha(merged_path),
              "segment_checkpoint_sha256": {tag: sha(ROOT / f"results/profile_calibration_checkpoint_{tag}.json")
                                             for tag in tags}}
    if args.mode == "control":
        single_ckpt_path = ROOT / "results/profile_calibration_checkpoint_exactsingle20.json"
        single_out_path = ROOT / "results/profile_calibration_exactsingle20.json"
        single_ckpt = json.loads(single_ckpt_path.read_text())
        single_out = json.loads(single_out_path.read_text())
        array_checks = {key: arrays[key] == single_ckpt[key] for key in arrays}
        stat_keys = ["N", "seed", "seeds", "theta", "OL_true", "grid_mode", "statistic_mode",
                     "covariance_scheme", "cmb_linbias_real", "OL_grid", "grid_spacing",
                     "coverage_dchi1", "coverage_dchi384", "dchi_true_quantiles",
                     "argmin_mean", "argmin_quantiles"]
        stat_checks = {key: stats[key] == single_out[key] for key in stat_keys}
        rng_match = checkpoints[-1]["rng_state"] == single_ckpt["rng_state"]
        passed = all(array_checks.values()) and all(stat_checks.values()) and rng_match
        report.update({"status": "EQUIVALENCE_CONTROL_PASS" if passed else "EQUIVALENCE_CONTROL_FAIL",
                       "array_checks": array_checks, "statistics_checks": stat_checks,
                       "final_rng_state_match": rng_match,
                       "single_checkpoint": str(single_ckpt_path),
                       "single_output": str(single_out_path),
                       "full_run_may_start": passed})
    else:
        report["final_rng_state"] = checkpoints[-1]["rng_state"]
    save_json(report_path, report)
    print(json.dumps(report, indent=2))
    return 0 if report["status"] in {"EQUIVALENCE_CONTROL_PASS", "MERGE_COMPLETE"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
