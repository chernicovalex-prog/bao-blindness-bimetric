"""Merge Monte-Carlo calibration chunks into the canonical JSON.

Two supported protocols:

A) Chunked: four runs, 500 realisations each,
       python3 profile_calibration.py 500 20261118 F0 fine
       python3 profile_calibration.py 500 20261119 F1 fine
       python3 profile_calibration.py 500 20261120 F2 fine
       python3 profile_calibration.py 500 20261121 F3 fine
   each writing results/profile_calibration_F{i}.json, then
       python3 merge_calibration.py
   which combines N, seeds and the coverages (count-weighted; the script
   stores quantiles, not per-realisation values, so counts are what merge).

B) Single run: python3 profile_calibration.py 2000 20261118 main fine
   writes results/profile_calibration_main.json; this script promotes it to
   results/profile_calibration.json with seeds=[20261118].

Both protocols are deterministic (fixed seeds; the spot-check stream of
cmb_linbias_check.py is independent and does not touch the calibration).
"""
import json, os, sys

R = os.path.join(os.path.dirname(__file__), "..", "results")
chunks = [os.path.join(R, f"profile_calibration_F{i}.json") for i in range(4)]
main = os.path.join(R, "profile_calibration_main.json")
outp = os.path.join(R, "profile_calibration.json")

if all(os.path.exists(p) for p in chunks):
    ds = [json.load(open(p)) for p in chunks]
    for c in ds:
        assert c["N"] == 500 and c["grid_mode"] == "fine", (c["N"], c["grid_mode"])
    d = dict(ds[0])
    d["N"] = sum(c["N"] for c in ds)
    d["seeds"] = [c["seed"] for c in ds]
    d["seed"] = ds[0]["seed"]
    d["coverage_dchi1"] = sum(c["coverage_dchi1"] * c["N"] for c in ds) / d["N"]
    d["coverage_dchi384"] = sum(c["coverage_dchi384"] * c["N"] for c in ds) / d["N"]
    d["quantiles_from"] = "chunk F0 (quantiles do not merge; coverages are count-weighted)"
elif os.path.exists(main):
    d = json.load(open(main))
    d["seeds"] = [d["seed"]]
else:
    sys.exit("no calibration outputs found (run profile_calibration.py first)")

json.dump(d, open(outp, "w"), indent=1)
print("wrote", outp, "N=", d["N"], "seeds=", d["seeds"],
      "cov1=%.4f cov3.84=%.4f" % (d["coverage_dchi1"], d["coverage_dchi384"]))
