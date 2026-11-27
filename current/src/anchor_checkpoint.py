#!/usr/bin/env python3
"""Resumable orchestration of anchor_scan.py; scientific functions are unchanged.

Writes one row at a time to results/anchor_scan_checkpoint.json and the final
results/anchor_scan.json. A checkpoint is accepted only when all scientific
source and public-data SHA-256 values match the current candidate.
"""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SRC / "probes"))
import anchor_scan as A

FINAL = ROOT / "results/anchor_scan.json"
CKPT = ROOT / "results/anchor_scan_checkpoint.json"
CODE_FILES = [
    "src/anchor_checkpoint.py", "src/anchor_scan.py", "src/core.py",
    "src/probes/cmb.py", "src/probes/bao.py", "src/probes/sne.py",
    "data/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "data/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    "data/Pantheon+SH0ES.dat", "data/Pantheon+SH0ES_STAT+SYS.cov",
]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def signature():
    return {rel: sha(ROOT / rel) for rel in CODE_FILES}


def save_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(value, f, indent=2)
    os.replace(tmp, path)


def main(verify_only=False):
    sig = signature()
    saved = json.loads(CKPT.read_text()) if CKPT.exists() else None
    if saved and saved.get("code_signature") != sig:
        raise RuntimeError("checkpoint code/data signature mismatch; refusing mixed-version resume")
    if verify_only:
        print(json.dumps({"status": "CHECKPOINT_DRIVER_VERIFIED",
                          "checkpoint_present": saved is not None,
                          "checkpoint_complete": bool(saved and saved.get("complete")),
                          "code_data_files": len(sig)}, indent=2))
        return 0

    A.bao.load_desi()
    A.sne.load_pantheon()
    out = saved["out"] if saved else {}
    done = set(saved.get("done", [])) if saved else set()
    complete = False

    def checkpoint():
        # The cache is only an accelerator. Clear it after each durable row so
        # hundreds of independent fits cannot exhaust memory; values are not changed.
        if hasattr(A.cmb, "_Y_CACHE"):
            A.cmb._Y_CACHE.clear()
        save_json(CKPT, {"complete": complete, "done": sorted(done), "out": out,
                         "code_signature": sig,
                         "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat()})

    def finish_section(name):
        done.add(name)
        checkpoint()
        print("checkpoint section=" + name, flush=True)

    if "refs" not in done:
        chi2_ref, _ = A.cmb.fit_lcdm_reference()
        refs = {"bao": A.bao.fit_lcdm()[2],
                "sne": A.sne.run(0.01)["full"]["chi2_min"], "cmb": chi2_ref}
        out["refs"] = refs
        finish_section("refs")
    refs = out["refs"]
    print("refs:", refs, flush=True)

    def row1d(OL, r_map="quartic"):
        cb = A.bao_chi2(6.1, OL)
        cs = A.sne_chi2(6.1, OL)
        cc = A.cmb.fit_bimetric(6.1, OL, r_map=r_map)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        return [OL] + vals + [None if None in vals else sum(vals)]

    def point(th, OL):
        cb = A.bao_chi2(th, OL)
        cs = A.sne_chi2(th, OL)
        cc = A.cmb.fit_bimetric(th, OL)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        joint = None if None in vals else sum(vals)
        return [th, OL] + vals + [joint]

    def run_rows(section, items, formatter):
        if section in done:
            print("resume complete section=" + section, flush=True)
            return
        rows = out.setdefault(section, [])
        for key, compute in items:
            marker = section + ":" + key
            if marker in done:
                continue
            row = compute()
            rows.append(row)
            done.add(marker)
            checkpoint()
            print(formatter(row), flush=True)
        finish_section(section)

    run_rows("profile_1D",
             [(str(ol), lambda ol=ol: row1d(ol)) for ol in A.OLs_1D],
             lambda r: "1D OL=%.3f joint=%.3f" % (r[0], r[4]))
    run_rows("profile_1D_lambda1",
             [(str(ol), lambda ol=ol: row1d(ol, 1.0)) for ol in A.OLs_1D],
             lambda r: "1D-l1 OL=%.3f joint=%.3f" % (r[0], r[4]))
    run_rows("profile_1D_fine",
             [(str(ol), lambda ol=ol: row1d(ol)) for ol in A.OLs_FINE],
             lambda r: "fine OL=%.3f joint=%.3f" % (r[0], r[4]))

    def theta_fine(th):
        cb = A.bao_chi2(th, 0.699)
        cs = A.sne_chi2(th, 0.699)
        cc = A.cmb.fit_bimetric(th, 0.699)[0]
        vals = [x - r for x, r in zip((cb, cs, cc),
                                      (refs["bao"], refs["sne"], refs["cmb"]))]
        return [th] + vals + [sum(vals)]

    run_rows("theta_fine",
             [(str(th), lambda th=th: theta_fine(th)) for th in A.THETAS_FINE],
             lambda r: "theta_fine th=%.1f joint=%.3f" % (r[0], r[4]))
    run_rows("grid_2D",
             [(f"{th}:{ol}", lambda th=th, ol=ol: point(th, ol))
              for th in A.thetas_2D for ol in A.OLs_2D],
             lambda r: "2D th=%.1f OL=%.3f joint=%s" %
                       (r[0], r[1], "%.3f" % r[5] if r[5] is not None else "None"))
    run_rows("grid_2D_ext",
             [(f"{th}:{ol}", lambda th=th, ol=ol: point(th, ol))
              for th in A.THETAS_EXT for ol in A.OLs_EXT],
             lambda r: "ext th=%.1f OL=%.3f joint=%s" %
                       (r[0], r[1], "%.3f" % r[5] if r[5] is not None else "None"))

    if "lcdm_shared" not in done:
        om_s, joint_s, cost = A.lcdm_shared_cost()
        out["lcdm_shared"] = {"Om": om_s, "joint": joint_s, "cost": cost,
                              "chi2_sum_min": joint_s - cost}
        finish_section("lcdm_shared")

    def continuous_theta(th):
        raw = A.theta_profile(thetas=(th,), refs=refs)[0]
        return [th, round(raw[1], 4), round(raw[2], 2)]

    if "theta_profile_continuous" not in done:
        block = out.setdefault("theta_profile_continuous", {
            "note": "Omega_Lambda minimised continuously at each theta", "rows": []})
        for th in A.theta_profile.__defaults__[0]:
            marker = "theta_profile_continuous:" + str(th)
            if marker in done:
                continue
            block["rows"].append(continuous_theta(th))
            done.add(marker)
            checkpoint()
            print("theta_profile th=%.1f OL*=%.4f joint=%.2f" %
                  tuple(block["rows"][-1]), flush=True)
        finish_section("theta_profile_continuous")

    if "profile_lambda1_continuous" not in done:
        l1 = A.lambda1_profile()
        out["profile_lambda1_continuous"] = {
            "OL_star": round(l1["OL_star"], 4), "joint_min": round(l1["joint_min"], 2),
            "note": "continuous Omega_Lambda profile at theta=6.1, lambda=1"}
        finish_section("profile_lambda1_continuous")

    save_json(FINAL, out)
    done.add("anchor_scan.json")
    complete = True
    checkpoint()
    print("ANCHOR SCAN DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(verify_only="--verify-only" in sys.argv))
