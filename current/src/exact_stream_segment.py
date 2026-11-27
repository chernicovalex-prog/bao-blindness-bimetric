#!/usr/bin/env python3
"""Run one exact-stream segment of profile_calibration.py with durable resume."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
CODE_FILES = [
    "src/exact_stream_segment.py", "src/profile_calibration.py", "src/core.py",
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


def save_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(value, f, indent=2)
    os.replace(tmp, path)


def signature():
    return {rel: sha(ROOT / rel) for rel in CODE_FILES}


def rng_state(seed, count, zcut):
    rng = np.random.default_rng(seed)
    for _ in range(count):
        rng.standard_normal(13)
        rng.standard_normal(zcut)
        rng.standard_normal(3)
    return rng.bit_generator.state


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--grid", choices=["fine", "offset"], default="fine")
    args = parser.parse_args()
    if not (0 <= args.start < args.end):
        raise SystemExit("require 0 <= start < end")

    raw = [line.split() for line in (ROOT / "data/Pantheon+SH0ES.dat").read_text().splitlines()
           if line.strip() and not line.startswith("#")]
    zcut = sum(float(row[2]) > 0.01 for row in raw)
    ckpt = ROOT / f"results/profile_calibration_checkpoint_{args.tag}.json"
    provenance = ROOT / f"results/exact_stream_{args.tag}_provenance.json"
    final = ROOT / f"results/profile_calibration_{args.tag}.json"
    sig = signature()
    if provenance.exists():
        old = json.loads(provenance.read_text())
        if old["code_data_signature"] != sig:
            raise SystemExit("segment code/data signature mismatch; refusing mixed-version resume")
        if old["start"] != args.start or old["end"] != args.end or old["seed"] != args.seed:
            raise SystemExit("segment identity mismatch")
    elif ckpt.exists():
        raise SystemExit("checkpoint exists without provenance; refusing unknown resume")

    if not ckpt.exists():
        state = {"done": args.start, "prof_argmin": [], "cov1": [], "cov384": [],
                 "dchi_true": [], "rng_state": rng_state(args.seed, args.start, zcut)}
        save_json(ckpt, state)
        save_json(provenance, {
            "schema": "exact_stream_segment_v1", "start": args.start, "end": args.end,
            "seed": args.seed, "tag": args.tag, "grid_mode": args.grid,
            "zcut_count": zcut, "code_data_signature": sig,
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "checkpoint": ckpt.name, "random_stream": "single contiguous stream"})

    command = [sys.executable, "profile_calibration.py", str(args.end), str(args.seed),
               args.tag, args.grid]
    print("EXACT_SEGMENT_COMMAND=" + " ".join(command), flush=True)
    result = subprocess.run(command, cwd=SRC)
    if result.returncode != 0:
        return result.returncode
    saved = json.loads(ckpt.read_text())
    expected = args.end - args.start
    if saved["done"] != args.end or len(saved["dchi_true"]) != expected:
        print(json.dumps({"status": "SEGMENT_INCOMPLETE", "done": saved["done"],
                          "expected": args.end, "records": len(saved["dchi_true"])}))
        return 3
    output = json.loads(final.read_text())
    output.update({"N": expected, "seed": args.seed, "seeds": [args.seed],
                   "segment_start": args.start, "segment_end": args.end,
                   "exact_stream": True, "exact_stream_tag": args.tag})
    save_json(final, output)
    print(json.dumps({"status": "EXACT_SEGMENT_COMPLETE", "tag": args.tag,
                      "start": args.start, "end": args.end, "records": expected}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
