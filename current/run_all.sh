#!/usr/bin/env bash
# Current staged reproduction route for the corrected candidate.
# Independent gates, raw artifacts, exact-stream MC and final QA; no monolithic
# production claim, no tied-r_d/H&M reconstruction, no MCMC.
set -euo pipefail
cd "$(dirname "$0")/src"

# 1. Independent and regression controls.
python3 cmb_spec_gates.py
python3 independent_bg.py
python3 scalarsec.py
python3 scalar_mfp_scan.py
python3 selfsim.py
python3 probes/sne.py
python3 probes/bao.py

# 2. CMB, grids, lambda sweeps and tables.
python3 probes/cmb.py
python3 scan_cmb_grids.py
python3 _lam.py
python3 lambda_theta_sweep.py
python3 make_tables.py

# 3. Anchor tree, scalar-at-minima and extra numbers.
python3 anchor_checkpoint.py
python3 scalar_at_minima.py
python3 extra_numbers.py

# 4. Exact-stream equivalence control; full MC is blocked if this fails.
python3 exact_stream_segment.py --start 0 --end 10 --tag exactctrl0b
python3 exact_stream_segment.py --start 10 --end 20 --tag exactctrl1b
python3 merge_exact_stream.py --mode control

# 5. Contiguous exact-stream fine segments, offset and CMB bias.
python3 exact_stream_segment.py --start 0 --end 500 --tag exactF0
python3 exact_stream_segment.py --start 500 --end 1000 --tag exactF1
python3 exact_stream_segment.py --start 1000 --end 1500 --tag exactF2
python3 exact_stream_segment.py --start 1500 --end 2000 --tag exactF3
python3 profile_calibration.py 1000 20261119 offset offset
python3 cmb_linbias_check.py 25
python3 merge_exact_stream.py --mode full

# 6. Final assertions and figures.
python3 checks.py
python3 independent_bg.py
python3 make_figures.py

echo "STAGED REPRODUCTION COMPLETE -> results/ and figures/"
