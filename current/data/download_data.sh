#!/usr/bin/env bash
# Download the official public Pantheon+ files used by the current pipeline.
set -euo pipefail
cd "$(dirname "$0")"

BASE="https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_data/4_distances_covariance"
echo "Downloading the official Pantheon+ distance table ..."
curl --retry 3 --retry-delay 5 -fL "$BASE/Pantheon%2BSH0ES.dat" -o "Pantheon+SH0ES.dat"
echo "Downloading the official Pantheon+ STAT+SYS covariance ..."
curl --retry 3 --retry-delay 5 -fL "$BASE/Pantheon%2BSH0ES_STAT%2BSYS.cov" -o "Pantheon+SH0ES_STAT+SYS.cov"

# Pin the reviewed covariance bytes; the numerical gates in src/checks.py
# verify the fields used from the official distance table and all quoted SNe results.
python3 - <<'PY'
import hashlib
from pathlib import Path
p = Path('Pantheon+SH0ES_STAT+SYS.cov')
got = hashlib.sha256(p.read_bytes()).hexdigest()
want = 'abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc'
if got != want:
    raise SystemExit(f'Pantheon+ covariance SHA mismatch: {got}')
rows = [line.split() for line in Path('Pantheon+SH0ES.dat').read_text().splitlines()
        if line.strip() and not line.startswith('#')]
if len(rows) != 1701 or any(len(row) < 10 for row in rows):
    raise SystemExit('Unexpected Pantheon+ distance-table structure')
print('Pantheon+ inputs: 1701 rows; reviewed covariance SHA-256 PASS')
PY
ls -lh "Pantheon+SH0ES.dat" "Pantheon+SH0ES_STAT+SYS.cov"
