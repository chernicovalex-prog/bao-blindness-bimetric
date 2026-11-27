#!/usr/bin/env bash
# Download the official public Pantheon+ files used by the current pipeline.
set -euo pipefail
cd "$(dirname "$0")"

BASE="https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR"
echo "Downloading the official Pantheon+ distance table ..."
curl --retry 3 --retry-delay 5 -fL "$BASE/Pantheon%2BSH0ES.dat" -o "Pantheon+SH0ES.dat"
echo "Downloading the official Pantheon+ STAT+SYS covariance ..."
curl --retry 3 --retry-delay 5 -fL "$BASE/Pantheon%2BSH0ES_STAT%2BSYS.cov" -o "Pantheon+SH0ES_STAT+SYS.cov"

# Pin the reviewed covariance bytes. Normalize only the official table's
# column-name line to a comment; the 1701 data rows remain unchanged.
python3 - <<'PY'
import hashlib
from pathlib import Path
p = Path('Pantheon+SH0ES_STAT+SYS.cov')
got = hashlib.sha256(p.read_bytes()).hexdigest()
want = 'abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc'
if got != want:
    raise SystemExit(f'Pantheon+ covariance SHA mismatch: {got}')
table = Path('Pantheon+SH0ES.dat')
lines = table.read_text().splitlines()
for i, line in enumerate(lines):
    if line.strip() and not line.startswith('#'):
        if not line.startswith('CID '):
            raise SystemExit('Unexpected first non-comment Pantheon+ line')
        lines[i] = '# ' + line
        break
rows = [line.split() for line in lines if line.strip() and not line.startswith('#')]
if len(rows) != 1701 or any(len(row) < 10 for row in rows):
    raise SystemExit('Unexpected Pantheon+ distance-table structure')
table.write_text('\n'.join(lines) + '\n')
print('Pantheon+ inputs: 1701 rows; reviewed covariance SHA-256 PASS')
PY
ls -lh "Pantheon+SH0ES.dat" "Pantheon+SH0ES_STAT+SYS.cov"
