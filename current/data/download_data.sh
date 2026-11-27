#!/usr/bin/env bash
# Download the official public Pantheon+ covariance used by the current pipeline.
# The curated Pantheon+SH0ES.dat used by this candidate is preserved in current/data/.
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -f "Pantheon+SH0ES.dat" ]]; then
  echo "The reviewed curated Pantheon+SH0ES.dat is missing from current/data/" >&2
  exit 1
fi

URL="https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_data/4_distances_covariance/Pantheon%2BSH0ES_STAT%2BSYS.cov"
echo "Downloading the official Pantheon+ STAT+SYS covariance ..."
curl --retry 3 --retry-delay 5 -fL "$URL" -o "Pantheon+SH0ES_STAT+SYS.cov"
ls -lh "Pantheon+SH0ES.dat" "Pantheon+SH0ES_STAT+SYS.cov"
