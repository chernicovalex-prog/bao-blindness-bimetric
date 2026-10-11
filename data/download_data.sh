#!/usr/bin/env bash
# Download the public datasets used by the probes.
# Pantheon+ (Brout et al. 2022; Scolnic et al. 2022), GitHub DataRelease.
set -euo pipefail
cd "$(dirname "$0")"

BASE="https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_data/4_distances_covariance"

echo "Downloading Pantheon+SH0ES.dat ..."
curl -fL "$BASE/Pantheon%2BSH0ES.dat" -o "Pantheon+SH0ES.dat"

echo "Downloading Pantheon+SH0ES_STAT+SYS.cov (~23 MB) ..."
curl -fL "$BASE/Pantheon%2BSH0ES_STAT%2BSYS.cov" -o "Pantheon+SH0ES_STAT+SYS.cov"

echo "Done. Files:"
ls -lh "Pantheon+SH0ES.dat" "Pantheon+SH0ES_STAT+SYS.cov"
