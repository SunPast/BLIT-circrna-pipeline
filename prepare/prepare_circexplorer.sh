#!/usr/bin/env bash
set -euo pipefail

REF_BASE=/data/home/dingjia/pipeline/
ANN_REF=${REF_BASE}/hg38_ref_all.txt

echo "===> Circexplorer2 annotation"

if [ ! -f "${ANN_REF}" ]; then
  fetch_ucsc.py hg38 > "${ANN_REF}"
fi

echo "Circexplorer2 annotation ready."
