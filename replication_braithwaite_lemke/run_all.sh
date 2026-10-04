#!/usr/bin/env bash
# Rebuild all data and results for the Braithwaite & Lemke (2011) replication.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output/logs

python3 src/validate_heckprob.py
python3 src/build_mids.py
python3 src/build_dyadyears.py
for s in main corrected jd_correct jd_stata_all fatal_strict polity2 contig_24 contig_land terr_any_participant; do
  echo "estimating spec: $s"
  python3 src/estimate.py --spec "$s" > "output/logs/estimate_$s.log" 2>&1
done
python3 src/summarize.py
python3 src/make_figure.py
